
sim.IBD <- function(GIS_data,dispersal=c(100,1000),Ne=c(1000,10000), hyp.matrix=NULL, mig.prior = NULL, pop.assign, neighbours,
                  mutation_rate=c(1e-11,1e-10),reduced_mig=c(1,1), bp=77,nsim.blocks, nloci, block.size){

  # This will preprocess the input. Generate the geographic distance matrix from GPS data, plus other objects.
  proc_data <- preprocess(GIS_data, hyp.matrix)
  
  ## gets the msABC binary from PipeMaster
  msABC <- PipeMaster:::get.msABC()
 
  ## this is the simulation loop
  IBD <- NULL
  thou <- 0
  for(z in 1:nsim.blocks){
    
    ibd <- NULL
  
    TIM <- system.time(
    
    for(j in 1:block.size){
      
      ## make the stepping stone model with a fixed number of neighbors.
      ss.model <- make.ssm(proc_data, neighbours)
    
    if(is.null(mig.prior)){  
      # sample dispersal rate per generation
    
      disp <- runif(1,dispersal[1],dispersal[2])
    
      #
      redu_mig_factor <- runif(1,reduced_mig[1],reduced_mig[2])
    
      # generates the migration matrix
      migmatrix <- (disp/as.numeric(ss.model$locality_migrations[,2]))*4
      
      # generate some noise 
      amount <- runif(1,0,1)
      noize <- migmatrix*amount
      migmatrix <- migmatrix + (noize*replicate(length(migmatrix),sample(c(1,-1),size=1)))
    
      for(i in 1:length(migmatrix)) if(migmatrix[i]<=0.1) migmatrix[i] <- migmatrix[i] + 0.1
    
    } else {
      
      migmatrix <- rep(runif(1, mig.prior[1], mig.prior[2]),nrow(ss.model$locality_migrations))*4
      
      disp <- (mean(migmatrix)/4) * mean(as.numeric(ss.model$locality_migrations[,2]))
      amount <- runif(1,0,1)
      noize <- migmatrix*amount
      migmatrix <- migmatrix + (noize*replicate(length(migmatrix),sample(c(1,-1),size=1)))
      
      for(i in 1:length(migmatrix)) if(migmatrix[i]<=0.1) migmatrix[i] <- migmatrix[i] + 0.1
      
        redu_mig_factor <- runif(1,reduced_mig[1],reduced_mig[2])
      
    }
    
    # calculate mean number of migrants
    mean.Nm <- NULL
    for(w in 1:length(ss.model$Nm)){
      mean.Nm <- c(mean.Nm, sum(migmatrix[ss.model$Nm[[w]]]))
    }
    mean.Nm <- mean(mean.Nm)/4
    mig_mean <- mean(migmatrix)/4
    
    ## generates the ms string (migrations)
    lm <- paste(ss.model$locality_migrations[,1], round(migmatrix, 4))
    mig_string <- paste(paste("-m",lm), collapse = " ")

    ## sample Ne
    Ne_pop <- runif(1,Ne[1],Ne[2])
    ## sample mutation rate per site per generation 
    mi <- runif(1,mutation_rate[1],mutation_rate[2])
    
    # theta scale
    theta <- round(4*Ne_pop*mi*bp,4)
    
    # sample growth rate
    G <- rtnorm(1, 0, 3, 0)
    
    # ms string
    ms_string <- paste("-t",theta,"-I", proc_data$n_localities, paste(proc_data$n_samples, collapse=" "), collapse = " ")
    ms_string <- paste(ms_string, mig_string)
    ms_string <- paste(ms_string, mig_string,"-G", G)
  
    ## simulate data
    X <- read.table(text = system(paste(msABC, sum(proc_data$n_samples), nloci, ms_string), intern = T), sep="\t", header=T)
    
    ## select summary stats
    X <- cbind(X[, grep("segs",colnames(X))], X[, grep("pi",colnames(X))], X[,grep("fst", colnames(X))],X[,grep("perc", colnames(X))], X[match("s_tajimasD", colnames(X))])
    
    ## mean of summary stats
    X <- colMeans(X, na.rm=T)
    
    # Resuting simulation
    DNA_dist <- X
   
    ## model parameters
    param <- c(Ne_pop, disp, amount, migmatrix/4)
    names(param) <- c("Ne","dispersal","noise", paste("mig", ss.model$locality_migrations[,1], sep="_Nm_"))
    
    ibd <- rbind(ibd, c(param, DNA_dist))
    print(j)
    }
  )[3]

  thou <- thou+block.size
  Total.time<-round((((TIM*nsim.blocks)-(TIM*z))/60)/60,3)
  cat(paste("PipeMaster:: ",thou," (",round(block.size/(TIM),1)," sims/sec) | ",Total.time," hours remaining",sep=""),"\n")
  IBD <- rbind(IBD,ibd)
  }
  return(IBD)
}
