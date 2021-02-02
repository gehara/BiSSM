preprocess<-function(GIS_data, hyp.matrix){
  if(is.null(hyp.matrix)){ 
    GEO_distance <- distm(cbind(unique(GIS_data[,1:2])))
  } else {
    GEO_distance <- as.matrix(hyp.matrix)
    colnames(GEO_distance) <- NULL
    }
    n_local <- length(unique(GIS_data[,4]))
    populations <- unique(GIS_data$species)
    n_pops <- length(populations)

    #################

    pops <- list()
    for(i in 1:n_pops){
      pops[[i]]<-grep(populations[i],GIS_data$species)
      }

      combi_pops <- t(combn(1:n_pops,2))

  local_info <- NULL
  for(i in 1:nrow(GIS_data)){
    local_info <- c(local_info,rep(GIS_data[i,]$site,GIS_data[i,]$samp))
  }
  
  n_samples<-NULL
  for(i in 1:n_local) n_samples <-c(n_samples, length(which(local_info==i)))
  
  out <- list(GEO_distance, n_local, populations, n_pops, combi_pops, local_info,n_samples, pops)
  names(out) <- c("GEO_distance","n_localities","Populations","N_populations","comb_pops","samples_local","n_samples","local_assign2pop") 
 return(out) 
}

