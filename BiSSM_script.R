# load required packages
library(PipeMaster)
#library(gdm)
library(reshape)
library(geosphere)
#library(adegenet)
#library(FNN)
library(keras)

#save working directory
path <- getwd()

## load functions
source("IBD.R")
source("preprocess.R")
source("make.ssm.R")
source("plotPCS2.R")

### read inputs
samp <- read.table("North.txt", header=T, stringsAsFactors = F)
samp[1,6] <- "B"
pop.assign <- read.table("N_pop_assign.txt", header=T, stringsAsFactors = F)
nloci <- as.numeric(strsplit(list.files(pattern = "N_observed"),"_")[[1]][3])

## you can use the loop below to get the average number of base pairs across n loci (fasta alignments)
# seqs <- list.files("./N_fasta", pattern = ".fas")
# for(i in 1:length(seqs)) bp <- c(bp, ncol(read.dna(paste("./N_fasta/",seqs[i],sep=""), format = "fasta")))
# bp <- round(mean(bp))
bp <- 89


# define number of neighbors based on the total number of localities.
if(length(unique(samp$site)) >= 3) neighbours <- 2
if(length(unique(samp$site)) == 2) neighbours <- 1

## this specifies the number of simulations you want to run for the training data
# the total number of simulations is equal to block.size x nsim.blocks
# block.size is the number of simulations run before writing to file. 
# R gets slow when storing a lot of memory. So don't go too high on block.size. The highest I would go is 1000.
# You can increase the number of simulations by increasing nsim.blocks.
#
block.size = 100
nsim.blocks = 1

# same as above but for the testing data
test.block.size = 10
test.nsim.blocks = 1


# simulate training data
#' @param GIS_data data frame with 5 columns. long; lat; number of samples for the locality; site number; species 
#' @param block.size size of the simulation block. Total number of simulations = block.size * nsim.blocks
#' @param nsim.blocks total number of simulation blocks.
#' @param nloci number of loci to simulate
#' @param dispersal min - max dispersal rate in meters per generation
#' @param neighbours number of neighbors of each locality in the bidimentional stepping stone model.
#' @param hyp.matrix distance matrix hypothesis, different than linear geographic distance
#' @param Ne min - max Ne
#' @param reduced_mig multiplier to reduce migration between species. It should be 1 if a single species is analyzed.
#' @param pop.assign data frame with two columns. sample name; locality number.
#' @param mutation_rate min max mutation rate
#' @param bp average number of base pairs across loci 
IBD <- sim.IBD(GIS_data = samp[,2:6],
               block.size = block.size,
               nsim.blocks = nsim.blocks, nloci = nloci,
               dispersal = c(1000, 200000), neighbours = neighbours, hyp.matrix = NULL,
               Ne = c(20000, 400000), reduced_mig = c(1, 1),  pop.assign = pop.assign,
               mutation_rate = c(2.5e-9, 2.5e-9),
               bp = bp)

write.table(IBD, "IBD_N", col.names = T, row.names = F, quote = F, sep="\t")

Island <- sim.IBD(GIS_data = samp[,2:6],
                  block.size = block.size,
                  nsim.blocks = nsim.blocks, nloci = nloci, mig.prior = c(0.25, 5),
                  dispersal = NULL, neighbours = neighbours, hyp.matrix = NULL,
                  Ne = c(20000, 400000), reduced_mig = c(1, 1),  pop.assign = pop.assign,
                  mutation_rate = c(2.5e-9, 2.5e-9),
                  bp = bp)

write.table(Island, "Island_N", col.names = T, row.names = F, quote = F, sep="\t")

## simulate test data
IBD <- sim.IBD(GIS_data = samp[,2:6],
               block.size = test.block.size,
               nsim.blocks =  test.nsim.blocks, nloci = nloci,
               dispersal = c(1000, 200000), neighbours = neighbours, hyp.matrix = NULL,
               Ne = c(20000, 400000), reduced_mig = c(1, 1),  pop.assign = pop.assign,
               mutation_rate = c(2.5e-9, 2.5e-9),
               bp = bp)

write.table(IBD, "test_IBD_N", col.names = T, row.names = F, quote = F, sep="\t")

Island <- sim.IBD(GIS_data = samp[,2:6],
               block.size = test.block.size,
               nsim.blocks = test.nsim.blocks, nloci = nloci, mig.prior = c(0.25, 5),
               dispersal = NULL, neighbours = neighbours, hyp.matrix = NULL,
               Ne = c(20000, 400000), reduced_mig = c(1, 1),  pop.assign = pop.assign,
               mutation_rate = c(2.5e-9, 2.5e-9),
               bp = bp)

write.table(Island, "test_Island_N", col.names = T, row.names = F, quote = F, sep="\t")

## read simulations
IBD <- read.table("IBD_N", sep="\t", header = T)
Island <- read.table("Island_N", sep="\t", header = T)

## read observed data
obs <- read.table(list.files(pattern = "_observed"), header=T)

## create index
index <- c(rep(0, nrow(IBD)), rep(1, nrow(Island)))
## bind simulations in a single data frame
sims <- rbind(IBD, Island)
sims <- sims[,colnames(sims)%in% colnames(obs)]


## Pdf of prior predictive checks, PCA and model fit.
pdf("PCA.pdf")
  plotPCS2(models = sims, index, obs, subsample = 1)

  plot(gfit(target = obs, 
          sumstat = Island[,colnames(Island)%in% colnames(obs)],
          nb.replicate = 50),
          main = "Island", breaks = 20)

  plot(gfit(target = obs, 
          sumstat = IBD[,colnames(IBD)%in% colnames(obs)],
          nb.replicate = 50),
     main = "IBD", breaks = 20)
dev.off()



###########################
######################
############### machine learning IBD

epochs <- 10 ## how long do you want to run the training

  ## This is the training data  
  sims <- IBD[,-ncol(IBD)] ## IBD dataframe without the last column, TajD.
  sims <- sims[, colnames(sims)%in% colnames(obs)]
  ## normalize the outcome by the mean
  means <- colMeans(IBD)[1:3]
  sims <- cbind(apply(IBD, 2, function(x) x/mean(x))[,1:3], sims)
  
  ## This is the testing data
  temp.test <- read.table("test_IBD_N", sep="\t", header = T)[, -ncol(IBD)]
  temp.test <- temp.test[,colnames(temp.test)%in% colnames(obs)]
  x <- read.table("test_IBD_N", sep="\t", header = T)[, 1:3]
  for(i in 1:3) x[,i] <- x[,i]/means[i]
  temp.test <- cbind(x,  temp.test)
  
  Res <- NULL
  for(i in 1:3){
    
    test <- list(temp.test[, 4:ncol(sims)], temp.test[, i])
    train <- list(sims[, 4:ncol(sims)], sims[, i])
    names(train) <- c("x","y")
    names(test) <- c("x","y")
    
    proc.data <- list(train, test)
    names(proc.data) <- c("train", "test")
    
    
    c(train_images, train_labels) %<-% proc.data$train
    c(test_images, test_labels) %<-% proc.data$test
    
    train_images <- scale(train_images)
   
    # Use means and standard deviations from training set to normalize test set
    col_means_train <- attr(train_images, "scaled:center")
    col_stddevs_train <- attr(train_images, "scaled:scale")
    test_images <- scale(test_images, center = col_means_train, scale = col_stddevs_train)
    
    model <- keras_model_sequential() %>%
      layer_dense(units = 32, activation = "relu",
                  input_shape = dim(train_images)[2]) %>%
      layer_dense(units = 32, activation = "relu") %>%
      layer_dense(units = 1, activation = "relu")
    
    model %>% compile(
      loss = "mae",
      optimizer = optimizer_rmsprop(),
      metrics = list("mean_absolute_error"))
    
    
    model %>% fit(
      train_images,
      train_labels,
      epochs = epochs,
      #batch_size = 10000,
      validation_split = 0.1)
    
    mae %<-% (model %>% evaluate(test_images, test_labels, verbose = 0))
    
    test_pred <- model %>% predict(test_images)
    cor.coef <- cor.test(test_pred, test_labels)$estimate
    
    pdf(paste(colnames(temp.test)[i],"_IBD.pdf",sep=""))
    plot(test_pred ~ test_labels, xlab = "true", ylab="estimated", main = colnames(temp.test)[i])
    lines(c(0:1000000),(0:1000000), col=2)
    dev.off()
    
    obs <- read.table(list.files(pattern ="observed"), header=T)[,-38]
    
    TAB2 <- scale(obs, center = col_means_train, scale = col_stddevs_train)
    predictions <- model %>% predict(TAB2)
    predictions <- predictions*means[i]
    
    z <- c(predictions, cor.coef, mae[1])
    names(z)<- c("estimate", "r2", "MAE")
    Res <- rbind(Res, z)
    
    
  }
    
  rownames(Res) <- colnames(temp.test)[1:3]
  
  write.table(Res, "IBD_result.txt")


  
  ###########################
  ######################
  ############### machine learning Island
  
  # training data
  sims <- Island[,-ncol(Island)] ## Island dataframe without the last column, TajD.
  sims <- sims[, colnames(sims)%in% colnames(obs)]
  means <- colMeans(Island)[1:3]
  sims <- cbind(apply(Island, 2, function(x) x/mean(x))[,1:3], sims)
  
  # testing data
  temp.test <- read.table("test_Island_N", sep="\t", header = T)[,-ncol(Island)]
  temp.test <- temp.test[,colnames(temp.test)%in% colnames(obs)]
  x <- read.table("test_Island_N", sep="\t", header = T)[, 1:3]
  for(i in 1:3) x[,i] <- x[,i]/means[i]
  temp.test <- cbind(x,  temp.test)
  
  Res <- NULL
  for(i in 1:3){
    
    test <- list(temp.test[, 4:ncol(sims)], temp.test[, i])
    train <- list(sims[, 4:ncol(sims)], sims[, i])
    names(train) <- c("x","y")
    names(test) <- c("x","y")
    
    proc.data <- list(train, test)
    names(proc.data) <- c("train", "test")
    
    
    c(train_images, train_labels) %<-% proc.data$train
    c(test_images, test_labels) %<-% proc.data$test
    
    train_images <- scale(train_images)
    
    # Use means and standard deviations from training set to normalize test set
    col_means_train <- attr(train_images, "scaled:center")
    col_stddevs_train <- attr(train_images, "scaled:scale")
    test_images <- scale(test_images, center = col_means_train, scale = col_stddevs_train)
    
    model <- keras_model_sequential() %>%
      layer_dense(units = 32, activation = "relu",
                  input_shape = dim(train_images)[2]) %>%
      layer_dense(units = 32, activation = "relu") %>%
      layer_dense(units = 1, activation = "relu")
    
    model %>% compile(
      loss = "mae",
      optimizer = optimizer_rmsprop(),
      metrics = list("mean_absolute_error"))
    
    
    model %>% fit(
      train_images,
      train_labels,
      epochs = epochs,
      #batch_size = 10000,
      validation_split = 0.1)
    
    mae %<-% (model %>% evaluate(test_images, test_labels, verbose = 0))
    
    test_pred <- model %>% predict(test_images)
    cor.coef <- cor.test(test_pred, test_labels)$estimate
    
    pdf(paste(colnames(temp.test)[i],"_Island.pdf",sep=""))
    plot(test_pred ~ test_labels, xlab = "true", ylab="estimated", main = colnames(temp.test)[i])
    lines(c(0:1000000),(0:1000000), col=2)
    dev.off()
    
    obs <- read.table(list.files(pattern ="observed"), header=T)[,-38]
    
    TAB2 <- scale(obs, center = col_means_train, scale = col_stddevs_train)
    predictions <- model %>% predict(TAB2)
    predictions <- predictions*means[i]
    
    z <- c(predictions, cor.coef, mae[1])
    names(z)<- c("estimate", "r2", "MAE")
    Res <- rbind(Res, z)
    
    
  }
  
  rownames(Res) <- colnames(temp.test)[1:3]
  
  write.table(Res, "Island_result.txt")
  