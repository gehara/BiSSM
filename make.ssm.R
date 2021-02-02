
make.ssm <- function(proc_data, neighbours){

  GEO_distance <- proc_data$GEO_distance
  n_local <- proc_data$n_localities
  
  #### make the steping stone model
  if(length(neighbours)>1){
    k <- sample(neighbours,1)
  } else {k <- neighbours}
  
  stepstone <- FNN::get.knn(GEO_distance, k=k)
  #stepstone <- as.vector(stepstone$nn.index)
  
  temp <- NULL
  for(i in 1: nrow(stepstone$nn.index))
  {
    temp <- rbind(temp, cbind(i,stepstone$nn.index[i,]))
  }
  stepstone$nn.index <- rbind(temp,temp[,2:1])
  stepstone$nn.index <- unique(stepstone$nn.index)
  
  Nm <- list()
    for(i in 1:proc_data$n_localities){
    Nm[[i]] <-  which(stepstone$nn.index[,1] %in% i)
    }
  
  SS_dist <- NULL
  locality_migs <- NULL
  dis <- GEO_distance[lower.tri(GEO_distance)]
  temp.comb <- t(combn(1:n_local,2))
  temp.comb <- cbind(temp.comb, temp.comb[,1])
  if(nrow(temp.comb)>1)
    {
    temp.names <- apply(temp.comb[,1:2],1,paste,collapse=" ")
    temp.names <- c(temp.names, apply(temp.comb[,2:3],1,paste,collapse=" "))
  } else {
    temp.names <- paste(temp.comb[1,1:2],collapse=" ")
    temp.names <- c(temp.names, paste(temp.comb[1,2:3],collapse=" "))
  }
  
  dis <- c(dis, dis)
  names(dis) <- temp.names 

for(i in 1:nrow(stepstone$nn.index)){
  locality_migs <- rbind(locality_migs, paste(stepstone$nn.index[i,1:2],collapse = " "))
  #locality_migs <- rbind(locality_migs, paste(stepstone[i,2:3],collapse = " "))
  SS_dist <- rbind(SS_dist, dis[which(names(dis) == paste(stepstone$nn.index[i,1:2], collapse = " "))])
  #SS_dist <- rbind(SS_dist, dis[which(names(dis) == paste(stepstone[i,2:3], collapse = " "))])
}
  out <- list(Nm, stepstone$nn.index , cbind(locality_migs, round(SS_dist)), k)
  names(out) <- c("Nm","nn.index","locality_migrations","knn")
  return(out)
}
#############