sim.obs <- function (sim, obs) 
{
  mylabels <- colnames(sim)
  for (i in 1:ncol(sim)) {
    hist(sim[, i], breaks = 20, xlab = mylabels[i], main = "", 
         xlim = c(min(c(na.omit(sim[, i]), obs[,i])), max(c(na.omit(sim[,i]), obs[,i]))))
    
    for(j in 1:nrow(obs)){
    abline(v = obs[j,i], col = 2)
    }
    
    }
}
