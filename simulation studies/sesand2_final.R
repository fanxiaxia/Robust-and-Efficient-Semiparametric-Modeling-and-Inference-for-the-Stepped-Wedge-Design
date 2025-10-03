se.sand2 <- function(dataset, deltahat,gammahat, W.list, diffsize, adj.var.prop)
{
  dataset <- merge(dataset, gammahat, by = c("cluster","time","ind"), all = TRUE)
  dataset <- dataset %>% mutate(y0 = y - x*deltahat - gamma)
  dataset <- dataset %>% arrange(cluster,time,ind)
  
  datacat <- propx.cal(dataset,ex.type = ex.type.set,adj.var = adj.var.prop)
    
    num <- denom <-  0  
    for(i in 1:Ncluster)
    {
      cluster.index <- i
      data.i <- datacat %>% filter(cluster == cluster.index)
      
      num.i <- t(data.i$L) %*% W.list[[cluster.index]] %*% matrix(data.i$y0,ncol = 1) %*% 
        matrix(data.i$y0,nrow = 1) %*% W.list[[cluster.index]] %*% data.i$L
      denom.i <- t(data.i$L) %*% W.list[[cluster.index]] %*% data.i$x
      
      num <- num + num.i
      denom <- denom + denom.i
    }
  
  A <- denom/Ncluster
  B <- num/Ncluster
  
  se <-  sqrt((B/(A^2))/Ncluster)
  return(se)
}