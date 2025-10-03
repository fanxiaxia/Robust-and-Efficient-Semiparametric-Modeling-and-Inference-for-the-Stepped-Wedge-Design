se.perm3 <- function(dataset,deltahat,gammahat, W.list,diffsize,adj.var.prop,crosstime)
{
  dataset <- merge(dataset, gammahat, by = c("cluster","time","ind"), all = TRUE)
  dataset <- dataset %>% mutate(y0 = y - x*deltahat - gamma)
  dataset <- dataset %>% arrange(cluster,time,ind)
  
  #datacat <- propx.cal(dataset,ex.type = ex.type.set,adj.var = adj.var.prop)
  
  XLu <- data.frame(dataset$cluster,dataset$time, dataset$ind)
  colnames(XLu) <- c("cluster","time","ind")
  for(u in 1:Nperm)
  {
    L_name <- paste0("L_", u)
    X_name <- paste0("X_", u)
    
    crosstime.u <- sample(crosstime,size = Ncluster)
    crosscluster.u <- data.frame(1:Ncluster, crosstime.u)
    colnames(crosscluster.u) <- c("cluster","crosstime")
    dataset.u <- merge(dataset %>% dplyr::select(-x), crosscluster.u, by = "cluster", all = TRUE)
    dataset.u <- dataset.u %>%
      group_by(cluster) %>%
      mutate(x = ifelse(time < crosstime, 0, 1))
    
    datacat.u <- propx.cal(dataset.u,ex.type = ex.type.set,adj.var = adj.var.prop)
    
    selected_dataset.u <- datacat.u %>%
      dplyr::select(cluster,time,ind,L,x)
    
    XLu <- left_join(XLu,selected_dataset.u, by = c("cluster","time","ind"))
    XLu <- XLu %>% arrange(cluster,time,ind)
    XLu <- XLu %>% rename(!!L_name := L, !!X_name := x)
  }
  
  
  denom.vec <- matrix(NA, nrow = Ncluster, ncol = Nperm)
  for(i in 1:Ncluster)
  {
    cluster.index <- i
    data.i <- dataset %>% filter(cluster == cluster.index)
    XL.i <- XLu %>% filter(cluster == cluster.index)
    data.i <- left_join(data.i[,c("cluster","time","ind")],XL.i,by = c("cluster","time","ind"))
    L.i <- data.i %>% dplyr::select(starts_with("L_"))
    Xu.i <- data.i %>% dplyr::select(starts_with("X_"))
    denom.i <- diag(t(matrix(unlist(L.i),ncol = Nperm)) %*% W.list[[cluster.index]] %*% matrix(unlist(Xu.i),ncol = Nperm))
    denom.vec[i,] <- denom.i
  }
  
  denom.byperm <- apply(denom.vec,2,sum)
  
  term1 <- term2 <- 0
  for(i in 1:Ncluster)
  {
    cluster.index <- i
    data.i <- dataset %>% filter(cluster == cluster.index)
    XL.i <- XLu %>% filter(cluster == cluster.index)
    data.i <- left_join(data.i[,c("cluster","time","ind","y0")],XL.i,by = c("cluster","time","ind"))
    L.i <- matrix(unlist(data.i %>% dplyr::select(starts_with("L_"))),ncol = Nperm)
    
    num.i <- diag(t(L.i) %*% W.list[[cluster.index]] %*% matrix(data.i$y0,ncol = 1) %*% 
                    matrix(data.i$y0,nrow = 1) %*% W.list[[cluster.index]] %*% L.i)
    term1 <- term1 + mean(num.i/(denom.byperm^2))
  }
  
  for(k in 1:(Ncluster -1))
  {
    for(kprime in (k+1):Ncluster)
    {
      cluster.index.k <- k
      data.k <- dataset %>% filter(cluster == cluster.index.k)
      XL.k <- XLu %>% filter(cluster == cluster.index.k)
      data.k <- left_join(data.k[,c("cluster","time","ind","y0")],XL.k,by = c("cluster","time","ind"))
      L.jk <- matrix(unlist(data.k %>% dplyr::select(starts_with("L_"))),ncol = Nperm)
      
      cluster.index.kprime <- kprime
      data.kprime <- dataset %>% filter(cluster == cluster.index.kprime)
      XL.kprime <- XLu %>% filter(cluster == cluster.index.kprime)
      data.kprime <- left_join(data.kprime[,c("cluster","time","ind","y0")],XL.kprime,by = c("cluster","time","ind"))
      L.jkprime <- matrix(unlist(data.kprime %>% dplyr::select(starts_with("L_"))),ncol = Nperm)
      
      num.diff.mtrx <- t(L.jk) %*% W.list[[cluster.index.k]] %*% matrix(data.k$y0,ncol = 1)%*% 
        matrix(data.kprime$y0,nrow = 1) %*% W.list[[cluster.index.kprime]] %*% L.jkprime
      num.diff.k.kprime <- mean(diag(num.diff.mtrx)/(denom.byperm^2))
      # print(num.diff.k.kprime)
      term2 <- term2 + num.diff.k.kprime
    }
  }
  se <-  sqrt((term1 + 2*term2))
}
