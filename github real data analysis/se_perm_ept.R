se.perm3_ept <- function(dataset,deltahat,gammahat,var.cluster,var.error,diffsize,adj.var.prop,crosstime)
{
  dataset <- merge(dataset, gammahat, by = c("cluster","time","ind"), all = TRUE)
  dataset <- dataset %>% mutate(y0 = y - x*deltahat - gamma)
  dataset <- dataset %>% arrange(cluster,time,ind)
  
  rho <- var.cluster/(var.cluster + var.error)
  constant <- 1/((var.cluster + var.error)*(1-rho))
  
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
    
    rho.i <- rho/(1 + (nrow(data.i)-1)*rho)
    # W.temp <- constant*(diag(1,nrow = nrow(data.i)) - rho.i*matrix(1, nrow = nrow(data.i), ncol = nrow(data.i)))
    # denom.i <- diag(t(matrix(unlist(L.i),ncol = Nperm)) %*% W.temp %*% matrix(unlist(Xu.i),ncol = Nperm))
    
    L_mat <- matrix(unlist(L.i), ncol = Nperm)   # dimensions: n x Nperm
    X_mat <- matrix(unlist(Xu.i), ncol = Nperm)    # dimensions: n x Nperm
    
    denom.i <- constant * (colSums(L_mat * X_mat) - rho.i * colSums(L_mat) * colSums(X_mat))
    
    denom.vec[i,] <- denom.i
    
    rm(data.i, XL.i, L_mat, X_mat)  # remove large objects
    gc()  # request garbage collection
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
    
     rho.i <- rho/(1 + (nrow(data.i)-1)*rho)
    # W.temp <- constant*(diag(1,nrow = nrow(data.i)) - rho.i*matrix(1, nrow = nrow(data.i), ncol = nrow(data.i)))
    # num.i <- diag(t(L.i) %*% W.temp %*% matrix(data.i$y0,ncol = 1) %*% 
    #                 matrix(data.i$y0,nrow = 1) %*% W.temp %*% L.i)
    
    # Compute the "adjusted" y0 vector without forming an n x n matrix:
    v <- data.i$y0 - rho.i * sum(data.i$y0)
    v <- constant * v  # Multiply by the constant factor
    
    # Multiply by L.i (assuming L.i is a matrix with appropriate dimensions)
    w <- t(L.i) %*% v
    
    # The diagonal of w %*% t(w) is simply the elementwise square of w.
    num.i <- (w)^2  # This gives a vector (if w is a column vector)
    term1 <- term1 + mean(num.i/(denom.byperm^2))
    
    rm(data.i, XL.i, L.i)  # remove large objects
    gc()  # request garbage collection
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
      
       rho.k <- rho/(1 + (nrow(data.k)-1)*rho)
      # W.temp.k <- constant*(diag(1,nrow = nrow(data.k)) - rho.k*matrix(1, nrow = nrow(data.k), ncol = nrow(data.k)))
       rho.kprime <- rho/(1 + (nrow(data.kprime)-1)*rho)
      # W.temp.kprime <- constant*(diag(1,nrow = nrow(data.kprime)) - 
      #   rho.kprime*matrix(1, nrow = nrow(data.kprime), ncol = nrow(data.kprime)))
      # 
      # num.diff.mtrx <- t(L.jk) %*% W.temp.k %*% matrix(data.k$y0,ncol = 1)%*% 
      #   matrix(data.kprime$y0,nrow = 1) %*% W.temp.kprime %*% L.jkprime
      # num.diff.k.kprime <- mean(diag(num.diff.mtrx)/(denom.byperm^2))
      
      # Compute the “weighted” y0 vectors:
      v_k <- constant * ( data.k$y0 - rho.k * sum(data.k$y0) )
      v_kprime <- constant * ( data.kprime$y0 - rho.kprime * sum(data.kprime$y0) )
      
      # Now form the outer product (v_k %*% t(v_kprime)) – this is equivalent to applying the weights to the outer product of y0’s:
      outer_v <- v_k %*% t(v_kprime)
      
      # Then, form the product with L matrices:
      num.diff.mtrx <- t(L.jk) %*% outer_v %*% L.jkprime
      
      # Finally, compute your desired statistic:
      num.diff.k.kprime <- mean(diag(num.diff.mtrx)/(denom.byperm^2))
      # print(num.diff.k.kprime)
      term2 <- term2 + num.diff.k.kprime
      
      rm(data.k, XL.k, L.jk,data.kprime, XL.kprime, L.jkprime)  # remove large objects
      gc()  # request garbage collection
    }
  }
  se <-  sqrt((term1 + 2*term2))
}
