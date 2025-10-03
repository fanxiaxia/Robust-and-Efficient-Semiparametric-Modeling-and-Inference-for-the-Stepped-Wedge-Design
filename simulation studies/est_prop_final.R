propx.cal <- function(dataset,ex.type,adj.var.prop)
{
  if(ex.type == "cluster emp")
  {
    coldata <- dataset %>% group_by(cluster,time) %>% summarise(nij = n(),xij = mean(x))
    formula <- as.formula(paste("xij ~", adj.var.prop))
    fit.cluster <- glm(formula, data = coldata, family = binomial(link = "logit"))
    dataset$meanx <- predict(fit.cluster, type = "response")
  }
  if(ex.type == "individual emp")
  {
    formula <- as.formula(paste("x ~", adj.var.prop))
    fit.ind <- glm(formula, data = dataset, family = binomial(link = "logit"))
    dataset$meanx <- predict(fit.ind, type = "response")
  }
  
  dataset$L <- dataset$x - dataset$meanx
  
  return(dataset)
}

West <- function(dataset,rhoest)
{
  n.cluster <- length(unique(dataset$cluster))
  cluster.index <- unique(dataset$cluster)
  
  W.list <- NULL
  for(i in 1:n.cluster)
  {
    dataset.i <- dataset %>% filter(cluster == cluster.index[i])
    ni <- nrow(dataset.i)
    
    a <- (ni*(1+(ni-1)*rhoest))^(-1) -(ni*(1-rhoest))^{-1}
    b <- (1-rhoest)^(-1)
    
    Wi <- diag(rep(b),ni)+matrix(a,nrow =ni, ncol = ni)
    
    W.list[[length(W.list)+1]] <- Wi
    
  }
  
  names(W.list) <- cluster.index
  W.list
}

est_one_adj <- function(dataset, gamma.df,W.list, ex.type = "individual emp",adj.var.prop)
{
  n.cluster <- length(unique(dataset$cluster))
  cluster.index <- unique(dataset$cluster)
  
  datamerge <- merge(dataset, gamma.df, by = c("cluster", "time","ind"), all = TRUE)
  
  datacat <- propx.cal(datamerge,ex.type,adj.var.prop)
  deltadenom <- deltanumer <- NULL
  for(i in 1:n.cluster)
  {
    data.i <- datacat %>% filter(cluster == cluster.index[i])
    cluster.i <- as.character(cluster.index[i])
    deltadenom <- c(deltadenom, matrix(data.i$L,nrow = 1) %*% W.list[[cluster.i]] %*% matrix(data.i$x, ncol = 1))
    deltanumer <- c(deltanumer, matrix(data.i$L,nrow = 1)  %*% W.list[[cluster.i]] %*%  matrix(data.i$y - data.i$gamma, ncol = 1))
  }
  est <- sum(deltanumer)/sum(deltadenom)
  
  return(est)
}


est_crc_newW <- function(dataset,delta.ini,it.num,ex.type = "individual emp",gamma.type = "zero",adj.var.out,adj.var.prop)
{
  it <- 0
  deltahat <- delta.ini
  while(it < it.num)
  {
    res.it <- est_gammaW_crc_adj(deltahat,dataset,gamma.type = gamma.type,adj.var = adj.var.out)
    W <- res.it$W
    gammahat <- res.it$gammahat
    deltahat <- est_one_adj(dataset,gamma.df=gammahat, W.list= W, ex.type = ex.type,adj.var.prop)
    
    gammafit <- res.it$gammafit
    it <- it + 1
  }
  
  res <- list(deltahat,gammahat,W,gammafit)
  names(res) <- c("est","gamma","W","gammafit")
  return(res)
}


est_ctime_newW <- function(dataset,delta.ini,it.num,ex.type = "individual emp",gamma.type = "zero",adj.var.out,adj.var.prop)
{
  it <- 0
  deltahat <- delta.ini
  
  while(it < it.num)
  {
    res.it <- est_gammaW_ctime_adj(deltahat,dataset,gamma.type = gamma.type,adj.var.out)
    W <- res.it$W
    gammahat <- res.it$gammahat
    deltahat <- est_one_adj(dataset,gamma.df=gammahat, W.list= W, ex.type = ex.type,adj.var.prop)
    
    gammafit <- res.it$gammafit
    it <- it + 1
  }
  
  res <- list(deltahat,gammahat,W,gammafit)
  names(res) <- c("est","gamma","W","gammafit")
  return(res)
  }
  


