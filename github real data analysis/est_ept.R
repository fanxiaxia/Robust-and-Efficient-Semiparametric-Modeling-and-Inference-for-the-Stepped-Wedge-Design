est_ept <- function(dataset,delta.ini=delta.0,it.num,ex.type,adj.var.prop,adj.var.out)
{    
  it <- 0
  n.cluster <- length(unique(dataset$cluster))
  cluster.index <- unique(dataset$cluster)
  deltahat <- delta.ini
  while(it < it.num)
  {
    dataset$y_t <- dataset$y - dataset$x*deltahat
    
    # Initialize variables
    fitlmm <- NULL
    model_used <- NULL
    
    ###
    if(is.null(adj.var.out))
    {
      formula_3 <- as.formula(paste("y_t ~ -1 + factor(time) + (1 | cluster)"))
    }
    formula_3 <- as.formula(paste("y_t ~ -1 + factor(time) +", adj.var.out, "+ (1 | cluster)"))
    fitlmm <- try(glmmTMB(formula_3, data = dataset), silent = TRUE)
    #fitlmm <- try(glmmTMB(y_t ~ -1 + factor(time)  + (1 | cluster), data = dataset), silent = TRUE)
    
    {
      # Second model was successful, proceed with variance-covariance calculation
      model_used <- "Model 2: (1 | cluster)"
      
      vc <- VarCorr(fitlmm)$cond
      residual_var <- attr(vc, "sc")^2
      rho <- as.numeric(vc)/(as.numeric(vc) + residual_var)
      
    }
    # Extract gamma
    gamma.df <- data.frame(dataset$cluster,dataset$time, dataset$ind,predict(fitlmm, newdata = dataset))
    colnames(gamma.df) <- c("cluster","time","ind","gamma")
    
    gammafit <- fitlmm
    datamerge <- merge(dataset, gamma.df, by = c("cluster", "time","ind"), all = TRUE)
      
    datacat <- propx.cal(datamerge,ex.type,adj.var.prop)
    deltadenom <- deltanumer <- NULL
    for(i in 1:n.cluster)
    {
      data.i <- datacat %>% filter(cluster == cluster.index[i])
      cluster.i <- as.character(cluster.index[i])
      rho.i <-  rho/(1 + (nrow(data.i)-1)*rho)
      deltadenom <- c(deltadenom,sum(data.i$L * data.i$x) - rho.i * sum(data.i$L) * sum(data.i$x))
      deltanumer <- c(deltanumer,sum(data.i$L * (data.i$y - data.i$gamma)) - rho.i * sum(data.i$L) * sum(data.i$y - data.i$gamma))
    }
      deltahat <- sum(deltanumer)/sum(deltadenom)
      
      it <- it + 1
    }
    
    res <- list(deltahat,gamma.df,rho,gammafit)
    names(res) <- c("est","gamma","rho","gammafit")
    return(res)
}