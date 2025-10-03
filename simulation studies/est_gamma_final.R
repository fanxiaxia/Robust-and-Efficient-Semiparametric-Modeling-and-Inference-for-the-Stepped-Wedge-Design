est_gammaW_crc_adj <- function(deltahat,dataset,gamma.type,adj.var)
{
  n.cluster <- length(unique(dataset$cluster))
  cluster.index <- unique(dataset$cluster)
  if(gamma.type == "cat")
  {
    dataset$y_t <- dataset$y - dataset$x*deltahat
    
    # Initialize variables
    fitlmm <- NULL
    model_used <- NULL
    cluster_vcov_list <- vector("list", n.cluster)  # Preallocate a list of length n.cluster
    W <- NULL
    
    ###
    if(is.null(adj.var))
    {
      formula_3 <- as.formula(paste("y_t ~ -1 + factor(time) + (1 | cluster)"))
    }
    formula_3 <- as.formula(paste("y_t ~ -1 + factor(time) +", adj.var, "+ (1 | cluster)"))
    fitlmm <- try(glmmTMB(formula_3, data = dataset), silent = TRUE)
    #fitlmm <- try(glmmTMB(y_t ~ -1 + factor(time)  + (1 | cluster), data = dataset), silent = TRUE)
    
    {
      # Second model was successful, proceed with variance-covariance calculation
      model_used <- "Model 2: (1 | cluster)"
      
      vc <- VarCorr(fitlmm)$cond
      residual_var <- attr(vc, "sc")^2
      
      # Loop through each cluster (1 to n.cluster) to calculate W
      for (i in 1:n.cluster) {
        data.i <- dataset %>% filter(cluster == cluster.index[i])
        z.i <- as.matrix(rep(1, nrow(data.i)))  # Only intercept
        Vu <- as.matrix(vc$cluster)
        
        # Store the matrix for this specific cluster in the i-th position
        cluster_vcov_list[[i]] <- solve(z.i %*% Vu %*% t(z.i) + diag(residual_var, nrow = nrow(data.i)))
      }
      names(cluster_vcov_list) <- cluster.index
      W <- cluster_vcov_list
    }
    # Extract gamma
    gamma <- data.frame(dataset$cluster,dataset$time, dataset$ind,predict(fitlmm, newdata = dataset))
    colnames(gamma) <- c("cluster","time","ind","gamma")
    
    gamma.fit <- fitlmm
  }
  
  if(gamma.type == "cont")
  {
    dataset$y_t <- dataset$y - dataset$x*deltahat
    
    # Initialize variables
    fitlmm <- NULL
    model_used <- NULL
    cluster_vcov_list <- vector("list", n.cluster)  # Preallocate a list of length n.cluster
    W <- NULL
    
    ###
    if(is.null(adj.var))
    {
      formula_3 <- as.formula(paste("y_t ~ I(time - 1) + (1 | cluster)"))
    }
    formula_3 <- as.formula(paste("y_t ~ I(time - 1) +", adj.var, "+ (1 | cluster)"))
    fitlmm <- try(glmmTMB(formula_3, data = dataset), silent = TRUE)
    #fitlmm <- try(glmmTMB(y_t ~ I(time - 1)  + (1 | cluster), data = dataset), silent = TRUE)
    
    {
      # Second model was successful, proceed with variance-covariance calculation
      model_used <- "Model 2: (1 | cluster)"
      
      vc <- VarCorr(fitlmm)$cond
      residual_var <- attr(vc, "sc")^2
      
      # Loop through each cluster (1 to n.cluster) to calculate W
      for (i in 1:n.cluster) {
        data.i <- dataset %>% filter(cluster == cluster.index[i])
        z.i <- as.matrix(rep(1, nrow(data.i)))  # Only intercept
        Vu <- as.matrix(vc$cluster)
        
        # Store the matrix for this specific cluster in the i-th position
        cluster_vcov_list[[i]] <- solve(z.i %*% Vu %*% t(z.i) + diag(residual_var, nrow = nrow(data.i)))
      }
      names(cluster_vcov_list) <- cluster.index
      W <- cluster_vcov_list
    }
    # Extract gamma
    gamma <- data.frame(dataset$cluster,dataset$time, dataset$ind,predict(fitlmm, newdata = dataset))
    colnames(gamma) <- c("cluster","time","ind","gamma")
    
    gamma.fit <- fitlmm
  }
  
  if(gamma.type == "zero")
  {
    dataset$y_t <- dataset$y - dataset$x*deltahat
    # Initialize variables
    fitlmm <- NULL
    model_used <- NULL
    cluster_vcov_list <- vector("list", n.cluster)  # Preallocate a list of length n.cluster
    W <- NULL
    
    if(is.null(adj.var))
    {
      formula_3 <- as.formula(paste("y_t ~ 1 + (1 | cluster)"))
    }
    formula_3 <- as.formula(paste("y_t ~ 1 +", adj.var, "+ (1 | cluster)"))
    fitlmm_2 <- try(glmmTMB(formula_3, data = dataset), silent = TRUE)
    #fitlmm_2 <- try(glmmTMB(y_t ~ 1 + (1 | cluster), data = dataset), silent = TRUE)
    
    {
      # Second model was successful, proceed with variance-covariance calculation
      fitlmm <- fitlmm_2
      model_used <- "Model 2: (1 + time | cluster)"
      
      vc <- VarCorr(fitlmm)$cond
      residual_var <- attr(vc, "sc")^2
      
      # Loop through each cluster (1 to n.cluster) to calculate W
      for (i in 1:n.cluster) {
        data.i <- dataset %>% filter(cluster == cluster.index[i])
        z.i <- as.matrix(rep(1, nrow(data.i)))  # Only intercept
        Vu <- as.matrix(vc$cluster)
        
        # Store the matrix for this specific cluster in the i-th position
        cluster_vcov_list[[i]] <- solve(z.i %*% Vu %*% t(z.i) + diag(residual_var, nrow = nrow(data.i)))
      }
      names(cluster_vcov_list) <- cluster.index
      W <- cluster_vcov_list
    } 
    # Extract gamma
    gamma <- data.frame(dataset$cluster,dataset$time, dataset$ind,predict(fitlmm, newdata = dataset))
    colnames(gamma) <- c("cluster","time","ind","gamma")
    
    gamma.fit <- fitlmm
  }
  res <- list(gamma,W,model_used,gamma.fit)
  names(res) <- c("gammahat","W","model_used","gammafit")
  return(res)
}

est_gammaW_ctime_adj <- function(deltahat,dataset,gamma.type,adj.var)
{
  n.cluster <- length(unique(dataset$cluster))
  cluster.index <- unique(dataset$cluster)
  if(gamma.type == "cat")
  {
    dataset$y_t <- dataset$y - dataset$x*deltahat
    
    # Initialize variables
    fitlmm <- NULL
    model_used <- NULL
    cluster_vcov_list <- vector("list", n.cluster)  # Preallocate a list of length n.cluster
    W <- NULL
    
    ###
    if(is.null(adj.var))
    {
      formula_2 <- as.formula(paste("y_t ~ -1 + factor(time) + (1 + time | cluster)"))
    }
    formula_2 <- as.formula(paste("y_t ~ -1 + factor(time) +", adj.var, "+ (1 + time | cluster)"))
    fitlmm_2 <- try(glmmTMB(formula_2, data = dataset), silent = TRUE)
    #fitlmm_2 <- try(glmmTMB(y_t ~ -1 + factor(time)  + (1 + time | cluster), data = dataset), silent = TRUE)
    
    {
      # Second model was successful, proceed with variance-covariance calculation
      fitlmm <- fitlmm_2
      model_used <- "Model 2: (1 + time | cluster)"
      
      vc <- VarCorr(fitlmm)$cond
      residual_var <- attr(vc, "sc")^2
      
      # Loop through each cluster (1 to n.cluster) to calculate W
      for (i in 1:n.cluster) {
        data.i <- dataset %>% filter(cluster == cluster.index[i])
        z.i <- as.matrix(cbind(1, data.i$time))  # Only intercept and time
        Vu <- as.matrix(vc$cluster)
        
        # Store the matrix for this specific cluster in the i-th position
        cluster_vcov_list[[i]] <- solve(z.i %*% Vu %*% t(z.i) + diag(residual_var, nrow = nrow(data.i)))
      }
      names(cluster_vcov_list) <- cluster.index
      W <- cluster_vcov_list
    }
    # Extract gamma
    gamma <- data.frame(dataset$cluster,dataset$time, dataset$ind,predict(fitlmm, newdata = dataset))
    colnames(gamma) <- c("cluster","time","ind","gamma")
    
    gamma.fit <- fitlmm
  }
  
  if(gamma.type == "cont")
  {
    dataset$y_t <- dataset$y - dataset$x*deltahat
    
    # Initialize variables
    fitlmm <- NULL
    model_used <- NULL
    cluster_vcov_list <- vector("list", n.cluster)  # Preallocate a list of length n.cluster
    W <- NULL
    
    ###
    if(is.null(adj.var))
    {
      formula_2 <- as.formula(paste("y_t ~ I(time - 1) + (1 + time | cluster)"))
    }
    formula_2 <- as.formula(paste("y_t ~ I(time - 1) +", adj.var, "+ (1 + time | cluster)"))
    fitlmm_2 <- try(glmmTMB(formula_2, data = dataset), silent = TRUE)
    #fitlmm_2 <- try(glmmTMB(y_t ~ I(time - 1)  + (1 + time | cluster), data = dataset), silent = TRUE)
    
    {
      # Second model was successful, proceed with variance-covariance calculation
      fitlmm <- fitlmm_2
      model_used <- "Model 2: (1 + time | cluster)"
      
      vc <- VarCorr(fitlmm)$cond
      residual_var <- attr(vc, "sc")^2
      
      # Loop through each cluster (1 to n.cluster) to calculate W
      for (i in 1:n.cluster) {
        data.i <- dataset %>% filter(cluster == cluster.index[i])
        z.i <- as.matrix(cbind(1, data.i$time))  # Only intercept and time
        Vu <- as.matrix(vc$cluster)
        
        # Store the matrix for this specific cluster in the i-th position
        cluster_vcov_list[[i]] <- solve(z.i %*% Vu %*% t(z.i) + diag(residual_var, nrow = nrow(data.i)))
      }
      names(cluster_vcov_list) <- cluster.index
      W <- cluster_vcov_list
    }
    # Extract gamma
    gamma <- data.frame(dataset$cluster,dataset$time, dataset$ind,predict(fitlmm, newdata = dataset))
    colnames(gamma) <- c("cluster","time","ind","gamma")
    
    gamma.fit <- fitlmm
  }
  
  if(gamma.type == "zero")
  {
    dataset$y_t <- dataset$y - dataset$x*deltahat
    # Initialize variables
    fitlmm <- NULL
    model_used <- NULL
    cluster_vcov_list <- vector("list", n.cluster)  # Preallocate a list of length n.cluster
    W <- NULL
    
    # Fit the first model
    #    fitlmm_1 <- try(glmmTMB(y_t ~ 1 + (1 + factor(time) | cluster), data = dataset), silent = TRUE)
    
    if(is.null(adj.var))
    {
      formula_2 <- as.formula(paste("y_t ~ -1 + (1 + time | cluster)"))
    }
    formula_2 <- as.formula(paste("y_t ~ -1  +", adj.var, "+ (1 + time | cluster)"))
    fitlmm_2 <- try(glmmTMB(formula_2, data = dataset), silent = TRUE)
    #fitlmm_2 <- try(glmmTMB(y_t ~ 1 + (1 + time | cluster), data = dataset), silent = TRUE)
    
    {
      # Second model was successful, proceed with variance-covariance calculation
      fitlmm <- fitlmm_2
      model_used <- "Model 2: (1 + time | cluster)"
      
      vc <- VarCorr(fitlmm)$cond
      residual_var <- attr(vc, "sc")^2
      
      # Loop through each cluster (1 to n.cluster) to calculate W
      for (i in 1:n.cluster) {
        data.i <- dataset %>% filter(cluster == cluster.index[i])
        z.i <- as.matrix(cbind(1, data.i$time))  # Only intercept and time
        Vu <- as.matrix(vc$cluster)
        
        # Store the matrix for this specific cluster in the i-th position
        cluster_vcov_list[[i]] <- solve(z.i %*% Vu %*% t(z.i) + diag(residual_var, nrow = nrow(data.i)))
      }
      names(cluster_vcov_list) <- cluster.index
      W <- cluster_vcov_list
    } 
    # Extract gamma
    gamma <- data.frame(dataset$cluster,dataset$time, dataset$ind,predict(fitlmm, newdata = dataset))
    colnames(gamma) <- c("cluster","time","ind","gamma")
    
    gamma.fit <- fitlmm
  }
  res <- list(gamma,W,model_used,gamma.fit)
  names(res) <- c("gammahat","W","model_used","gammafit")
  return(res)
}