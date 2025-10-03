generate_dataset <- function(mu, tau, theta, n_clusters, n_time_points,
                             n_ind_per_cluster, sigma,
                             delay_model, 
                             time_trend="none",rdtmsd,
                             same.trend = TRUE,
                             time_size_int = FALSE) 
{
  
  # Generate data frame
  data <- data.frame(
    "i" = integer(), # cluster
    "j" = integer(), # time; 1=baseline, J=endline
    "k" = integer(), # individual
    "l" = integer(), # time since intervention
    "x_ij" = integer(), # treatment state indicator
    "c_i" = integer(), # the start time of the treatment
    "y" = integer() # binary outcome
  )
  
  # Generate crossover times (assumes a "as balanced as possible and complete" design)
  crossover_times <-  #sample(rep(2:Ntime, length.out=Ncluster))
    randomize_clusters(Ncluster, 2:Ntime)$cross
  
  #time trend 1
  beta_js1 <- sapply(1:n_time_points, function(j) {
    if (time_trend == "incr") (-j + 1)
    else if (time_trend == "none") 0
    else stop("time_trend trend misspecified")
  })
  
  #time trend 2
  beta_js2 <- sapply(1:n_time_points, function(j) {
    if (time_trend == "incr") 0
    else if (time_trend == "none") 0
    else stop("time_trend trend misspecified")
  })
  
  # Create theta_ls (intervention effects) based on continuous fn "delay_model"
  theta_ls <- theta * effect_curve(
    x = 1:(n_time_points-1),
    type = delay_model$type,
    params = delay_model$params
  )
  
  # Loop through clusters, time, and individuals
  for (i in 1:n_clusters) {
    #time trend latent variable
    if(same.trend)
    {
      Ui <- 1
    }
    if(!same.trend)
    {
      Ui <- ifelse(i <= 4,1,0)
    }
    
    Ki <- n_ind_per_cluster[i,]
    #   Vi <- diffsize2[i,]
    rdcleff_i <- rnorm(1, mean = 0, sd = tau)
    rdtmeff.ij <- rnorm(1, mean=0, sd=rdtmsd)
    
    for (j in 1:(n_time_points)) {
      
      Kij <- Ki[j]
      
      x_ij <- ifelse(j<crossover_times[i], 0, 1)
      l <- ifelse(j<crossover_times[i], 0, (j-crossover_times[i])+1)
      theta_l <- ifelse(l>0, theta_ls[l], 0)
      
      beta_js <- Ui*beta_js1[j] + (1-Ui)*beta_js2[j]
      
      mu_ij <- mu +  4*(beta_js)^2 + (theta_l)*x_ij +  
      ifelse(time_size_int, ifelse(Ki[1]<=15, -4*(beta_js)^2,0),0) +
        rdcleff_i + rdtmeff.ij*j
      
      
      y <- mu_ij + rnorm(Kij, mean=0, sd=sigma) 
      
      
      data <- rbind(data, data.frame(cbind(
        i=rep(i, Kij), j=rep(j, Kij), k=c(1:Kij), l=rep(l, Kij), # alpha_i=rep(alpha_i,k)
        x_ij=rep(x_ij, Kij), y=y,size = rep(Kij, Kij), basesize = rep(Ki[1],Kij) # mu_ij=rep(mu_ij,k)
      )))
      
    }
  }
  
  return (list(
    "params" = list(
      n_clusters = n_clusters,
      n_time_points = n_time_points,
      crossover_times = crossover_times
    ),
    "beta_js" = cbind(beta_js1,beta_js2),
    "theta_ls" = theta_ls,
    "data" = data
  ))
  
}

