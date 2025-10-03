do.one <- function(trials)
{
  set.seed(trials+2333)
  diffsize <- generate_increasing_N_size(Ncluster = Ncluster, Ntime = Ntime)
  print(diffsize)
  data <- generate_dataset(
    mu = mu_size,
    n_clusters = Ncluster,
    n_time_points = Ntime,
    n_ind_per_cluster = diffsize,
    theta = effectsize,
    tau = sqrt(rdcluster.var), 
    sigma = sqrt(rderr.var),
    rdtmsd = sqrt(rdtimevar),
    # delay_model = list(type="spline", params=list(knots=c(0,1),slopes=1)),
    # delay_model = list(type="exp", params=list(d=1.5)),
    delay_model = instant.model,
    # n_extra_time_points = 2,
    # rte = list(type="height", rho=-0.2, nu=0.4)
    # rte = list(type="height+shape", rho1=-0.1, rho2=0.6, nu=0.4)
    time_trend = "incr",
    same.trend = T,
    time_size_int = T) 
  
  truetate <- mean(data$theta_ls)
  
  sim <- data$data
  sim$cluster <- sim$i
  sim$time <- sim$j
  sim$x <- sim$x_ij
  sim$y <- sim$y
  sim$ind <- sim$k
  dataset <- sim[,c("cluster","time","x","y","l","k","ind","size","basesize")]
  
  par.set <- data$params
  cor_size_x <- cor(par.set$crossover_times,diffsize[,1])
  sum.denom <- NA #par.set$sum.denom
  #cor_size_x <- cor(sim$xdb,sim$x)
  #cor_size_x
  #cor(par.set$crossover_times,diffsize[,1])
  
  # Initialize result placeholders with NA
  lmm_bias_s <-  lmmcover_s<-   lmm_bias_info<-   lmmcover_info <-  rep(NA, 4)
  prop_bias_simple<- prop_coverp_simple<- prop_covers_simple <- rep(NA, 5)
  prop_bias_info<- prop_coverp_info<- prop_covers_info <- rep(NA, 5)
  
  
  W.list.0 <- West(dataset,rhoest =0)
  gamma.ini <- data.frame(dataset$cluster,dataset$time, dataset$ind,0)
  colnames(gamma.ini) <- c("cluster","time","ind","gamma")
  delta.0 <- est_one_adj(dataset,gamma = gamma.ini, W.list = W.list.0, ex.type = ex.type.set,adj.var.prop = adjusted.var.prop)
  
  if("lmm_simple" %in% est)
  {
    #LMM.it without random intercpet (lm)
    lmmfit.a <- glmmTMB(y ~  x_ij + I(j-1) + (1 |i) , data = data$data)
    lmmfit.b <- glmmTMB(y ~  x_ij + I(j-1) + (1 + j |i) , data = data$data)
    lmmfit.c <- glmmTMB(y ~  x_ij + factor(j) + (1 |i) , data = data$data)
    lmmfit.d <- glmmTMB(y ~  x_ij + factor(j) + (1 + j |i) , data = data$data)
    
    lmmest_s <- c(as.numeric(summary(lmmfit.a)$coefficients$cond[,"Estimate"]["x_ij"]),
                  as.numeric(summary(lmmfit.b)$coefficients$cond[,"Estimate"]["x_ij"]),
                  as.numeric(summary(lmmfit.c)$coefficients$cond[,"Estimate"]["x_ij"]),
                  as.numeric(summary(lmmfit.d)$coefficients$cond[,"Estimate"]["x_ij"]))
    
    lmm_bias_s <- lmmest_s - truetate
    
    #cover.lm
    lmmci_s <- rbind(confint(lmmfit.a)["x_ij",1:2],
                     confint(lmmfit.b)["x_ij",1:2],
                     confint(lmmfit.c)["x_ij",1:2],
                     confint(lmmfit.d)["x_ij",1:2]
    )
    lmmcover_s <- apply(lmmci_s, 1, function(x) truetate >= x[1] && truetate <= x[2])
  }
  
  if("lmm_info" %in% est)
  {
    #LMM.it without random intercpet (lm)
    lmmfit.e0 <- glmmTMB(y ~  x_ij + j + basesize +  (1 |i), data = data$data)
    lmmfit.e <- glmmTMB(y ~  x_ij + j+ basesize + (1 + j|i), data = data$data)
    lmmfit.f0 <- glmmTMB(y ~  x_ij + factor(j) + basesize + (1 |i), data = data$data)
    lmmfit.f <- glmmTMB(y ~  x_ij + factor(j) + basesize + (1 + j|i), data = data$data)
    
    lmmest_info <- c(as.numeric(summary(lmmfit.e0)$coefficients$cond[,"Estimate"]["x_ij"]),
                     as.numeric(summary(lmmfit.e)$coefficients$cond[,"Estimate"]["x_ij"]),
                     as.numeric(summary(lmmfit.f0)$coefficients$cond[,"Estimate"]["x_ij"]),
                     as.numeric(summary(lmmfit.f)$coefficients$cond[,"Estimate"]["x_ij"]))
    
    lmm_bias_info <- lmmest_info - truetate
    
    #cover.lm
    lmmci_info <- rbind(confint(lmmfit.e0)["x_ij",1:2],
                        confint(lmmfit.e)["x_ij",1:2],
                        confint(lmmfit.f0)["x_ij",1:2],
                        confint(lmmfit.f)["x_ij",1:2]
    )
    lmmcover_info <- apply(lmmci_info, 1, function(x) truetate >= x[1] && truetate <= x[2])
  }
  
  if("prop_info" %in% est)
  {
    prop.list_info <- list(
      est_crc_newW(dataset,delta.ini=delta.0,it.num,ex.type = ex.type.set,gamma.type = "cont",adj.var.out = adjusted.var.out, adj.var.prop = adjusted.var.prop),
      est_ctime_newW(dataset,delta.ini=delta.0,it.num,ex.type = ex.type.set,gamma.type = "cont",adj.var.out = adjusted.var.out, adj.var.prop = adjusted.var.prop),
      est_crc_newW(dataset,delta.ini=delta.0,it.num,ex.type = ex.type.set,gamma.type = "cat",adj.var.out = adjusted.var.out, adj.var.prop = adjusted.var.prop),
      est_ctime_newW(dataset,delta.ini=delta.0,it.num,ex.type = ex.type.set,gamma.type = "cat",adj.var.out = adjusted.var.out, adj.var.prop = adjusted.var.prop)
    )
    
    prop.est_info <- unlist(lapply(prop.list_info, function(model) model$est))
    propgamma_info <- lapply(prop.list_info, function(model) model$gamma)
    propW_info <- lapply(prop.list_info, function(model) model$W)
    
    
    prop_bias_info <- prop.est_info - truetate
    
    
    # ##coverage_ctime
    prop_seperm_info <- c(
      se.perm3(dataset, deltahat =  prop.est_info[1], gammahat = propgamma_info[[1]],
               W.list = propW_info[[1]], diffsize = diffsize, adj.var.prop = adjusted.var.prop, crosstime = data$params$crossover_times),
      se.perm3(dataset, deltahat = prop.est_info[2], gammahat =  propgamma_info[[2]],
               W.list = propW_info[[2]], diffsize = diffsize, adj.var.prop = adjusted.var.prop,crosstime = data$params$crossover_times),
      se.perm3(dataset, deltahat = prop.est_info[3], gammahat =  propgamma_info[[3]],
               W.list = propW_info[[3]], diffsize = diffsize, adj.var.prop = adjusted.var.prop,crosstime = data$params$crossover_times),
      se.perm3(dataset, deltahat = prop.est_info[4], gammahat =  propgamma_info[[4]],
               W.list = propW_info[[4]], diffsize = diffsize, adj.var.prop = adjusted.var.prop,crosstime = data$params$crossover_times)
    )
    # # 
    # 
    # 
     prop_ciperm_info <- t(apply(rbind(prop.est_info, prop_seperm_info),2, function(x){calculate_ci(x[1], x[2], conf_level = 0.95)}))
     prop_coverp_info <- apply(prop_ciperm_info, 1, function(x) truetate >= x[1] && truetate <= x[2])      
    # 
    prop_sesand_info <-   c(
      se.sand2(dataset, deltahat =  prop.est_info[1], gammahat = propgamma_info[[1]],
               W.list = propW_info[[1]], diffsize = diffsize, adj.var.prop = adjusted.var.prop),
      se.sand2(dataset, deltahat = prop.est_info[2], gammahat = propgamma_info[[2]],
               W.list = propW_info[[2]], diffsize = diffsize, adj.var.prop = adjusted.var.prop),
      se.sand2(dataset, deltahat = prop.est_info[3], gammahat = propgamma_info[[3]],
               W.list = propW_info[[3]], diffsize = diffsize, adj.var.prop = adjusted.var.prop),
      se.sand2(dataset, deltahat = prop.est_info[4], gammahat = propgamma_info[[4]],
               W.list = propW_info[[4]], diffsize = diffsize, adj.var.prop = adjusted.var.prop)
    )
    # 
    # 
     prop_cisand_info <- t(apply(rbind(prop.est_info, prop_sesand_info),2, function(x){calculate_ci(x[1], x[2], conf_level = 0.95)}))
     prop_covers_info <- apply(prop_cisand_info, 1, function(x) truetate >= x[1] && truetate <= x[2])  
     
     if(loo == T)
     {
       ##leave_one_out
       {
         
         deltahat.list1 <- deltahat.list2 <- deltahat.list3 <- deltahat.list4 <- NULL
         gamma.list1 <- gamma.list2 <- gamma.list3 <- gamma.list4 <- NULL
         for(l in 1:Ncluster)
         {
           dataset.l <- dataset[which(dataset$cluster != l),]
           datacluster.l <- dataset[which(dataset$cluster == l),]
           
           prop.list_info.l <- list(
             est_crc_newW(dataset.l,delta.ini=delta.0,it.num,ex.type = ex.type.set,gamma.type = "cont",adj.var.out = adjusted.var.out, adj.var.prop = adjusted.var.prop),
             est_ctime_newW(dataset.l,delta.ini=delta.0,it.num,ex.type = ex.type.set,gamma.type = "cont",adj.var.out = adjusted.var.out, adj.var.prop = adjusted.var.prop),
             est_crc_newW(dataset.l,delta.ini=delta.0,it.num,ex.type = ex.type.set,gamma.type = "cat",adj.var.out = adjusted.var.out, adj.var.prop = adjusted.var.prop),
             est_ctime_newW(dataset.l,delta.ini=delta.0,it.num,ex.type = ex.type.set,gamma.type = "cat",adj.var.out = adjusted.var.out, adj.var.prop = adjusted.var.prop)
           )
           
           prop.est_info.l <- unlist(lapply(prop.list_info.l, function(model) model$est))
           propgamma_fit.l <- lapply(prop.list_info.l, function(model) model$gammafit)
           #propW_info.l <- lapply(prop.list_info, function(model) model$W)
           
           deltahat.list1[[l]] <- prop.est_info.l[1]
           deltahat.list2[[l]] <- prop.est_info.l[2]
           deltahat.list3[[l]] <- prop.est_info.l[3]
           deltahat.list4[[l]] <- prop.est_info.l[4]
           
           gamma.list1[[l]] <- data.frame(datacluster.l$cluster,datacluster.l$time, datacluster.l$ind,predict(propgamma_fit.l[[1]], newdata = datacluster.l))
           gamma.list2[[l]] <- data.frame(datacluster.l$cluster,datacluster.l$time, datacluster.l$ind,predict(propgamma_fit.l[[2]], newdata = datacluster.l))
           gamma.list3[[l]] <- data.frame(datacluster.l$cluster,datacluster.l$time, datacluster.l$ind,predict(propgamma_fit.l[[3]], newdata = datacluster.l))
           gamma.list4[[l]] <- data.frame(datacluster.l$cluster,datacluster.l$time, datacluster.l$ind,predict(propgamma_fit.l[[4]], newdata = datacluster.l))
           
           colnames(gamma.list1[[l]]) <- c("cluster","time","ind","gamma")
           colnames(gamma.list2[[l]]) <- c("cluster","time","ind","gamma")
           colnames(gamma.list3[[l]]) <- c("cluster","time","ind","gamma")
           colnames(gamma.list4[[l]]) <- c("cluster","time","ind","gamma")
           
         }
         
         seperm.loo <- c(se.perm3.loo(dataset, deltahat.list = deltahat.list1, gammahat.list = gamma.list1,
                                         W.list = propW_info[[1]], diffsize = diffsize, adj.var.prop = adjusted.var.prop,crosstime = data$params$crossover_times),
                            se.perm3.loo(dataset, deltahat.list = deltahat.list2, gammahat.list = gamma.list2,
                                         W.list = propW_info[[2]], diffsize = diffsize, adj.var.prop = adjusted.var.prop,crosstime = data$params$crossover_times),
                            se.perm3.loo(dataset, deltahat.list = deltahat.list3, gammahat.list = gamma.list3,
                                         W.list = propW_info[[3]], diffsize = diffsize, adj.var.prop = adjusted.var.prop,crosstime = data$params$crossover_times),
                            se.perm3.loo(dataset, deltahat.list = deltahat.list4, gammahat.list = gamma.list4,
                                         W.list = propW_info[[4]], diffsize = diffsize, adj.var.prop = adjusted.var.prop,crosstime = data$params$crossover_times)
         )
         
         ciperm.loo <- t(apply(rbind(prop.est_info, seperm.loo),2, function(x){calculate_ci(x[1], x[2], conf_level = 0.95)}))
         coverp.loo <- apply(ciperm.loo, 1, function(x) truetate >= x[1] && truetate <= x[2])

       }
     }
  }
  
  cor.list <- c(cor_size_x,sum.denom)
  
  res.new.list <- list(cor.list,lmm_bias_s, lmmcover_s, lmm_bias_info, lmmcover_info,
                       prop_bias_simple, #prop_coverp_simple, prop_covers_simple,
                       prop_bias_info, prop_coverp_info, prop_covers_info,coverp.loo
  )

  return(res.new.list)
}