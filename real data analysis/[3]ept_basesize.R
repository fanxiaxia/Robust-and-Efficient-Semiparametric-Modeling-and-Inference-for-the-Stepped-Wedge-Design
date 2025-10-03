source("est_prop_final.R")
source("ci.R")
source("libraries.R")
source("est_ept.R")
source("se_perm_ept.R")
source("se_perm_loo_ept.R")

Ncluster <- nrow(diffsize)
Ntime <- ncol(diffsize)

set.seed(2333)
#estimation paramters
it.num <- 10
Nperm <- 50
ex.type.set <- "individual emp"
adjusted.var.out <- "log(basesize)"
adjusted.var.prop <- "log(basesize)*factor(time)"

#adjusted.var.out <- NULL
#adjusted.var.prop <- "factor(time)"

est <- c("lmm_simple",
         "prop_info"
) 
loo <- T


cor_size_x <- cor(df_crosstime$cross_time,df_crosstime$basesize)

# Initialize result placeholders with NA
lmm_bias_s <-  lmmcover_s<-   lmm_bias_info<-   lmmcover_info <-  rep(NA, 4)
prop_bias_simple<- prop_coverp_simple<- prop_covers_simple <- rep(NA, 5)
prop_bias_info<- prop_coverp_info<- prop_covers_info <- rep(NA, 5)


W.list.0 <- West(dataset,rhoest =0)
gamma.ini <- data.frame(dataset$cluster,dataset$time, dataset$ind,0)
colnames(gamma.ini) <- c("cluster","time","ind","gamma")
delta.0 <- est_one_adj(dataset,gamma = gamma.ini, W.list = W.list.0, ex.type = ex.type.set,adj.var.prop = adjusted.var.prop)

#delta.0 <- -0.008978432
if("lmm_simple" %in% est)
{
  
  lmmfit.c <- glmmTMB(y ~  x + factor(time) + log(basesize) + (1|cluster) , data = dataset)
  lmmfit.f <- glmmTMB(y ~  x + factor(time) + log(basesize) +  (1 + factor(time)|cluster), data = dataset)
  # Warning messages:
  #   1: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
  #   Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')
  # 2: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
  #   Model convergence problem; false convergence (8). See vignette('troubleshooting'), help('diagnose')
  
  lmmest_s <- as.numeric(c(summary(lmmfit.c)$coefficients$cond[,"Estimate"]["x"],
                         summary(lmmfit.f)$coefficients$cond[,"Estimate"]["x"]))
  
  # lmm_bias_s <- lmmest_s - truetate
  
  #cover.lm
  lmmci_s <- rbind(
    confint(lmmfit.c)["x",1:2],
    confint(lmmfit.f)["x",1:2])
}



if("prop_info" %in% est)
{
  
  prop.list_info <- est_ept(dataset,delta.ini=delta.0,it.num,ex.type = ex.type.set,adj.var.out = adjusted.var.out, adj.var.prop = adjusted.var.prop)
  
  prop.est_info <- prop.list_info$est
  propgamma_info <- prop.list_info$gamma
  
  prop_seperm_info <- se.perm3_ept(dataset, deltahat = prop.est_info, gammahat =  propgamma_info,
                                   var.cluster = as.numeric(VarCorr(prop.list_info$gammafit)$cond), 
                                   var.error = attr(VarCorr(prop.list_info$gammafit)$cond, "sc")^2, diffsize = diffsize, 
                                   adj.var.prop = adjusted.var.prop,crosstime = df_crosstime$cross_time)
  
  #write.csv(prop.est_info, "est1.csv")
  #write.csv(prop_seperm_info, "seperm1.csv")
  
  # 
  prop_ciperm_info <- calculate_ci(prop.est_info, prop_seperm_info, conf_level = 0.95)
  
  if(loo == T)
  {
    ##leave_one_out
    {
      
      deltahat.list <-  NULL
      gamma.list <- NULL
      # var.cluster.list <- var.error.list <- NULL
      for(l in 1:Ncluster)
      {
        dataset.l <- dataset[which(dataset$cluster != l),]
        datacluster.l <- dataset[which(dataset$cluster == l),]
        
        prop.list_info.l <- est_ept(dataset.l,delta.ini=delta.0,it.num,ex.type = ex.type.set,
                                    adj.var.out = adjusted.var.out, adj.var.prop = adjusted.var.prop)
        
        
        prop.est_info.l <- prop.list_info.l$est
        propgamma_fit.l <- prop.list_info.l$gammafit
        # var.cluster.l <- as.numeric(VarCorr(propgamma_fit.l)$cond)
        # var.error.l <- attr(VarCorr(propgamma_fit.l)$cond, "sc")^2
        
        deltahat.list[[l]] <- prop.est_info.l
        
        gamma.list[[l]] <- data.frame(datacluster.l$cluster,datacluster.l$time, 
                                      datacluster.l$ind,predict(propgamma_fit.l, newdata = datacluster.l))
        
        colnames(gamma.list[[l]]) <- c("cluster","time","ind","gamma")
        
        # var.cluster.list[[l]] <- var.cluster.l
        # var.error.list[[l]] <- var.error.l
        
      }
      #set.seed(122)
      seperm.loo <- se.perm3.loo_ept(dataset, deltahat.list = deltahat.list, gammahat.list = gamma.list,
                                     var.cluster = as.numeric(VarCorr(prop.list_info$gammafit)$cond), 
                                     var.error = attr(VarCorr(prop.list_info$gammafit)$cond, "sc")^2,
                                     diffsize = diffsize, 
                                     adj.var.prop = adjusted.var.prop,crosstime = df_crosstime$cross_time)
      
      
      ciperm.loo <- calculate_ci(prop.est_info, seperm.loo, conf_level = 0.95)
      
    }
  }
}


res <- rbind(c(prop.est_info, prop_seperm_info, prop_ciperm_info),
             c(NA, seperm.loo, ciperm.loo))

colnames(res) <- c("est", "SE", "lower", "upper")
rownames(res) <- c("permuation-based","Leave-one-out")

res.lmm <- data.frame(lmmest_s, c(as.numeric(summary(lmmfit.c)$coefficients$cond[,"Std. Error"]["x"]),
                                  as.numeric(summary(lmmfit.f)$coefficients$cond[,"Std. Error"]["x"])),lmmci_s[,1],lmmci_s[,2])

colnames(res.lmm) <- c("est", "SE", "lower", "upper")

write.csv(res,"res_basesize.csv")
write.csv(res.lmm,"res_lmm_basesize.csv")



