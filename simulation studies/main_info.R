source("effectmodels.R")

source("datagen_final.R")
source("clustersize.R")
source("doone_final.R")
source("est_prop_final.R")
source("rdclusters.R")
source("seperm3_final.R")
source("seperm3loo_final.R")
source("sesand2_final.R")
source("ci.R")
source("est_gamma_final.R")
source("libraries.R")

#cluster parameters
set.seed(233)
Ncluster <- 10
Ntime <- 5
#matrix(20, nrow = Ncluster, ncol = Ntime)

#effect parameters
mu_size <- 3
effectsize <- 4
rdcluster.var <- 0.5^2
rderr.var <- 4
rdtimevar <- 0.5^2

#estimation paramters
it.num <- 10
Nperm <- 20
ex.type.set <- "individual emp"
adjusted.var.out <- "basesize"
#adjusted.var.prop <- "factor(time)"
adjusted.var.prop <- "factor(basesize)*factor(time)"

est <- c("lmm_info",
         "prop_info"
) 
loo <- T
do.one(1)

library(parallel)

res1 <- mclapply(1:1000,do.one,mc.cores = detectCores())

save(res1, file = "reslist_final.RData")

