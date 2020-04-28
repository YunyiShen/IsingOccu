source("./R/est_logLik.R")
source("./R/misc.R")
source("./R/posterior_prediction.R")
source("./R/Simu_data_Sampling.R")
Rcpp::sourceCpp('./src/IsingCpp_CFTP_sparse.cpp')
require(Matrix)
require(jsonlite)

link = "./data/APIS/"
detmat = detmat = list(as.matrix(read.csv(paste0(link,"Coyote_Fox_Bobcat_60dfull_by_islands.csv"),header = F)[1:310,]))

CF_MI_full = read_json("./Main analysis/Results/imperfect_obs/Gibbs_norm/CF/CF_MI_60d_0.1.json",T)

CF_MI_full$means$spp_mat <- matrix(CF_MI_full$means$spp_mat,2,2)
CF_MI_full$means$spp_mat_det <- matrix(CF_MI_full$means$spp_mat_det,2,2)

CF_MI_full$linkmap <- lapply(CF_MI_full$linkmap,function(w){Matrix(w,sparse = T)})

CF_MI <- make_list_version_mcmc( CF_MI_full$theta.mcmc,CF_MI_full$means)


set.seed(42)
n_sample <- 2000

samples_ind <- sample(length(CF_MI),n_sample)


CF_MI_sampled <- lapply(samples_ind,function(i,w){w[[i]]},w=CF_MI)

PP_CF_MI <- Posterior_predicting(CF_MI_sampled, CF_MI_full$distM, CF_MI_full$dist_mainland,
                                CF_MI_full$envX, CF_MI_full$linkmap,
                                CF_MI_full$link_mainland,
                                CF_MI_full$interaction.range$intra,CF_MI_full$interaction.range$inter,
                                detX=NULL, nperiod = 17 ,nIter = 50,method = "CFTP")


dets <- abs( detmat[[1]])

PP_CF_MI <- lapply(PP_CF_MI,'*',dets)

mean_det_coyote <- sapply(PP_CF_MI,function(w){
  mean(w[1:155,],na.rm = T)
})  

mean_det_fox <- sapply(PP_CF_MI,function(w){
  mean(w[1:155+155,],na.rm = T)
}) 


########### start checking

# check average in detection

#plot(mean_det_coyote,mean_det_fox)

obs_mean_det_coyote <- mean(detmat[[1]][1:155,],na.rm = T)
obs_mean_det_fox <- mean(detmat[[1]][1:155+155,],na.rm = T)


#points(obs_mean_det_coyote,obs_mean_det_fox
#       ,col = "red")

p_val_grand_coyote <- mean(mean_det_coyote<obs_mean_det_coyote)
p_val_grand_fox <- mean(mean_det_fox>obs_mean_det_fox)


# naive occupancy 

abs_det_coyote <- sapply(PP_CF_MI,function(w){
  abs_det <- rowSums(w[1:155,]==1,na.rm = T)>0
  mean(abs_det)
})

abs_det_fox <- sapply(PP_CF_MI,function(w){
  abs_det <- rowSums(w[1:155+155,]==1,na.rm = T)>0
  mean(abs_det)
})


obs_abs_det_coyote <- mean(rowSums(detmat[[1]][1:155,]==1,na.rm = T)>0)
obs_abs_det_fox <- mean(rowSums(detmat[[1]][1:155+155,]==1,na.rm = T)>0)


#plot(abs_det_coyote,abs_det_fox)

#points(obs_abs_det_coyote,obs_abs_det_fox
#       ,col = "red")


p_val_abs_coyote <- mean(abs_det_coyote<obs_abs_det_coyote)
p_val_abs_fox <- mean(abs_det_fox<obs_abs_det_fox)



hist(abs_det_coyote,main = "",xlab="coyote average of site had detection")
abline(v=obs_abs_det_coyote,col = "red")
text(x = .3,y = 300,labels = paste0("p=",signif( p_val_abs_coyote,3)))

hist(abs_det_fox,main = "",xlab="fox average of site had detection")
abline(v=obs_abs_det_fox,col = "red")
text(x = .5,y = 400,labels = paste0("p=",signif( p_val_abs_fox,3)))




## correlation between absolute

cor_absolute_det <- sapply(PP_CF_MI,function(w){
  cor( rowSums(w[1:155,]==1,na.rm = T)>0, rowSums(w[1:155+155,]==1,na.rm = T)>0)
  
})

obs_cor_absolute_det <- sapply(detmat,function(w){
  cor( rowSums(w[1:155,]==1,na.rm = T)>0, rowSums(w[1:155+155,]==1,na.rm = T)>0)
  
})

p_val_cor <- mean(cor_absolute_det>obs_cor_absolute_det)
hist(cor_absolute_det,main = "",xlab = "Cor of naive occupancy")
abline(v=obs_cor_absolute_det,col = "red")
text(x = .2,y = 300,labels = paste0("p=",signif( p_val_cor,3)))

## number of site coexist

coex <- sapply(PP_CF_MI,function(w){
  sum( rowSums(w[1:155,]==1,na.rm = T)>0 & rowSums(w[1:155+155,]==1,na.rm = T)>0)
  
})


obs_coex <- sapply(detmat,function(w){
  sum( rowSums(w[1:155,]==1,na.rm = T)>0 & rowSums(w[1:155+155,]==1,na.rm = T)>0)
  
})

p_val_coex <- mean(coex>obs_coex)



png("posterior_predictive_check_CF.png",width = 6.25,height=4.8,unit = "in",res=600)
par(mfrow=c(2,3))

hist(mean_det_coyote,main = "",xlab="Coyote grand average of detection")
abline(v=obs_mean_det_coyote,col = "red")
text(x = -0.75,y = 400,
     labels = paste0("p=", 
                      signif(min(p_val_grand_coyote,
                                 1-p_val_grand_coyote),
                             3)
                     ))

hist(abs_det_coyote,main = "",xlab="Coyote naive occupancy")
abline(v=obs_abs_det_coyote,col = "red")
text(x = 0.5,y = 550,
     labels = paste0("p=", 
                     signif(min(p_val_abs_coyote,
                                1-p_val_abs_coyote),
                            3)
     ))

hist(cor_absolute_det,main = "",xlab = "Corr of naive occupancy")
abline(v=obs_cor_absolute_det,col = "red")
text(x = 0.5,y = 310,
     labels = paste0("p=", 
                     signif(min(p_val_cor,
                                1-p_val_cor),
                            3)
     ))


hist(mean_det_fox,main = "",xlab="Fox grand average of detection")
abline(v=obs_mean_det_fox,col = "red")
text(x = -0.9,y = 400,
     labels = paste0("p=", 
                     signif(min(p_val_grand_fox,
                                1-p_val_grand_fox),
                            3)
     ))


hist(abs_det_fox,main = "",xlab="Fox naive occupancy")
abline(v=obs_abs_det_fox,col = "red")
text(x = .22,y = 400,
     labels = paste0("p=", 
                     signif(min(p_val_abs_fox,
                                1-p_val_abs_fox),
                            3)
     ))


hist(coex,main = "",xlab = "Number of confirmed coexistence")
abline(v=obs_coex,col = "red")
text(x = 21,y = 300,
     labels = paste0("p=", 
                     signif(min(p_val_coex,
                                1-p_val_coex),
                            3)
     ))




dev.off()




