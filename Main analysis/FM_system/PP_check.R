source("./R/est_logLik.R")
source("./R/misc.R")
source("./R/posterior_prediction.R")
source("./R/Simu_data_Sampling.R")
Rcpp::sourceCpp('./src/IsingCpp_CFTP_sparse.cpp')
require(Matrix)
require(jsonlite)

link = "./data/APIS/"
detmat = list(as.matrix(read.csv(paste0(link,"Fisher_Marten_60dfull_by_islands.csv"),header = F)))

FM_MI_full = read_json("./Main analysis/Results/imperfect_obs/Gibbs_norm/FM/FM_MI_60d_1k.json",T)

FM_MI_full$means$spp_mat <- matrix(FM_MI_full$means$spp_mat,2,2)
FM_MI_full$means$spp_mat_det <- matrix(FM_MI_full$means$spp_mat_det,2,2)

FM_MI_full$linkmap <- lapply(FM_MI_full$linkmap,function(w){Matrix(w,sparse = T)})

FM_MI <- make_list_version_mcmc( FM_MI_full$theta.mcmc,FM_MI_full$means)


set.seed(42)
n_sample <- 2000

samples_ind <- sample(length(FM_MI),n_sample)


FM_MI_sampled <- lapply(samples_ind,function(i,w){w[[i]]},w=FM_MI)

PP_FM_MI <- Posterior_predicting(FM_MI_sampled, FM_MI_full$distM, FM_MI_full$dist_mainland,
                                FM_MI_full$envX, FM_MI_full$linkmap,
                                FM_MI_full$link_mainland,
                                FM_MI_full$interaction.range$intra,FM_MI_full$interaction.range$inter,
                                detX=NULL, nperiod = 17 ,nIter = 100,method = "CFTP")

dets <- abs( detmat[[1]])

PP_FM_MI <- lapply(PP_FM_MI,"*",dets)

mean_det_fisher <- sapply(PP_FM_MI,function(w){
  mean(w[1:155,],na.rm = T)
})  

mean_det_marten <- sapply(PP_FM_MI,function(w){
  mean(w[1:155+155,],na.rm = T)
}) 


########### start checking

# check average in detection

#plot(mean_det_fisher,mean_det_marten)

obs_mean_det_fisher <- mean(detmat[[1]][1:155,],na.rm = T)
obs_mean_det_marten <- mean(detmat[[1]][1:155+155,],na.rm = T)


#points(obs_mean_det_fisher,obs_mean_det_marten
#       ,col = "red")

p_val_grand_fisher <- mean(mean_det_fisher<obs_mean_det_fisher,na.rm = T)
p_val_grand_marten <- mean(mean_det_marten>obs_mean_det_marten,na.rm = T)


  
# naive occupancy 

abs_det_fisher <- sapply(PP_FM_MI,function(w){
  abs_det <- rowSums(w[1:155,]==1,na.rm = T)>0
  mean(abs_det)
})

abs_det_marten <- sapply(PP_FM_MI,function(w){
  abs_det <- rowSums(w[1:155+155,]==1,na.rm = T)>0
  mean(abs_det)
})


obs_abs_det_fisher <- mean(rowSums(detmat[[1]][1:155,]==1,na.rm = T)>0)
obs_abs_det_marten <- mean(rowSums(detmat[[1]][1:155+155,]==1,na.rm = T)>0)


#plot(abs_det_fisher,abs_det_marten)

#points(obs_abs_det_fisher,obs_abs_det_marten
#       ,col = "red")


p_val_abs_fisher <- mean(abs_det_fisher<obs_abs_det_fisher)
p_val_abs_marten <- mean(abs_det_marten<obs_abs_det_marten)



hist(abs_det_fisher,main = "",xlab="Fisher average of site had detection")
abline(v=obs_abs_det_fisher,col = "red")
text(x = .3,y = 300,labels = paste0("p=",signif( p_val_abs_fisher,3)))

hist(abs_det_marten,main = "",xlab="Marten average of site had detection")
abline(v=obs_abs_det_marten,col = "red")
text(x = .5,y = 400,labels = paste0("p=",signif( p_val_abs_marten,3)))




## correlation between absolute

cor_absolute_det <- sapply(PP_FM_MI,function(w){
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

coex <- sapply(PP_FM_MI,function(w){
  sum( rowSums(w[1:155,]==1,na.rm = T)>0 & rowSums(w[1:155+155,]==1,na.rm = T)>0)
  
})


obs_coex <- sapply(detmat,function(w){
  sum( rowSums(w[1:155,]==1,na.rm = T)>0 & rowSums(w[1:155+155,]==1,na.rm = T)>0)
  
})

p_val_coex <- mean(coex>obs_coex)



png("posterior_predictive_check_FM.png",width = 6.25,height=4.8,unit = "in",res=600)
par(mfrow=c(2,3))

hist(mean_det_fisher,main = "",xlab="Fisher grand average of detection")
abline(v=obs_mean_det_fisher,col = "red")
text(x = -0.85,y = 300,labels = paste0("p=",signif(p_val_grand_fisher,3)))

hist(abs_det_fisher,main = "",xlab="Fisher naive occupancy")
abline(v=obs_abs_det_fisher,col = "red")
text(x = .23,y = 300,labels = paste0("p=",signif( p_val_abs_fisher,3)))

hist(cor_absolute_det,main = "",xlab = "Corr of naive occupancy")
abline(v=obs_cor_absolute_det,col = "red")
text(x = .2,y = 300,labels = paste0("p=",signif( p_val_cor,3)))


hist(mean_det_marten,main = "",xlab="Marten grand average of detection")
abline(v=obs_mean_det_marten,col = "red")
text(x = -0.64,y = 700,labels = paste0("p=",p_val_grand_marten))


hist(abs_det_marten,main = "",xlab="Marten naive occupancy")
abline(v=obs_abs_det_marten,col = "red")
text(x = .4,y = 400,labels = paste0("p=",signif( p_val_abs_marten,3)))


hist(coex,main = "",xlab = "Number of confirmed coexistence")
abline(v=obs_coex,col = "red")
text(x = 15,y = 400,labels = paste0("p=",signif( p_val_coex,3)))

dev.off()
  
  


