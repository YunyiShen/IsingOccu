
par_name_occu = c("env_spp1","env_spp2","J_intra_spp1","J_intra_spp2","J_inter_spp1","J_inter_spp2","eta12")
par_occu = data.frame(matrix(0,nrow = 1,ncol = 4))
par_occu = par_occu[-1,]
names(par_occu) = c("mean","low","high","point")
par_occu_real = par_occu

# for occu
for(i in (1:length(kk$theta.mcmc))[c(-2,-5)]){
	mean.par = apply(kk$theta.mcmc[[i]],2,mean)
	BI.low.par = apply(kk$theta.mcmc[[i]],2,quantile,probs = .025)
	BI.high.par = apply(kk$theta.mcmc[[i]],2,quantile,probs = .975)
	temp = data.frame(mean.par,BI.low.par,BI.high.par,"Posterior Samples")
	names(temp) = names(par_occu)
	par_occu = rbind(par_occu,temp)
	
	mean.par = (theta[[i]])
	temp = data.frame(mean.par,mean.par,mean.par,"Simulation")
	names(temp) = names(par_occu)
	par_occu_real = rbind(par_occu_real,temp)
}

mean.par = mean(kk$theta.mcmc[[5]][,2])
BI.low.par = quantile(kk$theta.mcmc[[5]][,2],probs = .025)
BI.high.par = quantile(kk$theta.mcmc[[5]][,2],probs = .975)
temp = data.frame(mean.par,BI.low.par,BI.high.par,"Posterior Samples")
names(temp) = names(par_occu)
par_occu = rbind(par_occu,temp)

mean.par = (theta[[5]][2])
temp = data.frame(mean.par,mean.par,mean.par,"Simulation")
names(temp) = names(par_occu)
par_occu_real = rbind(par_occu_real,temp)


par_occu$name = par_name_occu
par_occu_real$name = par_name_occu

par_occu_plot = rbind(par_occu,par_occu_real)



require(ggplot2)
ggplot(par_occu_plot,aes(x=name, y=mean, colour = point)) + 
  geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
  #geom_line() +
  geom_point()+
  theme(axis.text.x = element_text(size = 10, 
                                   color = "black", 
                                   vjust = 1, 
                                   hjust = 1, 
                                   angle = 45))+
  theme(axis.text.y = element_text(size = 10,
                                   color = 'black',
                                   vjust = 0.5,
                                   hjust = 0)) + 
  ylab("value")+
  xlab("parameter")





# for det

mean.par = apply(kk$theta.mcmc[[2]],2,mean)
BI.low.par = apply(kk$theta.mcmc[[2]],2,quantile,probs = .025)
BI.high.par = apply(kk$theta.mcmc[[2]],2,quantile,probs = .975)


par_name = c("det1_spp1","det2_spp1","det3_spp1","det1_spp2","det2_spp2","det3_spp2")
temp1 = data.frame(point = "Posterior Samples"
                   , mean = mean.par
                   , low=BI.low.par
                   , high = BI.high.par
                   , name = par_name )
temp2 = data.frame(point="Simulation"
                   ,mean=theta$beta_det
                   ,low=theta$beta_det
                   ,high=theta$beta_det
                   ,name=par_name)

par_det_plot = rbind(temp1,temp2)


ggplot(par_det_plot,aes(x=name, y=mean, colour = point)) + 
  geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
  #geom_line() +
  geom_point()+
  theme(axis.text.x = element_text(size = 10, 
                                   color = "black", 
                                   vjust = 1, 
                                   hjust = 1, 
                                   angle = 45))+
  theme(axis.text.y = element_text(size = 10,
                                   color = 'black',
                                   vjust = 0.5,
                                   hjust = 0)) + 
  ylab("value")+
  xlab("parameter")



