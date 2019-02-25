mean.par = apply(kk$theta.mcmc,2,mean)
BI.low.par = apply(kk$theta.mcmc,2,quantile,probs = .025)
BI.high.par = apply(kk$theta.mcmc,2,quantile,probs = .975)

# par_name = c("env_spp1","env_spp2"
#              ,"det1_spp1","det2_spp1","det3_spp1"
#              ,"det1_spp2","det2_spp2","det3_spp2"
#              ,"J_spp1","useless"
#              ,"J_spp2","useless2"
#              ,"eta_12")

par_name = c("env","J","log(d)")
temp1 = data.frame(point = "model fit"
                   , mean = mean.par
                   , low=BI.low.par
                   , high = BI.high.par
                   , name = par_name )
temp2 = data.frame(point="simulation"
                   ,mean=theta
                   ,low=theta
                   ,high=theta
                   ,name=par_name)

temp = rbind(temp1[-c(10,12),],temp2[-c(10,12),])


require(ggplot2)
ggplot(temp,aes(x=name, y=mean, colour = point)) + 
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
  
