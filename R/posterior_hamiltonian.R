source("misc_island.R")
require(coda)
H = Hamiltonian_posterior(kk$means,envX,distM_full,link_map,distM_mainland,link_mainland =  link_mainland * exp(-distM_mainland),int_range_intra="nn",int_range_inter="nn",Z = Z_sample)
mcmc_iter = nrow(kk$theta.mcmc$beta_occu)
sample_temp = as.list(1:mcmc_iter)

post_para = lapply(sample_temp,function(k,posterior){
  lapply(posterior,function(post,k){post[k,]},k=k)
},posterior = kk$theta.mcmc)

contri = sapply(post_para,Hamiltonian_posterior,
                envX,distM_full,link_map,
                distM_mainland,link_mainland =  link_mainland * exp(-distM_mainland),
                int_range_intra="nn",int_range_inter="nn",Z = Z_sample)

contri = mcmc(t(contri))

mean.par = apply(contri,2,mean)
BI.low.par = apply(contri,2,quantile,probs = .025)
BI.high.par = apply(contri,2,quantile,probs = .975)
Parameter = colnames(contri)

temp1 = data.frame(point = "model fit"
                   , mean = mean.par
                   , low=BI.low.par
                   , high = BI.high.par
                   , name = Parameter )

require(ggplot2)
ggplot(temp1[-c(1,5,9),],aes(x=name, y=mean, colour = point)) + 
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
