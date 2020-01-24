Z1_map
ggsave("flip_Z1.jpg",width = 18,height = 12,dpi=500,units = "cm")
Z2_map
ggsave("flip_Z2.jpg",width = 18,height = 12,dpi=500,units = "cm")

plotdata = data.frame(matrix(NA,12,5))

colnames(plotdata) = c("parameters","Value","lw","hi","point")
plotdata$parameters = rep(c("mainland_1","mainland_2","intra_1","intra_2","gamma_oc","gamma_det"),2)

plotdata$point = c(rep("Estimation(90%CI)",6),rep("Ground_truth",6))

plotdata$Value[1:2+6] = theta$eta_inter
plotdata$Value[3:4+6] = theta$eta_intra

plotdata$Value[5+6] = theta$spp_mat[1,2]
plotdata$Value[6+6] = theta$spp_mat_det[1,2]

plotdata$lw[1:6+6] <- plotdata$hi[1:6+6] <- plotdata$Value[1:6+6]

plotdata$Value[1:2] = apply(kk$theta.mcmc$eta_inter,2,quantile,0.5)
plotdata$lw[1:2] = apply(kk$theta.mcmc$eta_inter,2,quantile,0.05)
plotdata$hi[1:2] = apply(kk$theta.mcmc$eta_inter,2,quantile,0.95)

plotdata$Value[3:4] = apply(kk$theta.mcmc$eta_intra,2,quantile,0.5)
plotdata$lw[3:4] = apply(kk$theta.mcmc$eta_intra,2,quantile,0.05)
plotdata$hi[3:4] = apply(kk$theta.mcmc$eta_intra,2,quantile,0.95)


plotdata$Value[5] = quantile(kk$theta.mcmc$spp_mat[,2],0.5)
plotdata$lw[5] = quantile(kk$theta.mcmc$spp_mat[,2],0.05)
plotdata$hi[5] = quantile(kk$theta.mcmc$spp_mat[,2],0.95)

plotdata$Value[6] = quantile(kk$theta.mcmc$spp_mat_det[,2],0.5)
plotdata$lw[6] = quantile(kk$theta.mcmc$spp_mat_det[,2],0.05)
plotdata$hi[6] = quantile(kk$theta.mcmc$spp_mat_det[,2],0.95)



ggplot(data = plotdata,mapping = aes(x=parameters,y=Value,color=point))+
  geom_point(size=2.5)+
  geom_errorbar(aes(ymin=lw, ymax=hi) , width=.1) +
  theme(text = element_text(size=14), axis.text.x = element_text(angle=45, hjust=1,size = 12))

ggsave("competition_para.jpg",width = 18,height = 12,dpi = 500,units = "cm")
