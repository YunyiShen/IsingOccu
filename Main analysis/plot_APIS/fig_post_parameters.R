source("./R/misc.R")
require(coda)

Fisher = sapply(kk$theta.mcmc,function(kk){kk[,1]})
Fisher = Fisher[,-c(1,2,5,6)]
colnames(Fisher) = paste0("Fisher_",colnames(Fisher))

Marten = sapply(kk$theta.mcmc,function(kk){kk[,2]})

Inter = Marten[,c(5,6)]
Marten = Marten[,-c(1,2,5,6)]
colnames(Marten) = paste0("Marten_",colnames(Marten))
colnames(Inter) = c("Association_occu","Association_det")

parameters = data.frame(Fisher,Marten,Inter)

require(reshape2)
require(ggplot2)
temp = melt(parameters,value.name = "value")
ggplot(data = temp,aes(x=variable,y=value))+geom_boxplot()+
  theme(axis.text.x = element_text(size = 10, 
                                   color = "black", 
                                   vjust = 1, 
                                   hjust = 1, 
                                   angle = 45))+
  theme(axis.text.y = element_text(size = 10,
                                   color = 'black',
                                   vjust = 0.5,
                                   hjust = 0)) + 
  ylab("model fit")+
  xlab("parameter")
