vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)


## not used
myboxplot <- function(x, data, col = NULL, xlab) {
    boxplot(x, data, axes = FALSE, col = col,xlab = NULL)
    #axis(1, at = c(1,2), labels =FALSE)
    text(1.5, y=par()$usr[3]-0.08*(par()$usr[4]-par()$usr[3]),
         srt=60, xpd=T, adj=1, labels = xlab)
}

test_data <- data.frame(value=rnorm(100),
                        id = sample(c("T","F"),100,replace = TRUE),
                        group = sample(c("T","F"),100,replace = TRUE))

par(mfrow = c(2,2))

layout(t(1:5))
par(oma=c(1, 1, 1, 1), mar=c(5,2,1,1), cex=1)

myboxplot(value~id,test_data,xlab = "foo")
axis(2, las=1)
myboxplot(value~id,test_data,xlab = "foo2")
myboxplot(value~id,test_data,xlab = "foo3")
myboxplot(value~id,test_data,xlab = "foo4")
myboxplot(value~id,test_data,xlab = "foo5")

# use this
library(ggplot2)
true_setting_noAPIS <- read.csv("./Main analysis/Tests/100_datasets/true_setting_noAPIS.csv",row.names = 1)
true_setting_APIS <- read.csv("./Main analysis/Tests/100_datasets/true_setting_APIS.csv",row.names = 1)

big_simulation_root <- "./Main analysis/Results/Big_simulation"
grid_sys <- c(10,15,20,25,"APIS")
setting <- c("C","N","S")
parameters <- c("beta_1","beta_2","eta_intra_1","eta_intra_2","gamma_oc","gamma_de")

all_settings <- expand.grid(big_simulation_root,grid_sys,setting,parameters,stringsAsFactors = F)
all_settings$para <- sub("_2","" ,sub("_1","" ,all_settings$Var4))
all_settings$spp <- sub("beta_","" ,sub("eta_intra_","" ,all_settings$Var4))

all_settings$dir <- sapply(1:nrow(all_settings),function(i,all_settings){
  paste0(paste0(as.character(all_settings[i,1:4]),collapse = "/"),".csv")
  
},all_settings)

all_posterior_medians <- lapply(1:nrow(all_settings), function(i,all_settings){
  temp <- read.csv(all_settings$dir[i],row.names = 1)
  medians <- apply(as.matrix(temp),1,median)
  data.frame(medians = medians,
             size = all_settings$Var2[i],
             setting = all_settings$Var3[i],
             para_full = all_settings$Var4[i],
             para = all_settings$para[i],
             spp = all_settings$spp[i])
},all_settings)

all_posterior_medians <- Reduce(rbind,all_posterior_medians)

all_posterior_medians_APIS <- all_posterior_medians[all_posterior_medians$size=="APIS",]
all_posterior_medians_other <- all_posterior_medians[all_posterior_medians$size!="APIS",]

all_posterior_medians_APIS$medians <- sapply(1:nrow(all_posterior_medians_APIS),
                                             function(i,all_posterior_medians_APIS,true_setting_APIS){
                                               all_posterior_medians_APIS$medians[i]-
                                                 true_setting_APIS[as.character(all_posterior_medians_APIS$setting[i]),as.character(all_posterior_medians_APIS$para_full[i])]
                                               
                                               
                                             },all_posterior_medians_APIS,true_setting_APIS) 

all_posterior_medians_other$medians <- sapply(1:nrow(all_posterior_medians_other),
                                             function(i,all_posterior_medians_APIS,true_setting_APIS){
                                               all_posterior_medians_APIS$medians[i]-
                                                 true_setting_APIS[as.character(all_posterior_medians_APIS$setting[i]),as.character(all_posterior_medians_APIS$para_full[i])]
                                             },all_posterior_medians_other,true_setting_noAPIS) 

all_posterior_medians <- rbind(all_posterior_medians_APIS,all_posterior_medians_other)


########### make plots #############

##### no interactions

makeaplot <- function(para,setting,all_posterior_medians){
  ggplot(all_posterior_medians[all_posterior_medians$para==para & 
                                                all_posterior_medians$setting==setting,],
                aes(size, medians, fill=factor(spp))) + 
          geom_boxplot(width=0.6) + 
          ylim(-1,1)+
          xlab("")+
          ylab("")+
          theme(legend.position="none") + 
          scale_fill_brewer()+
          theme(text = element_text(size=14), 
                axis.text.x = element_text(angle=0,size = 12),
                plot.margin = margin(.15, .15, .15, .15, "cm"))+
          geom_hline(yintercept = 0) 
}

all_plots <- expand.grid(as.character(unique(all_posterior_medians$setting)),
                         as.character(unique(all_posterior_medians$para))[c(3,4,1,2)]
                         )

all_ggplots <- lapply(1:nrow(all_plots),function(i,all_plots,all_posterior_medians){
  makeaplot(all_plots[i,2],all_plots[i,1],all_posterior_medians)
},all_plots,all_posterior_medians)


library(grid)

cowplot::plot_grid(plotlist=all_ggplots, ncol=3, align='vh')
ggsave(filename = "bigsimulation_temp.png",width = 9,height = 8,unit="in",dpi=500)


