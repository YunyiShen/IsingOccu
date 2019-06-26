source("pPCA.R")
link = "E:/UW Lab jobs/2. ISING Occupancy model/4. DATA/TJH/"
design = read.csv(paste0(link,"TJH_grid_env_coverted_cordi.csv")
                  ,header = T,row.names = 1
                  ,fileEncoding = "UTF-8")

pPCA_res = pPCA( design[,c(1)],design[,2:12],ML=T)

corr = cor(design)

ww = eigen(corr)

design[,1] = (design[,1]-mean(design[,1]))/sd(design[,1])

FP = prcomp(design)
