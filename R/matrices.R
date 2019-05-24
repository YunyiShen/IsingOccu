nspp = 2
ncov = ncol(envX)
hs = matrix(0,nrow(envX),nspp)

for(i in 1:nspp){
  beta_temp = kk$means$beta_occu[1:ncov + (i-1)*ncov]
  hs[,i] = envX %*% beta_temp
  
}

library(RColorBrewer)
library(export)
# E = matrix(0,nspp,nspp)
# for(i in 1:nspp){
#   for(j in 1:nspp){
#     if(i==j) next
#     E[i,j] = mean(hs[,i]*hs[,j])
#     
#   }
#   
# }
E = cov(hs)
diag(E)=0
cuts=0.2*(-100:100)/200
pal <- colorRampPalette(c("pink","white","darkgreen"))
raster::plot(raster::raster(E),col = pal(200))
graph2ppt(file = "TJHcarn_envmat.pptx")

raster::plot(raster::raster(matrix(kk$means$spp_mat,nspp,nspp)),col = pal(200))
graph2ppt(file = "FM_intermat.pptx")

raster::plot(raster::raster(matrix(kk$means$spp_mat_det,nspp,nspp)),col =pal(200))
graph2ppt(file = "FM_behamat.pptx")
