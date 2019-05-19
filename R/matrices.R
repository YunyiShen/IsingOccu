nspp = 3
ncov = ncol(envX)

for(i in 1:nspp){
  beta_temp = kk$means$beta_occu[1:ncov + (i-1)*ncov]
  hs[,i] = envX %*% beta_temp
  
}

library(RColorBrewer)

# E = matrix(0,nspp,nspp)
# for(i in 1:nspp){
#   for(j in 1:nspp){
#     if(i==j) next
#     E[i,j] = mean(hs[,i]*hs[,j])
#     
#   }
#   
# }

cuts=0.2*(-100:100)/200
pal <- colorRampPalette(c("pink","white","darkgreen"))
raster::plot(raster::raster(cov(hs)),col = pal(200))
graph2ppt(file = "TJHcarn_envmat.pptx")

raster::plot(raster::raster(matrix(kk$means$spp_mat,3,3)),col = pal(200))
graph2ppt(file = "TJHcarn_intermat.pptx")

raster::plot(raster::raster(matrix(kk$means$spp_mat_det,3,3)),col =pal(200))
graph2ppt(file = "TJHcarn_behamat.pptx")
