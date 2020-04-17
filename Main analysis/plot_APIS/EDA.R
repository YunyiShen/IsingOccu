require(ggplot2)
require(ggmap)

require(ggplot2)
require(reshap2)
require(grid)
require(dplyr)

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

set.seed(42)

## prepare the data

Z_abs = read.csv("./data/APIS/PA_all_full.csv",stringsAsFactors = F)
dist = read.csv("./data/APIS/dist_to_mainland.csv",stringsAsFactors = F)
Posi = read.csv("./data/APIS/CT_posi_only_island.csv")
camera_days = read.csv("./data/APIS/cameradays.csv")

islands_list = unique(substr(Z_abs$X,1,nchar(Z_abs$X)-2))

Z_abs_01 = Z_abs
Z_abs_01$X = substr(Z_abs$X,1,nchar(Z_abs$X)-2)
Z_abs_01[,2:30]=(Z_abs_01[,2:30]+1)/2

Z_abs_01$FM = 1-(1-Z_abs_01$Fisher)*(1-Z_abs_01$Marten)
Z_abs_01$CF = 1-(1-Z_abs_01$Coyote)*(1-Z_abs_01$Fox_red)
Z_abs_01$dist = dist$dist

names(Z_abs_01[,1])="site"
Z_abs_01$site = Z_abs_01$X

Z_abs_01$X = Posi$X
Z_abs_01$Y = Posi$Y

Z_abs_01$lat = Posi$Lat
Z_abs_01$long = Posi$Long

Z_abs_01$cameradays = camera_days$V1


## check the map
map_data = Z_abs_01[,c("long","lat","Fisher","Marten","FM","Coyote","Fox_red","CF","cameradays")]
map_data[,3:8] = lapply(map_data[,3:8],as.factor)


APIS_map = get_stamenmap(bbox = c(
  left = -91.05, bottom = 46.75, 
  right =-90.35, top = 47.10),zoom = 12)

Fi_map = ggmap(APIS_map) + 
  geom_point(aes(x=long,y=lat,color= Fisher),data = map_data,size = 1.2)+
  theme(text = element_text(size=15))

#ggsave("Fisher_map.png",dpi=500)

Ma_map = ggmap(APIS_map) + 
  geom_point(aes(x=long,y=lat,color=  Marten),data = map_data,size = 1.2)+
  theme(text = element_text(size=15))
#ggsave("Marten_map.png",dpi=500)

FM_map = ggmap(APIS_map) + 
  geom_point(aes(x=long,y=lat,color=  FM),data = map_data,size = 1.2)+
  theme(text = element_text(size=15))
#ggsave("FM_together_map.png",dpi=500)


Co_map = ggmap(APIS_map) + 
  geom_point(aes(x=long,y=lat,color=  Coyote),data = map_data,size = 1.2)+
  theme(text = element_text(size=15))
#ggsave("Coyote_map.png",dpi=500)

Fo_map = ggmap(APIS_map) + 
  geom_point(aes(x=long,y=lat,color=  Fox_red),data = map_data,size = 1.2)+
  theme(text = element_text(size=15))
#ggsave("Fox_map.png",dpi=500)

CF_map = ggmap(APIS_map) + 
  geom_point(aes(x=long,y=lat,color= CF),data = map_data,size = 1.2)+
  theme(text = element_text(size=15))

camera_map = ggmap(APIS_map) + 
  geom_point(aes(x=long,y=lat,color= cameradays),data = map_data,size = 1.2)+
  theme(text = element_text(size=15))

ggsave("camera_map.png",dpi=500)

png("maps.png",height=1800,width=4500,res = 300)


grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 3)))
print(Fi_map, vp = vplayout(1, 1))
print(Ma_map, vp = vplayout(1, 2))
print(FM_map, vp = vplayout(1, 3))

print(Co_map, vp = vplayout(2, 1))
print(Fo_map, vp = vplayout(2, 2))
print(CF_map, vp = vplayout(2, 3))

dev.off()


png("together.png",height=900,width=3000,res = 300)


grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(CF_map, vp = vplayout(1, 1))
print(FM_map, vp = vplayout(1, 2))


dev.off()

## naive occupancy map

naive_occu_rate = sapply(islands_list,function(island,Z){
  temp = Z[Z$X==island,2:32]
  colSums(temp)/nrow(temp)
  
},Z_abs_01)

naive_occu_rate = t(naive_occu_rate)


dist_island = dist
dist_island$sites = substr(dist_island$sites,1,nchar(dist_island$sites)-2)

dist_island_min = sapply(islands_list,function(island,dist){
  min( dist[dist$site==island,2])
  
  
},dist_island)

naive_occu_rate = as.data.frame(naive_occu_rate)
naive_occu_rate$dist = dist_island_min

plot(dist_island_min,naive_occu_rate$Marten)
plot(dist_island_min,naive_occu_rate$Fisher)

plot(dist_island_min,naive_occu_rate[,"fisher_marten"])
plot(dist_island_min,naive_occu_rate[,"coyote_fox"])

### simple logitstics regression on distance
logistic_fm = glm(fisher_marten~exp(-dist/max(dist)),data = naive_occu_rate,family = binomial(link = "logit"))
logistic_cf = glm(coyote_fox~exp(-2*(dist-min(dist))/(max(dist)-min(dist))),data = naive_occu_rate,family = binomial(link = "logit"))

logistic_site_fm = glm(fisher_marten~exp(-2*(dist-min(dist))/(max(dist)-min(dist))),data = Z_abs_01,family = binomial)
logistic_site_fi = glm(Fisher~exp(-2*(dist-min(dist))/(max(dist)-min(dist))),data = Z_abs_01,family = binomial)
logistic_site_ma = glm(Marten~exp(-2*(dist-min(dist))/(max(dist)-min(dist))),data = Z_abs_01,family = binomial)

logistic_site_cf = glm(coyote_fox~exp(-2*(dist-min(dist))/(max(dist)-min(dist))),data = Z_abs_01,family = binomial)
logistic_site_fo = glm(Coyote~exp(-2*(dist-min(dist))/(max(dist)-min(dist))),data = Z_abs_01,family = binomial)
logistic_site_co = glm(Fox_red~exp(-2*(dist-min(dist))/(max(dist)-min(dist))),data = Z_abs_01,family = binomial)



png( file="dethist.png", height=6, width=8,units = "in" ,res = 600)

detmat1 = list(as.matrix(read.csv(paste0(link,"Fisher_Marten_60dfull_by_islands.csv"),header = F)))

par(mfrow=c(2,2)) 
image(x=1:17,y=1:155,z=t(detmat1[[1]][1:155,]),xlab = "repeats",ylab = "site",main = "Fisher")
image(x=1:17,y=1:155,z=t(detmat1[[1]][1:155+155,]),xlab = "repeats",ylab = "site",main = "Marten")

detmat2 = list(as.matrix(read.csv(paste0(link,"Coyote_Fox_Bobcat_90dfull_by_islands.csv"),header = F)[1:310,]))

image(x=1:17,y=1:155,z=t(detmat2[[1]][1:155,]),xlab = "repeats",ylab = "site",main = "Coyote")
image(x=1:17,y=1:155,z=t(detmat2[[1]][1:155+155,]),xlab = "repeats",ylab = "site",main = "Fox,red")


dev.off()













