source("misc_island.R")

###### graph data ######
island = read.csv("C:/Users/yshen99/Documents/GitHub/RFIBM_MCMC/Island/CT_posi_only_island.csv")

link_inner = as.matrix( read.csv("C:/Users/yshen99/Documents/GitHub/RFIBM_MCMC/Island/link_inner.csv",row.names = 1))
#link_outer = as.matrix( read.csv("link_outer.csv",row.names = 1))
link_outer = as.matrix( read.csv("C:/Users/yshen99/Documents/GitHub/RFIBM_MCMC/Island/link_outer_full.csv",row.names = 1))
link_mainland = as.matrix( read.csv("C:/Users/yshen99/Documents/GitHub/RFIBM_MCMC/Island/link_mainland.csv"))
#link_outer = 0 * link_outer # this makes it a mainland-island system
#link_mainland = matrix(0,155,1)

distM_full = as.matrix( read.csv("C:/Users/yshen99/Documents/GitHub/RFIBM_MCMC/Island/distM_full.csv",row.names = 1))
distM_mainland = as.matrix( read.csv("C:/Users/yshen99/Documents/GitHub/RFIBM_MCMC/Island/dist_to_mainland.csv",row.names = 1))

intcd = min(min((distM_mainland*link_mainland)[(distM_mainland*link_mainland)>0]),
            min((link_outer*distM_full)[(link_outer*distM_full)>0]))
normd = max(max(distM_mainland*link_mainland),max(link_outer*distM_full))-intcd


distM_full = (distM_full-intcd)/normd # normalizing the distance
distM_mainland = (distM_mainland-intcd)/normd

###### simulation ######

spp_mat = matrix(c(0,1,1,0),2,2)
envX = matrix(1,155,1)

theta = list(beta_occu = c(0,0),
             eta_intra = c(.15,.15),
             eta_inter = c(.3,.3),
             d_inter = c(.2,.2),
             spp_mat = -0.1 * spp_mat)

link_map = 
  list(inter = link_outer,
       intra = link_inner)

set.seed(42)
Z_sample = rIsingOccu_multi(theta,
                            envX,
                            distM_full,link_map ,
                            distM_mainland,link_mainland,
                            int_range_intra="nn",int_range_inter="exp",
                            n=4,method = "CFTP",nIter = 100)

require(ggplot2)

tempdata = data.frame(island[,6:7],
                      Z_1 = Z_sample[1:155,1],
                      Z_2 = Z_sample[156:310,1])

ggplot(data = tempdata,aes(x=X,y=Y,color = Z_1))+
  geom_point()

Hamiltonian(theta,envX,distM_full,link_map,distM_mainland,link_mainland,int_range_intra="nn",int_range_inter="exp",Z_sample)

