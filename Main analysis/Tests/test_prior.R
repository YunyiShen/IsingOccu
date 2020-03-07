xs = expand.grid(seq(2,-2,-0.01),seq(-2,2,.01))

ps_co = (exp(xs$Var1+xs$Var1+xs$Var2))/
  (exp(xs$Var1+xs$Var1+xs$Var2)+
     exp(-xs$Var1+xs$Var1-xs$Var2)+
     exp(xs$Var1-xs$Var1-xs$Var2)+
     exp(-xs$Var1-xs$Var1+xs$Var2))

ps_co_m = matrix(ps_co,sqrt(length(ps_co)))

raster::plot(raster::raster(ps_co_m))

ps_no = (exp(-xs$Var1-xs$Var1+xs$Var2))/
  (exp(xs$Var1+xs$Var1+xs$Var2)+
     exp(-xs$Var1+xs$Var1-xs$Var2)+
     exp(xs$Var1-xs$Var1-xs$Var2)+
     exp(-xs$Var1-xs$Var1+xs$Var2))

ps_no_m = matrix(ps_no,sqrt(length(ps_no)))

raster::plot(raster::raster(ps_no_m))
  
