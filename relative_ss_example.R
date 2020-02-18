#------------------------------#
#--- Relative Sink-Strength ---#
#------------------------------#

#--- Assumptions:
#---  For the sake of this example we assume our crop is constituted by 3 Organs: 
#---  1 Root, 1 Stem and 1 Leaf

#--- Load libraries
library(ggplot2)
library(reshape2)

#--- Flexible sigmoid function
flex.sigm.fun = function(d,y.ini,y.end,tb,tm,te,t,delta){
  
  #--- A flexible sigmoid function of determinate growth.
  #--- https://doi.org/10.1093/aob/mcg029
  
  #--- Few modification were done:
  #---  (i)     y.ini and y.end were included in place of wmax to account for initial conditions
  #---  (ii)    option for cm = 1, for relative rate curve
  #---  (iii)   parameters and time is now scaled to cope with different initial tb and te
  
  
  #--- Input parameters
  # integer d               # Derivative type (Integral form [d=0], 1st derivative [d=1], 1st derivative but cm is set to 1 to return relative response [d=2])
  # real    y.ini           # Initial value,  when t = tb: y = y.ini
  # real    y.end           # Final value,    when t = te: y = y.end
  # real    tb              # Initial time of growth
  # real    tm              # Time of maximun growth rate
  # real    te              # Time for growth end
  # real    t               # Current time
  # real    delta           # Asymmetry shape factor
  
  #--- Variables
  # real    cm              # Maximum Growth Rate (inflexion point) [dw/dt]
  
  #--- Outputs
  # real    y               # Current Biomass in time t
  # real    dy              # Current Growth Rate in time t [dy/dt]
  
  fnm = "sigm.fun"
  
  if(missing(d)){d = 0}         # Integral Form as default
  if(missing(delta)){delta = 1} # No Asymmetry as default
  if(missing(y.ini) & d == 2){y.ini = 0}
  if(missing(y.end) & d == 2){y.end = 1}
  
  #--- parameter correction for growing period
  tb.cor  = 0
  tm.cor  = tm - tb
  te.cor  = te - tb
  
  #--- t correction for growing period
  t.cor             = t
  t.cor[t.cor<tb]   = 0
  t.cor[t.cor>te]   = 1
  t.cor[t.cor >= tb & t.cor <= te] = 1 - ((t.cor[t.cor >= tb & t.cor <= te] - te) / (tb - te))
  t.cor = t.cor * te.cor
  
  if(d == 0){
    
    #--- Integral Form
    y   = y.ini + (y.end - y.ini) * (1. + (te.cor-t.cor)/(te.cor-tm.cor)) * (t.cor/te.cor)^(te.cor/(te.cor-tm.cor))
    
    #--- Constraint to wmax (either forvector, single value or gridded data)
    if(length(y) == length(y.end)){
      y[t>=te] = y.end[t>=te]  
    }else{
      y[t>=te] = y.end  
    }
    
    return(y)
  }
  
  if(d == 1){
    
    #--- 1st Derivative Form
    cm  = ((2. * te.cor - tm.cor) / (te.cor * (te.cor-tm.cor)) * (tm.cor / te.cor) ^ (tm.cor/(te.cor-tm.cor))) * (y.end - y.ini)
    dy  = cm * (((te.cor-t.cor)/(te.cor-tm.cor)) * ((t.cor-tb.cor)/(tm.cor-tb.cor)) ^ ((tm.cor-tb.cor)/(te.cor-tm.cor))) ^ delta
    
    #--- constrait to time limits
    dy[t >= te] = 0
    dy[t <  tb] = 0
    
    return(dy)
  }
  
  if(d == 2){
    
    #--- 1st Derivative Form but cm = 1 to return relative rate response
    cm  = 1
    ry  = cm * (((te.cor-t.cor)/(te.cor-tm.cor)) * ((t.cor-tb.cor)/(tm.cor-tb.cor)) ^ ((tm.cor-tb.cor)/(te.cor-tm.cor))) ^ delta
    
    #--- constrait to time limits
    ry[t >= te] = 0
    ry[t <  tb] = 0
    
    return(ry)
    
  }
  
}

#--- Growth Parameters Under Optimun Conditions
root.ini.w = 0
root.end.w = 10
root.ini.tt= 0
root.mid.tt= 50
root.end.tt= 300

stem.ini.w = 0
stem.end.w = 25
stem.ini.tt= 150
stem.mid.tt= 400
stem.end.tt= 600

leaf.ini.w = 0
leaf.end.w = 20
leaf.ini.tt= 20
leaf.mid.tt= 200
leaf.end.tt= 600

#----------------------------------------#
#--- Simulation Setup No Water-stress ---#
#----------------------------------------#

n.days    = 200
sim.tt    = rnorm(n.days, 7,2) # assuming an average thermal time of 7 Cdays and 2 Cdays deviation

#--- Sensitivity to Water Stress
leaf.sens = 1.0 # 100%  sensitive
root.sens = 0.0 # 0%    sensitive
stem.sens = 0.2 # 20%   sensitive

h2o.stress= rep(1,n.days) # No Water Stress

#--- Daily Steps
for(day in 1:n.days){
  
  if(day == 1){
    #--- Initialize tt
    cum.tt = 0
  }
  
  #--- Thermal time today
  today.tt = max(0,sim.tt[day])
  
  #--- Calculate organs sink strength for today in g organ-1
  root.ss = flex.sigm.fun(1,
                          root.ini.w,
                          root.end.w,
                          root.ini.tt,
                          root.mid.tt,
                          root.end.tt,
                          cum.tt,
                          delta = 1) * today.tt
  
  stem.ss = flex.sigm.fun(1,
                          stem.ini.w,
                          stem.end.w,
                          stem.ini.tt,
                          stem.mid.tt,
                          stem.end.tt,
                          cum.tt,
                          delta = 1) * today.tt
  
  leaf.ss = flex.sigm.fun(1,
                          leaf.ini.w,
                          leaf.end.w,
                          leaf.ini.tt,
                          leaf.mid.tt,
                          leaf.end.tt,
                          cum.tt,
                          delta = 1) * today.tt
  
  #--- Include Water Stress
  root.ss = root.ss * (1 - ((1 - h2o.stress[day]) * root.sens))
  leaf.ss = leaf.ss * (1 - ((1 - h2o.stress[day]) * leaf.sens))
  stem.ss = stem.ss * (1 - ((1 - h2o.stress[day]) * stem.sens))
  
  #--- Total Crop Sink-Strength for today
  crop.ss = root.ss + stem.ss + leaf.ss
  
  #--- Compute partitioning based on the relative sink-strength
  if(crop.ss > 0 ){
    pfac.root = root.ss / crop.ss
    pfac.stem = stem.ss / crop.ss
    pfac.leaf = leaf.ss / crop.ss
  }else{
    pfac.root = 0
    pfac.stem = 0
    pfac.leaf = 0
  }
  
  #--- Update tt
  cum.tt = cum.tt + today.tt
  
  #--- Store states
  if(day == 1){
    
    df.out = data.frame(tt         = cum.tt,
                                     root.ss    = root.ss,     
                                     stem.ss    = stem.ss,
                                     leaf.ss    = leaf.ss,
                                     pfac.root  = pfac.root,
                                     pfac.stem  = pfac.stem,
                                     pfac.leaf  = pfac.leaf)
    
  }else{
    df.out = rbind(df.out,data.frame(tt         = cum.tt,
                                     root.ss    = root.ss,     
                                     stem.ss    = stem.ss,
                                     leaf.ss    = leaf.ss,
                                     pfac.root  = pfac.root,
                                     pfac.stem  = pfac.stem,
                                     pfac.leaf  = pfac.leaf))
  }
    
}

df.pfac = melt(df.out[,c("tt","pfac.root","pfac.stem","pfac.leaf")],  id.vars = "tt")
df.ss   = melt(df.out[,c("tt","root.ss","stem.ss","leaf.ss")],  id.vars = "tt")

#--- plot pfacs
gg.pfac.optimum.h2o = 
ggplot(df.pfac,aes(x = tt, y = value, colour  = variable)) +
  geom_line() + 
  ylab("Partitioning Factors") +
  xlab("Thermal time") +
  xlim(0,max(c(leaf.end.tt,
               root.end.tt,
               stem.end.tt))) + 
  theme_bw()

#--- plot ss (g)
gg.ss.optimum.h2o =
ggplot(df.ss,aes(x = tt, y = value, colour  = variable)) +
  geom_line() + 
  xlim(0,max(c(leaf.end.tt,
               root.end.tt,
               stem.end.tt))) + 
  theme_bw()


#------------------------------------------#
#--- Simulation Setup With Water-stress ---#
#------------------------------------------#

n.days    = 200
sim.tt    = rnorm(n.days, 7,2) # assuming an average thermal time of 7 Cdays and 2 Cdays deviation

#--- Sensitivity to Water Stress
leaf.sens = 1.0 # 100%  sensitive
root.sens = 0.0 # 0%    sensitive
stem.sens = 0.2 # 20%   sensitive

h2o.stress= rnorm(n.days, 0.5,0.1)  # Average of 0.5 dev 0.1
h2o.stress[h2o.stress < 0 ] = 0     # Contraint 0-1
h2o.stress[h2o.stress > 1 ] = 1     # Contraint 0-1

#--- Daily Steps
for(day in 1:n.days){
  
  if(day == 1){
    #--- Initialize tt
    cum.tt = 0
  }
  
  #--- Thermal time today
  today.tt = max(0,sim.tt[day])
  
  #--- Calculate organs sink strength for today in g organ-1
  root.ss = flex.sigm.fun(1,
                          root.ini.w,
                          root.end.w,
                          root.ini.tt,
                          root.mid.tt,
                          root.end.tt,
                          cum.tt,
                          delta = 1) * today.tt
  
  stem.ss = flex.sigm.fun(1,
                          stem.ini.w,
                          stem.end.w,
                          stem.ini.tt,
                          stem.mid.tt,
                          stem.end.tt,
                          cum.tt,
                          delta = 1) * today.tt
  
  leaf.ss = flex.sigm.fun(1,
                          leaf.ini.w,
                          leaf.end.w,
                          leaf.ini.tt,
                          leaf.mid.tt,
                          leaf.end.tt,
                          cum.tt,
                          delta = 1) * today.tt
  
  #--- Include Water Stress
  root.ss = root.ss * (1 - ((1 - h2o.stress[day]) * root.sens))
  leaf.ss = leaf.ss * (1 - ((1 - h2o.stress[day]) * leaf.sens))
  stem.ss = stem.ss * (1 - ((1 - h2o.stress[day]) * stem.sens))
  
  #--- Total Crop Sink-Strength for today
  crop.ss = root.ss + stem.ss + leaf.ss
  
  #--- Compute partitioning based on the relative sink-strength
  if(crop.ss > 0 ){
    pfac.root = root.ss / crop.ss
    pfac.stem = stem.ss / crop.ss
    pfac.leaf = leaf.ss / crop.ss
  }else{
    pfac.root = 0
    pfac.stem = 0
    pfac.leaf = 0
  }
  
  #--- Update tt
  cum.tt = cum.tt + today.tt
  
  #--- Store states
  if(day == 1){
    
    df.out = data.frame(tt         = cum.tt,
                        root.ss    = root.ss,     
                        stem.ss    = stem.ss,
                        leaf.ss    = leaf.ss,
                        pfac.root  = pfac.root,
                        pfac.stem  = pfac.stem,
                        pfac.leaf  = pfac.leaf)
    
  }else{
    df.out = rbind(df.out,data.frame(tt         = cum.tt,
                                     root.ss    = root.ss,     
                                     stem.ss    = stem.ss,
                                     leaf.ss    = leaf.ss,
                                     pfac.root  = pfac.root,
                                     pfac.stem  = pfac.stem,
                                     pfac.leaf  = pfac.leaf))
  }
  
}

df.pfac = melt(df.out[,c("tt","pfac.root","pfac.stem","pfac.leaf")],  id.vars = "tt")
df.ss   = melt(df.out[,c("tt","root.ss","stem.ss","leaf.ss")],  id.vars = "tt")

#--- plot pfacs
gg.pfac.h2o.stress = 
ggplot(df.pfac,aes(x = tt, y = value, colour  = variable)) +
  geom_line() + 
  ylab("Partitioning Factors") +
  xlab("Thermal time") +
  xlim(0,max(c(leaf.end.tt,
               root.end.tt,
               stem.end.tt))) + 
  theme_bw()

#--- plot ss (g)
gg.ss.h2o.stress = 
ggplot(df.ss,aes(x = tt, y = value, colour  = variable)) +
  geom_line() + 
  xlim(0,max(c(leaf.end.tt,
               root.end.tt,
               stem.end.tt))) + 
  theme_bw()

#--- Check Comparison No-stress vs h2o.stress
gg.pfac.optimum.h2o
gg.pfac.h2o.stress

#--- Save ggs
ggsave("pfac_optimum_h2o.png",
       gg.pfac.optimum.h2o,
       width  = 5,
       height = 4,
       dpi    = 300)

ggsave("pfac_h2o_stress.png",
       gg.pfac.h2o.stress,
       width  = 5,
       height = 4,
       dpi    = 300)

