source("code/parameters.R")
source("code/models/water_model.R")
source("code/W0_funcs.R")

param <- param_list_drake_W0_1000()

W0min(param)

# W0_maxS
W0_maxS <- W0max(param, simple_water_model_sim_time_breaks, c(4000,7000), deltaW0 = 100, deltat=0.1, kf=0)
print(W0_maxS)
W0_maxS_100 = W0_maxS
W0_maxS_100
W0_maxS <- W0max(param, simple_water_model_sim_time_breaks, c(W0_maxS-100,W0_maxS), deltaW0 = 10, deltat=0.1, kf=0)
print(W0_maxS)
W0_maxS_10 = W0_maxS
W0_maxS <- W0max(param, simple_water_model_sim_time_breaks, c(W0_maxS-10,W0_maxS), deltaW0 = 1, deltat=0.1, kf=0)
print(W0_maxS)

# W0_maxM+S
W0_maxMS <- W0max(param, simple_water_model_sim_time_breaks, c(4000,7000), deltat=0.1, kf=0.5)
kprint(W0_maxMS)
W0_maxMS_100 = W0_maxMS
W0_maxMS_100
W0_maxMS <- W0max(param, simple_water_model_sim_time_breaks, c(W0_maxMS-100,W0_maxMS), deltaW0 = 10, deltat=0.1, kf=0.5)
print(W0_maxMS)
W0_maxMS_10 = W0_maxMS
W0_maxMS <- W0max(param, simple_water_model_sim_time_breaks, c(W0_maxMS-10,W0_maxMS), deltaW0 = 1, deltat=0.1, kf=0.5)
print(W0_maxMS)

kf_list=seq(from=0, to=1, by=0.1)
seqlength = length(kf)
W0_max_res <- vector()

j = 1
for(kf in kf_list){
  print(kf)
  W0_max <- W0max(param, simple_water_model_sim_time_breaks, c(W0min(param),7000), deltaW0 = 100, deltat=0.1, kf=kf)
  W0_max <- W0max(param, simple_water_model_sim_time_breaks, c(W0_max-100,W0_max), deltaW0 = 10, deltat=0.1, kf=kf)
  W0_max <- W0max(param, simple_water_model_sim_time_breaks, c(W0_max-10,W0_max), deltaW0 = 1, deltat=0.1, kf=kf)
  W0_max_res[j] = W0_max
  j = j+1
}


