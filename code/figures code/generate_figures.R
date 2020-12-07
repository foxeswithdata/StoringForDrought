# 1) Figures for Storage and Biomass

plot_biomass_storage <- function(M,S){
  par(mfrow=c(1,2))
  plot(M, type='l', main='Change in biomass', xlab="Time (days)", ylab="Biomass (kg)")
  plot(S, type='l', main='Change in storage', xlab="Time (days)", ylab="Storage (kg)")
}


# 2) Figure for adjunct values 

# 3) Varying parameters (initial size of S and M, proportion of allocation) vs result on the ts


# 4) Results on S and M of varying ts manually
