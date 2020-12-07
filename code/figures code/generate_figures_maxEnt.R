library(tidyverse)
library(plyr)
library(dplyr)
library(ggpubr)

pk_dist <- lapply(pk_dist_list, function(x) return(l2df(x, maxts)))
pk_df <- ldply(pk_dist, rbind)
pk_df$Itot <- as.factor(pk_df$Itot)

p <- ggplot(pk_df, aes(x=ts, y=pk_distribiution, group=Itot)) +
  geom_vline(xintercept=30, color="#AAAAAA") +
  geom_line(aes(linetype=Itot, color = Itot)) + 
  xlab("growth switch day ts (days)") + 
  ylab("probability of growth switch day ts") + 
  ggtitle("Distribution of switch day probability \n for different Itot values (only live strategies)")
p

q <- ggplot(pk_df, aes(x=Ek, y=pk_distribiution, group=Itot)) +
  geom_line(aes(linetype=Itot, color = Itot)) + 
  xlab("Final total evapotranspiration") + 
  ylab("probability of equivalent growth switch day ts") +
  ggtitle("Relationship between probability of switch day \n and evapotranspiration \n for different Itot values (only live strategies)")
q

detach(package:plyr)
out <- pk_df %>% group_by(Ek,Itot) %>%
  summarise(pk_distribiution=sum(pk_distribiution))

t <- ggplot(out, aes(x=Ek, y=log(pk_distribiution), group=Itot)) +
  geom_line(aes(linetype=Itot, color = Itot)) + 
  xlab("Total evapotranspiration [units]") + 
  ylab("log-probability \n of equivalent growth switch day ts") +
  ggtitle("Relationship between probability of switch day (log transformed) \n and evapotranspiration \n for different Itot values (all strategies)")
t

figure <- ggarrange(q,t,
                    ncol = 2, nrow = 1)
figure

# Plot entropy and Itot 

library(plyr)
pk_dist_long <- lapply(pk_dist_list_long, function(x) return(l2df(x, maxts)))
pk_df_long <- ldply(pk_dist_long, rbind)

unique(pk_df_long$Itot)

r <- ggplot(pk_df_long, aes(x=Itot, y=entropy)) +
  geom_line() + 
  xlab("Total evapotranspiration per plant Itot [units]") + 
  ylab("Total entropy")+
  ggtitle("Relationship between entropy and Itot \n (only live strategies)")

r

o <- ggplot(pk_df_long, aes(x=Itot, y=variance)) +
  geom_line() + 
  xlab("Total evapotranspiration per plant Itot [units]") + 
  ylab("Variance of pk") + 
  ggtitle("Relationship between variance and Itot (only live strategies)")
o

s <- ggplot(pk_df_long, aes(x=variance, y=entropy)) +
  geom_point() + 
  xlab("variance of pk") + 
  ylab("Total entropy") + 
  ggtitle("Relationship between variance and entropy \n (only live strategies)")
s

#### NOW PLOT THE DEAD TREES ####
library("plyr")
pk_dist <- lapply(pk_dist_list_dead, function(x) return(l2df(x, T_end)))
pk_df <- ldply(pk_dist, rbind)


p <- ggplot(pk_df, aes(x=ts, y=pk_distribiution, group=Itot)) +
  geom_vline(xintercept=maxts, color="#AAAAAA") +
  geom_vline(xintercept=(tcrit_av), linetype="dashed", color="#FF0000") +
  geom_line(aes(linetype=Itot, color = Itot)) + 
  xlab("growth switch day ts (days)") + 
  ylab("probability of growth switch day ts")+ 
  ggtitle("Distribution of switch day probability \n for different Itot values (all strategies)")
p

#needs fixing for final Evap#
detach(package:plyr)
out <- pk_df %>% group_by(Ek,Itot) %>%
  summarise(pk_distribiution=sum(pk_distribiution))

out <- out[out$Itot != 1000, 1:3]
out$Itot <- as.factor(out$Itot)
out <- out[out$Ek != max(out$Ek),1:3]

q <- ggplot(out, aes(x=Ek, y=pk_distribiution, group=Itot)) +
  geom_line(aes(linetype=Itot, color = Itot)) + 
  xlab("Total evapotranspiration [units]") + 
  ylab("probability of equivalent growth switch day ts") +
  ggtitle("Relationship between probability of switch day \n and evapotranspiration \n for different Itot values (all strategies)")
q

q <- ggplot(out, aes(x=Ek, y=log(pk_distribiution), group=Itot)) +
  geom_line(aes(linetype=Itot, color = Itot)) + 
  xlab("Total evapotranspiration [units]") + 
  ylab("log-probability \n of equivalent growth switch day ts") +
  ggtitle("Relationship between probability of switch day (log transformed) \n and evapotranspiration \n for different Itot values (all strategies)")
q


# Plot entropy and Itot 

library(plyr)
pk_dist_long <- lapply(pk_dist_list_dead_long, function(x) return(l2df(x, T_end)))
pk_df_long <- ldply(pk_dist_long, rbind)

unique(pk_df_long$Itot)

r <- ggplot(pk_df_long, aes(x=Itot, y=entropy)) +
  geom_line() + 
  xlab("Total evapotranspiration per plant Itot [units]") + 
  ylab("Total entropy") +
  ggtitle("Relationship between entropy and Itot (all strategies)")
r

o <- ggplot(pk_df_long, aes(x=Itot, y=variance)) +
  geom_line() + 
  xlab("Total evapotranspiration per plant Itot [units]") + 
  ylab("Variance of pk") + 
  ggtitle("Relationship between variance and Itot (all strategies)")
o

s <- ggplot(pk_df_long, aes(x=variance, y=entropy)) +
  geom_point() + 
  xlab("variance of pk") + 
  ylab("Total entropy") + 
  ggtitle("Relationship between variance and entropy \n (only live strategies)")
s

