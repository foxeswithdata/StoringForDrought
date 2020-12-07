water_analysis_figures_w0 <- function(strategy_df, directory){
  
  if (!file.exists(directory)){
    dir.create(file.path(directory))
  }
  
  p <- ggplot(strategy_df, aes(x=w0)) +
    geom_line(aes(y=ts)) +
    xlab("Total water input [kgH20]") + 
    ylab("Optimal growth switch day ts")
  
  ggsave( "optimal_switch_t.png", plot = p, device = png(), path = directory, dpi = 400)
  
  p <- ggplot(strategy_df, aes(x=w0)) +
    geom_line(aes(y=mT, colour="MT")) +
    geom_line(aes(y=sT, colour="ST")) +
    xlab("Total water input [kgH20]") + 
    ylab("Final biomass")+ 
    scale_colour_manual("", 
                            values = c("MT"="#32A287", "ST"="#6C464E"))
  
  ggsave( "optimal_final_biomass_storage.png", plot = p, device = png(), path = directory, dpi = 400)
  
  
  p <- ggplot(strategy_df, aes(x=w0)) +
    geom_line(aes(y=wT)) +
    xlab("Total water input [kgH20]") + 
    ylab("Final water availability [kgH20]")
  
  ggsave( "optimal_final_water_available.png", plot = p, device = png(), path = directory, dpi = 400)
  
  p <- ggplot(strategy_df, aes(x=w0)) +
    geom_line(aes(y=a)) +
    xlab("Total water input [kgH20]") + 
    ylab("Total photosynthesis gCO2")
  
  ggsave( "optimal_total_A.png", plot = p, device = png(), path = directory, dpi = 400)
  
  
}

water_analysis_figures_tend <- function(strategy_df, directory){
  
  if (!file.exists(directory)){
    dir.create(file.path(directory))
  }
  
  p <- ggplot(strategy_df, aes(x=tend)) +
    geom_line(aes(y=ts)) +
    xlab("Length of simulation (days)") + 
    ylab("Optimal growth switch day ts")
  
  ggsave( "optimal_switch_t.png", plot = p, device = png(), path = directory, dpi = 400)
  
  p <- ggplot(strategy_df, aes(x=tend)) +
    geom_line(aes(y=mT, colour="MT")) +
    geom_line(aes(y=sT, colour="ST")) +
    xlab("Length of simulation (days)") + 
    ylab("Final biomass")+ 
    scale_colour_manual("", 
                        values = c("MT"="#32A287", "ST"="#6C464E"))
  
  ggsave( "optimal_final_biomass_storage.png", plot = p, device = png(), path = directory, dpi = 400)
  
  
  p <- ggplot(strategy_df, aes(x=tend)) +
    geom_line(aes(y=wT)) +
    xlab("Length of simulation (days)") + 
    ylab("Final water availability [kgH20]")
  
  ggsave( "optimal_final_water_available.png", plot = p, device = png(), path = directory, dpi = 400)
  
  p <- ggplot(strategy_df, aes(x=tend)) +
    geom_line(aes(y=a)) +
    xlab("Length of simulation (days)") + 
    ylab("Total photosynthesis gCO2")
  
  ggsave( "optimal_total_A.png", plot = p, device = png(), path = directory, dpi = 400)
  
  
}

water_analysis_figures_kp <- function(strategy_df, directory){
  
  if (!file.exists(directory)){
    dir.create(file.path(directory))
  }
  
  p <- ggplot(strategy_df, aes(x=kp)) +
    geom_line(aes(y=ts)) +
    xlab("photosynthesis parameter kp [gCg^-1Cday^-1]") + 
    ylab("Optimal growth switch day ts")
  
  ggsave( "optimal_switch_t.png", plot = p, device = png(), path = directory, dpi = 400)
  
  p <- ggplot(strategy_df, aes(x=kp)) +
    geom_line(aes(y=mT, colour="MT")) +
    geom_line(aes(y=sT, colour="ST")) +
    xlab("photosynthesis parameter kp [gCg^-1Cday^-1]") + 
    ylab("Final biomass")+ 
    scale_colour_manual("", 
                        values = c("MT"="#32A287", "ST"="#6C464E"))
  
  ggsave( "optimal_final_biomass_storage.png", plot = p, device = png(), path = directory, dpi = 400)
  
  
  p <- ggplot(strategy_df, aes(x=kp)) +
    geom_line(aes(y=wT)) +
    xlab("photosynthesis parameter kp [gCg^-1Cday^-1]") + 
    ylab("Final water availability [kgH20]")
  
  ggsave( "optimal_final_water_available.png", plot = p, device = png(), path = directory, dpi = 400)
  
  p <- ggplot(strategy_df, aes(x=kp)) +
    geom_line(aes(y=a)) +
    xlab("photosynthesis parameter kp [gCg^-1Cday^-1]") + 
    ylab("Total photosynthesis gCO2")
  
  ggsave( "optimal_total_A.png", plot = p, device = png(), path = directory, dpi = 400)
  
  
}
