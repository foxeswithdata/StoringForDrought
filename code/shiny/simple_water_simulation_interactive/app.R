#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

library(tidyverse)
library(plyr)
library(dplyr)
library(ggpubr)

source("helper/water_model.R")

    
    ui = fluidPage(
        
        titlePanel("Optimal Allocation Trajectory"),
        
        sidebarLayout(
            inputPanel(
                selectInput("goal_f", label = "Goal Function:",
                            choices = c("Maximise Biomass" = "biom", "Maximise Storage" = "stor", "Maximise Storage and Biomass" = "biost"), selected = "biom"),
                selectInput("delta_t", label = "Search Tolerance:",
                            choices = c("E0" = "1", "E-1" = "0.1", "E-2" = "0.01", "E-3" = "0.001"), selected = "0.1"),
            # ),
            
            # inputPanel(
                sliderInput("W0_adjust", label = "Initial Water Availability:",
                            min = 10, max = 500, value = 200, step = 10),
                
                sliderInput("M0_adjust", label = "Initial Biomass:",
                            min = 10, max = 300, value = 100, step = 10),
                
                sliderInput("S0_adjust", label = "Initial Storage:",
                            min = 0, max = 200, value = 10, step = 10),
                
                sliderInput("Tend_adjust", label = "Total run time:",
                            min = 32, max = 128, value = 64, step = 1),
            # ),
            
            # inputPanel(
                sliderInput("kp_adjust", label = "Photosynthesis Parameter:",
                            min = 0.05, max = 1, value = 0.04, step = 0.01),
                
                sliderInput("kr_adjust", label = "Respiration Parameter:",
                            min = 0.01, max = 0.5, value = 0.02, step = 0.01),
                
                sliderInput("ks_adjust", label = "Storage Parameter:",
                            min = 0, max = 1, value = 0.1, step = 0.02),
                
                sliderInput("kw_adjust", label = "Water Use Parameter:",
                            min = 0.1, max = 3, value = 1, step = 0.1)
            ),
            
            # Show a plot of the generated distribution
            mainPanel(
                plotOutput("plot")
            ), 
            position = "left"
        )
        
    )
    
server <- function(input, output) {
        output$plot = renderPlot({
            Tend <- as.numeric(input$Tend_adjust)
            W0 <- as.numeric(input$W0_adjust)
            S0 <- as.numeric(input$S0_adjust)
            M0 <- as.numeric(input$M0_adjust)
            kp <- as.numeric(input$kp_adjust)
            kw <- as.numeric(input$kw_adjust)
            ks <- as.numeric(input$ks_adjust)
            kr <- as.numeric(input$kr_adjust)
            deltat = 0.5
            if(input$delta_t == "0.1"){
                deltat = 0.05
            }
            else if(input$delta_t == "0.01"){
                deltat = 0.005
            }
            else if(input$delta_t == "0.001"){
                deltat = 0.0005
            }
            simlength = Tend/deltat + 1
            S_out <- vector(length=simlength)
            M_out <- vector(length=simlength)
            A_out <- vector(length=simlength)
            W_out <- vector(length=simlength)
            withProgress(message = 'Making plot', value = 0, {
                k = 1
                for (i in seq(from=0, to=Tend, by=deltat)){
                    out <- simple_water_model_sim_time_breaks(M0,S0,W0, Tend, kp, kr, ks, kw, i, deltat)
                    S_out[k] <- out$S[simlength]
                    M_out[k] <- out$M[simlength]
                    A_out[k] <- out$M[simlength]
                    W_out[k] <- out$W[simlength]
                    incProgress(deltat/Tend, detail = paste("Computing trajectory ts = ", i))
                    k = k + 1
                }
                ts <- which.max(M_out[S_out>=0])
                if(input$goal_f == "stor"){
                    ts <- which.max(S_out)
                }
                else if (input$goal_f == "biost"){
                    valid <- S_out>=0;
                    ts <- which.max(S_out[valid]+M_out[valid])
                }
                out <- simple_water_model_sim_time_breaks(M0,S0,W0, Tend, kp, kr, ks, kw, ts, deltat)
                tcrit <- out$t[which.min(out$W)]
                
                
                df <- data.frame(t = out$t, M = out$M, S = out$S, W = out$W)
                ts_label <- paste0(c("time of allocation switch = ", as.character(ts)), collapse="")
                tcrit_label <- paste0(c("time water runs out = ", as.character(tcrit)), collapse="")
                
                ggplot(df, aes(x=t)) +
                    geom_line(aes(y=M, colour="Biomass")) +
                    geom_line(aes(y=S,colour="Storage")) +
                    theme_bw() +
                    geom_vline(xintercept=ts, colour="grey") +
                    geom_text(aes(x=ts, label=ts_label, y=50), colour="blue", angle=90, text=element_text(size=11), vjust = 1,) +
                    geom_vline(xintercept=tcrit, colour="grey") +
                    geom_text(aes(x=tcrit, label=tcrit_label, y=50), colour="blue", angle=90, text=element_text(size=11), vjust = 1,) +
                    xlab("Time (days)") +
                    ylab("Pool size [kgC]") +
                    scale_colour_manual("",
                                        values = c("Biomass"="#0072B2", "Storage"="#D55E00"))
                
                # Number of times we'll go through the loop
                
                
                # Increment the progress bar, and update the detail text.
                # incProgress(1/n, detail = paste("Doing part", i))
                
            })
        })
    }

# Run the application 
shinyApp(ui = ui, server = server)
