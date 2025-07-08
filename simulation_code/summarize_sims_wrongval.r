## Code to summarize simulations with validation study with wrong params ##

# SET UP ENVIRONMENT ----
library(tidyverse)
library(ggthemr)
library(cowplot)


# read in all data ----
dfraw <-
  list.files(path = "final_code/sims/main_paper/r1/output/wvraw200", pattern = "\\.csv", 
  full.names = T) %>% 
  map_df(~read_csv(.))
dfraw

dfraw2 <-
  list.files(path = "final_code/sims/main_paper/r1/output/wvraw2500", pattern = "\\.csv", 
  full.names = T) %>% 
  map_df(~read_csv(.))
dfraw2




# make plots ----

# inputs to get truth 
tau <- 2                                     # end of follow up
p_event <- 0.2                               # risk at end
alpha <- 2                                   # weibull shape
lambda <- (-log(1-p_event))^(1/alpha)/tau    # weibull scale for main study

# arrange data 
pl200 <- dfraw  %>% 
            filter(time == tau) %>% 
            mutate(truerisk =  1 - exp(-(lambda*time)^alpha), 
                   bias = risk - truerisk) %>% 
            mutate(method = case_when(method == "naive" ~ "0_naive", method == "val" ~ "1_val", 
                   method == "npest" ~ "2_npest", 
                   method == "pest" ~ "3_pest"), 
                   grouping = paste0(sc, method))

pl2500 <- dfraw2  %>% 
            filter(time == tau) %>% 
            mutate(truerisk =  1 - exp(-(lambda*time)^alpha), 
                   bias = risk - truerisk) %>% 
            mutate(method = case_when(method == "naive" ~ "0_naive", method == "val" ~ "1_val", 
                   method == "npest" ~ "2_npest", 
                   method == "pest" ~ "3_pest"), 
                   grouping = paste0(sc, method)) #%>% 


plotsims_server <- function(pldat, guide = F){
  #ggthemr("fresh")
  test <- ggplot() +
    geom_hline(aes(yintercept = 0), color = "black") + 
    geom_boxplot(aes(x = sc, y = bias * 100, fill = factor(method), group = factor(grouping)),
                 size = .3,  data = pldat, alpha = .8, 
                 outliers = F, position = position_dodge((width = 0.8))
                 ) +     
    scale_x_continuous(name = "Scenario", breaks = c(1, 2, 3, 4, 5, 6), 
                       labels = c("A", "B", "C", "D", "E", "F"))+
    scale_y_continuous(name = "Bias", 
                        limits = c(-10, 30), 
                        # breaks = c(-10, -5, 0, 5, 10, 15, 20, 25)) 
    ) +
    theme_classic() +
    theme(text = element_text(size = 12))

  if(guide == T){
    test <- test +
      scale_fill_npg(name = "", labels = c("Naive", "Validation only", "Nonparametric estimator", 
                                                "Parametric estimator")) +
      theme(legend.position = c(0.65, 0.88),
            legend.background = element_rect(fill = NA, color = NA),
            legend.key = element_rect(fill = NA), 
            legend.text = element_text(size = 12))
  }

  else{
    test <- test +  
    scale_fill_npg(guide = F) +
        theme(legend.position = "none")
  }
  return(test)
}

library(ggsci)
plot200 <- plotsims_server(pl200, guide = F)
plot2500 <- plotsims_server(pl2500, guide = T)

# ggsave(plot2500, file = "final_code/sims/main_paper/r1/output/box2500.pdf", width = 22, height = 15, units = "cm")
# ggsave(plot200, file = "final_code/sims/main_paper/r1/output/box200.pdf", width = 22, height = 15, units = "cm")
# 
# bothpl <- plot_grid(plot200, plot2500, ncol = 2, labels = c("A", "B"))
# ggsave(bothpl, file = "final_code/sims/main_paper/r1/output/bothpl.pdf", width = 38, height = 15, units = "cm")


# plot rmst -----

truermst <- integrate(function(x) exp(-(lambda*x)^alpha), lower = 0, upper = tau)[[1]]

plotsims_rmst <- function(pldat, guide = F){
  #ggthemr("fresh")
  test <- ggplot() +
    geom_hline(aes(yintercept = 0), color = "black") + 
    geom_boxplot(aes(x = sc, y = (rmst - truermst) * 365.25, fill = factor(method), group = factor(grouping)),
                 size = .3,  data = pldat, alpha = .8, 
                 outliers = F, position = position_dodge((width = 0.8))
    ) +     
    scale_x_continuous(name = "Scenario", breaks = c(1, 2, 3, 4, 5, 6), 
                       labels = c("A", "B", "C", "D", "E", "F"))+
    scale_y_continuous(name = "Bias in RMST (in days)" ,
                       limits = c(-80, 30), 
                       # breaks = c(-10, -5, 0, 5, 10, 15, 20, 25)) 
    ) +
    theme_classic() +
    theme(text = element_text(size = 12))
  
  if(guide == T){
    test <- test +
      scale_fill_npg(name = "", labels = c("Naive", "Validation only", "Nonparametric estimator", 
                                           "Parametric estimator")) +
      theme(legend.position = c(0.83, 0.88),
            legend.background = element_rect(fill = NA, color = NA),
            legend.key = element_rect(fill = NA), 
            legend.text = element_text(size = 12))
  }
  
  else{
    test <- test +  
      scale_fill_npg(guide = F) +
      theme(legend.position = "none")
  }
  return(test)
}

library(ggsci)
plot_r_200 <- plotsims_rmst(pl200, guide = F)
plot_r_2500 <- plotsims_rmst(pl2500, guide = F)

plotmat <- plot_grid(plot200, plot2500, plot_r_200, plot_r_2500, nrow = 2, labels = "AUTO")
ggsave(plotmat, file = "final_code/sims/main_paper/r1/output/wrongvalpl.pdf", width = 25, height = 20, units = "cm")
