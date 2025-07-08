## Code to run sims on server ##
### Repeat sims with different shape function in validation dataset (nval = 200) ###

# SCENARIO CONTROL ----
scenario <- as.integer(commandArgs(trailingOnly=T))
print(scenario)
theta_s <- c(1, .7, 0.7, 1, .7, .7)
tp_s <- c(0, 0, .3, 0, 0, .3)
fp_s <- c(0, 0, 0, .1, .1, .1)

# SET UP ENVIRONMENT ----
library(tidyverse)
library(survival)
library(doParallel)
library(foreach)

# source("code/utils.R")
source("code/sim_utilsv5_r1.R")
source("code/heet_utils.R")
source("code/dgm_utils.R")

today <- Sys.Date()

# SET INPUTS ----

nval <- 200
n <- 5000                                    # sample size
tau <- 2                                     # end of follow up
p_event <- 0.2                               # risk at end
alpha <- 2                                   # weibull shape
lambda <- (-log(1-p_event))^(1/alpha)/tau    # weibull scale for main study

alpha_val <- 0.5
p_event_val <- 0.35                               # risk at end
lambda_val <- (-log(1-p_event_val))^(1/alpha_val)/tau 
taus <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)  # timepoints of interest

J <- 10000                                    # num sims
B <- 0                                     # num bootstraps

# initialize stream of random numbers
randos <- sample(1:(J*10000), J, replace = F)


# SET UP PROCESSING ----

cl<-makeCluster(64) # max should be 125, but cannot get to run with > 64?
registerDoParallel(cl)

start <- Sys.time()


# TRUTH ----

times <- seq(from = 0, to = 2, by = 0.001)
r0 <- data.frame(time = times, 
                 risk = 1 - exp(-(lambda*times)^alpha))

# RUN SCENARIO ----
theta <- theta_s[scenario]
fp_prob <- fp_s[scenario]
tp_prob <- tp_s[scenario]
fp_rate <- ifelse(fp_prob>0, -log((1 - (fp_prob)))/tau, 0)
tp_rate <- ifelse(tp_prob>0, -log((1 - (tp_prob)))/tau, 0)

res <- foreach(s=1:J, .combine='rbind',
             .packages = c("survival", "tidyr", "dplyr", "doParallel", "resample")) %dopar%{
               wrap(s, randos[s], taus, B, rmst = T)
             }


results <- res %>% 
    pivot_longer(cols = c(naive, val, npest, pest), values_to = "risk", names_to = "method") %>% 
    mutate(se = ifelse(method == "naive", naive_se, 
                       ifelse(method == "val", val_se, 
                              ifelse(method == "npest", se_np, se_p))),
           rmst = case_when(method == "naive" ~ rmst_naive, 
                            method == "val" ~ rmst_val, 
                            method == "npest" ~ rmst_np, 
                            method == "pest" ~ rmst_p),
           sc = scenario) %>% 
    select(k, time, method, risk, se, rmst, sc)
write.csv(results, file = paste0("output/altvalserver200", scenario, "_", today, ".csv"))

stopCluster(cl)
end <- Sys.time()
end-start

# SUMMARIZE ----

sumfun <- function(dat){
  sum_res <- dat %>% 
    mutate(grouping = paste0(method, time), 
           sc = scenario)
  sc_tab <- sumsims(lambda, alpha, sum_res, tau_=2)
  #sc <- plotsims(r0, sum_res, guide = labs)
  return(sc_tab)
}


sc_res <- sumfun(results)

# SAVE OUTPUT ----

write.csv(sc_res, file = paste0("output/altvaltab200_", today, "_", scenario, ".csv"))
