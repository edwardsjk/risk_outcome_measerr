## Apply package to paper analysis ##

pak::pak("edwardsjk/heet")

library(heet)
library(tidyverse)
library(ggsci)


tau <- 10

# read in data ----

datadir <- "/Volumes/Jess Edwards/UCHCC_data/deduplicated_data_1120"
dat <- read_csv(paste0(datadir, "/analysisdat.csv"))

## create times ----
dat_analyze <- dat %>%  
  filter(uncentrytocare < as.Date("2022-12-31")) %>% # limit to those enrolled in 2022 or earlier
  mutate(origin = uncentrytocare, 
         cen_dt = pmin(as.Date("2022-12-31")), 
         c = pmin(tau, as.numeric(cen_dt - origin)/365.25),
         w = as.numeric(dod_link - origin)/365.25, # tstar in years
         t = as.numeric(DOD - origin)/365.25, 
         w = replace_na(w, tau + 1), 
         t = replace_na(t, tau + 1),
         tstar = pmin(t, c), 
         wstar = pmin(w, c), 
         delta = ifelse(tstar == t, 1, 0), 
         eta = ifelse(wstar == w, 1, 0)) %>% 
  select(patient_key, origin, tstar, wstar, delta, eta, race, ethnicity, gender, riskfactor, birth_date) %>% 
  filter(origin > as.Date("2001-01-01") & wstar>=0) 

# truth ----

trueriskfxn <- est.riskfxn(dat_analyze$tstar, dat_analyze$delta, ci = T)
true_tau <- get.risk(trueriskfxn$time, trueriskfxn$risk, tau, se = trueriskfxn$se, output_se = T)



# naive ----
obsriskfxn <- est.riskfxn(dat_analyze$wstar, dat_analyze$eta, ci = T)
obs_tau <- get.risk(obsriskfxn$time, obsriskfxn$risk, tau, se = obsriskfxn$se, output_se = T)

# sample validation data ----
valdat <- read_csv(file="/Volumes/Jess Edwards/analysis_datasets/val.csv")

# validation study only ----
valriskfxn <- est.riskfxn(valdat$tstar, valdat$delta, ci = T)
val_tau <- get.risk(valriskfxn$time, valriskfxn$risk, tau, se = valriskfxn$se, output_se = T)

# nonparametric correction ----

###  empirical (np) a and b ----
emp_ab <- est.np.ab(obsriskfxn$time, valdat$wstar, valdat$tstar, valdat$eta, valdat$delta)
at <- ggplot() +
  geom_step(aes(x = time, y = a), data = emp_ab) +
  ylab("a(t)") +
  xlab("Time") +
  theme_bw()
bt <- ggplot() +
  geom_step(aes(x = time, y = b), data = emp_ab) +
  ylab("b(t)") +
  xlab("Time") +
  ylim(-0.05, 0.05)+
  theme_bw()
atbt <- cowplot::plot_grid(at, bt, labels = c("A", "B"))
ggsave(atbt, filename = "output/atbt.png", width = 20, height = 10, units = "cm")

### do empirical correction ----
est_np <- np.cor(obsriskfxn$time, obsriskfxn$risk, emp_ab$a, emp_ab$b) %>% 
  select(time, risk) %>% 
  mutate(method = "np")

point_est_np <- analysis.np(dat_analyze$wstar, dat_analyze$eta, valdat$wstar, valdat$tstar, valdat$eta, valdat$delta, 10)

# bootstrap to get se at end of study
boot_ests_np <- misclass.boot(500, dat_analyze, valdat, tau, analysis.np)
np_se <- sd(boot_ests_np)

# bootstrap for pointwise confidence intervals
boot_ests_np_all <- misclass.boot(500, dat_analyze, valdat, obsriskfxn$time, analysis.np)
np_se_all <- apply(boot_ests_np_all, 2, FUN = sd)

est_np$lcl <- est_np$risk - 1.96*np_se_all
est_np$ucl <- est_np$risk + 1.96*np_se_all


### parametric correction ----
params <- as.data.frame(est.mc.params(10, valdat$wstar, valdat$tstar, valdat$eta, valdat$delta))
est_p <- p.cor(obsriskfxn$time, obsriskfxn$risk, tail(params$lambda_fp, n=1), tail(params$lambda_d, n = 1), tail(params$theta, n = 1)) %>% 
  select(time, risk) %>% 
  mutate(method = "p")

point_est_p <- analysis.p(dat_analyze$wstar, dat_analyze$eta, valdat$wstar, valdat$tstar, valdat$eta, valdat$delta, 10)

#### bootstrap ----
# at end of follow up
boot_ests_p <- misclass.boot(500, dat_analyze, valdat, tau, analysis.p)
p_se <- sd(boot_ests_p)
p_se
lcl <- point_est_p - 1.96*p_se
ucl <- point_est_p + 1.96 * p_se
print(c(lcl, ucl))

# for entire curve
boot_ests_p_all <- misclass.boot(500, dat_analyze, valdat, obsriskfxn$time, analysis.p)
p_se_all <- apply(boot_ests_p_all, 2, FUN = sd)
est_p$lcl <- est_p$risk - 1.96*p_se_all
est_p$ucl <- est_p$risk + 1.96*p_se_all

### pseudo rmse
prmse <- sqrt((true_tau$risk - point_est_p)^2 + p_se^2)
prmse

# plot ----

est_obs <- obsriskfxn %>% mutate(method = "obs") %>% select(-se)
est_val <- valriskfxn %>% mutate(method = "v") %>% select(-se)
est_tru <- trueriskfxn %>% mutate(method = "t") %>% select(-se)

all <- bind_rows(est_tru, est_obs, est_val, est_p, est_np) %>% 
  mutate(method2 = factor(method, levels = c("obs", "v", "np", "p", "t"))) 

flat <- all %>% group_by(method) %>% slice(n()) %>%  mutate(time = tau)
all <- bind_rows(all, flat)

respl <- ggplot() + 
  geom_step(aes(x = time, y = risk, color = method2), data = all) +
  theme_classic() +
  scale_x_continuous(name = "Time (years) from entry into care", 
                     limits = c(0, tau), expand = c(0, 0), 
                     breaks = seq(0, 10, by = 2)) + 
  scale_y_continuous(name = "Risk", breaks = c(0, .05, .1, .15, .2, .25, .3, .35, .4), 
                     limits = c(0, 0.15), expand = c(0,0), minor_breaks = seq(0,.15, by = .025)) +
  scale_color_npg(breaks = c("t", "obs", "v", "np", "p"), labels = c("Gold Standard", "Naive", "Validation Data", "Proposed, Nonparametric", "Proposed, Parametric" )) + #"Proposed, Parametric, iterative",
  theme(legend.position = c(0.3, 0.8), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_blank(), 
        text = element_text(size = 14))
respl
#ggsave(respl, file = "output/death_link_4.png", width = 12, height = 9, units = "cm")
ggsave(respl, file = "output/death_link_5.pdf", width = 12.5, height = 9, units = "cm")

## with confidence intervals----

labs <- c("t" = "Gold Standard", "obs" = "Naive", "v" = "Validation Data", "np" = "Nonparametric", "p" = "Parametric" )
resplcis <- ggplot() + 
  geom_step(aes(x = time, y = risk, color = method2), data = all) +
  geom_step(aes(x = time, y = lcl, color = method2), data = all, alpha = 0.5, lty = 2) +
  geom_step(aes(x = time, y = ucl, color = method2), data = all, alpha = 0.5, lty = 2) +
  facet_grid( ~ method2, labeller = labeller(method2 = labs)) + 
  theme_classic() +
  scale_x_continuous(name = "Time (years) from entry into care", 
                     limits = c(0, tau), expand = c(0, 0), 
                     breaks = seq(0, 10, by = 2)) + 
  scale_y_continuous(name = "Risk", breaks = c(0, .05, .1, .15, .2, .25, .3, .35, .4), 
                     limits = c(0, 0.2), expand = c(0,0), minor_breaks = seq(0,.15, by = .025)) +
  scale_color_npg(breaks = c("t", "obs", "v", "np", "p"), labels = c("Gold Standard", "Naive", "Validation Data", "Proposed, Nonparametric", "Proposed, Parametric" ), 
                  guide = "none") + #"Proposed, Parametric, iterative",
  theme(legend.position = c(0.3, 0.8), 
        legend.key = element_blank(), 
        legend.title = element_blank(), 
        legend.background = element_blank(), 
        text = element_text(size = 14), 
        panel.spacing.x = unit(0.5, "cm"), 
        plot.margin = margin(.5,.5,.5,.5, "cm"))
resplcis
#ggsave(respl, file = "output/death_link_4.png", width = 12, height = 9, units = "cm")
ggsave(resplcis, file = "output/death_link_cis.pdf", width = 20, height = 9, units = "cm")

