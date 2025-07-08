# wrapper function for running sims #


##' wrapper function

wrap <- function(i, seed, taus, B=0, rmst = F, perfectval = T){
  
  if(i%%10 == 0) print(i)
  
  # generate main study data 
  dat <- gen(seed, n, lambda, alpha, tau = 2, cenrate = 0.5)  # generate data
  obsdat <- dat %>%                                # misclassify outcomes
    misclassify(., theta, fp_rate, tp_rate) %>% 
    select(t, j, tstar, jstar, id, t_to_fp)
  
  # generate validation data 
  valdattmp <- gen(seed, nval, lambda_val, alpha_val, tau = 2, cenrate = 0.5)  # generate data
  valdat <- valdattmp %>%                                # misclassify outcomes
    misclassify(., theta, fp_rate, tp_rate) %>% 
    select(t, j, tstar, jstar, id, t_to_fp)
  
  if(perfectval == F){
    valdat2 <- valdattmp %>%                               
      misclassify(., 0.95, 0.005, 0.2) %>% 
      select(tstar, jstar, id) %>% 
      rename(t2 = tstar, j2 = jstar)
    
    valdat <- valdat %>% 
      left_join(valdat2) %>% 
      select(t2, j2, tstar, jstar, id, t_to_fp) %>% 
      rename(t = t2, j = j2)
  }
  
  # naive estimator
  naiversk <- est.riskfxn(obsdat$tstar, obsdat$jstar)                   
  naive <- get.risk(naiversk$time, naiversk$risk, taus, se = naiversk$se, output_se = T)
  #naive_se <- get_risk_se(naiversk, taus)
  
  
  
  # val estimator
  valrsk <- est.riskfxn(valdat$t, valdat$j)             
  val <- get.risk(valrsk$time, valrsk$risk, taus, se = valrsk$se, output_se = T)

  # get empirical (np) a and b
  emp_ab <- est.np.ab(naiversk$time, valdat$tstar, valdat$t, valdat$jstar, valdat$j)
  
  # do empirical correction
  est_np <- np.cor(naiversk$time, naiversk$risk, emp_ab$a, emp_ab$b) %>% 
    select(time, risk) 
  est_npt <- get.risk(est_np$time, est_np$risk, taus)
  
  # get parameters for parametric cor
  params <- (est.mc.params(max(taus), valdat$tstar, valdat$t, valdat$jstar, valdat$j, suppress = T))
  est_fp_rate <- ifelse(is.na(params[1]), 0, params[1])
  est_d_rate <- ifelse(is.na(params[2]), 0, params[2])
  est_theta <- params[3]
  
  # do parametric correction
  est_p <- p.cor(naiversk$time, naiversk$risk, est_fp_rate, est_d_rate, est_theta) %>% 
    select(time, risk)
  est_pt <- get.risk(est_p$time, est_p$risk, taus)
  
  # bootstrap stderrs
  
  cor_risk_boots <- bootloop(B, obsdat, valdat, taus, get.cor.all)
  
  if(B>0){
    se_np <- apply(cor_risk_boots[,1:length(taus)], 2, FUN = sd)
    se_p <- apply(cor_risk_boots[,(length(taus) + 1) : ncol(cor_risk_boots)], 2, FUN = sd)
  }
  
  else{
    se_np <- 0
    se_p <- 0
  }
  
  # compile results
  if(rmst == F){
    res <- data.frame(time = taus, naive = naive$risk, val = val$risk, npest = est_npt,
                    pest = est_pt, k=rep(i, length(taus)), 
                    naive_se = naive$se, val_se = val$se, se_np = se_np, se_p = se_p)
  }
  if(rmst == T){
    n_r <- get.rmst(naiversk)
    v_r <- get.rmst(valrsk)
    np_r <- get.rmst(est_np)
    p_r <- get.rmst(est_p)
    res <- data.frame(time = taus, naive = naive$risk, val = val$risk, npest = est_npt,
                      pest = est_pt, k=rep(i, length(taus)), 
                      naive_se = naive$se, val_se = val$se, se_np = se_np, se_p = se_p, 
                      rmst_naive = n_r, rmst_val = v_r, rmst_np = np_r, rmst_p = p_r)
  }
  return(res)
}



plotsims <- function(r0, resultsdat, guide = F, metric = "bias"){

  if(metric == "bias") resultsdat$y <- (risk - truerisk)*100
  if(metric == "relbias") resultsdat$y <- ((risk - truerisk)*100)/(truerisk * 100)
  
  resultsdat <- resultsdat %>% left_join(r0 %>% rename(truerisk = risk))
  test <- ggplot() +
    geom_hline(aes(yintercept = 0), color = "black") + 
    geom_boxplot(aes(x = time, y = y, fill = factor(method), group = grouping),
                 size = .3,  data = resultsdat,
                 outlier.shape = NA, alpha = .8 #,position = position_dodge((width = 0.2))
                 ) + #
    #geom_step(aes(x = time, y = risk), data = r0) +
    
    scale_x_continuous(name = "Time")+
    scale_y_continuous(name = "Bias", limits = c(-10, 25), 
                       breaks = c(-10, -5, 0, 5, 10, 15, 20, 25))

  if(guide == T){
    test <- test +
      scale_fill_discrete(name = "", labels = c("Naive", "Nonparametric estimator", 
                                                "Parametric estimator", "Validation only")) +
      theme(legend.position = c(0.35, 0.85),
            legend.background = element_rect(fill = NA),
            legend.key = element_rect(fill = NA), 
            legend.text = element_text(size = 10))
  }

  else{
    test <- test +  theme(legend.position = "none")
  }
  return(test)
}

#' function to summarize simulation results
#'
sumsims <- function(lambda, alpha, results, tau_){
  truerisk <-  1 - exp(-(lambda*tau_)^alpha)
  truermst <- integrate(function(x) exp(-(lambda*x)^alpha), lower = 0, upper = tau_)[[1]]
  res_tau <- results %>%
    filter(time == tau_) %>%
    group_by(sc, method) %>%
    mutate(bias_i = (risk - truerisk)*100,
           lcl = risk - 1.96*se,
           ucl = risk + 1.96*se,
            cov = ifelse(truerisk<=ucl & truerisk>=lcl, 1, 0), 
           bias_rmst_i = rmst - truermst) %>%
    summarize(bias = mean(bias_i, na.rm = T), ese = sd(bias_i, na.rm = T), 
              ase = mean(se * 100), cov = mean(cov, na.rm = T), 
              biasrmst = mean(bias_rmst_i, na.rm = T), mcse = sqrt(ese^2/10000)) %>% #
    mutate(rmse = sqrt(bias^2 + ese^2)) %>%
    ungroup() %>%
    mutate(method = case_when(method == "naive" ~ "0_naive", method == "val" ~ "1_val", 
            method == "npest" ~ "2_npest", 
            method == "pest" ~ "3_pest")) %>% 
    arrange((method)) %>%
    mutate(bias = round(bias, 1), ese = round(ese, 1), rmse = round(rmse, 1),
           ase = round(ase, 1), ser = round(ase/ese, 1), cov = round(cov, 2), 
           biasrmst = round(biasrmst, 2)) %>% 
    select(sc, method, bias, biasrmst, ese, rmse, ase, ser, cov, mcse) #, ase, ser, cov
}


#' a function to apply all estimators
get.cor.all <- function(obsdat, valdat, taus_){
  # naive estimator
  naiversk <- est.riskfxn(obsdat$tstar, obsdat$jstar)                   

  # get empirical (np) a and b
  emp_ab <- est.np.ab(naiversk$time, valdat$tstar, valdat$t, valdat$jstar, valdat$j)
  
  # do empirical correction
  est_np <- np.cor(naiversk$time, naiversk$risk, emp_ab$a, emp_ab$b) %>% 
    select(time, risk) 
  est_npt <- get.risk(est_np$time, est_np$risk, taus_)
  
  # get parameters for parametric cor
  params <- (est.mc.params(max(taus_), valdat$tstar, valdat$t, valdat$jstar, valdat$j, suppress = T))
  est_fp_rate <- ifelse(is.na(params[1]), 0, params[1])
  est_d_rate <- ifelse(is.na(params[2]), 0, params[2])
  est_theta <- params[3]
  
  # do parametric correction
  est_p <- p.cor(naiversk$time, naiversk$risk, est_fp_rate, est_d_rate, est_theta) %>% 
    select(time, risk)
  est_pt <- get.risk(est_p$time, est_p$risk, taus_)
  allests <- c(est_npt, est_pt)
  return(allests)
}

##' a function to do bootstrap in a loop
##' Useful for bootstrapping within the simulation
bootloop <- function(B, bootdata, valdata, taus_, risk_func){
  require(resample)
  if(B == 0) return(0)
  if(B>0){
    bsr <- matrix(NaN, nrow = B, ncol = 2*length(taus_))
    datbi<-samp.bootstrap(nrow(bootdata), B)
    valbi<-samp.bootstrap(nrow(valdata), B)
    for(k in 1:B){
      dati <- bootdata[datbi[,k],]
      vali <- valdata[valbi[,k],]
      bsr[k,] <- risk_func(dati, vali, taus_)
    } 
    return(bsr)
  }
}



get.rmst <- function(riskdat){
  times <- c(0, riskdat$time, max(taus))
  survival <- c(1, 1 - riskdat$risk)
  sum(diff(times)*survival)
}

