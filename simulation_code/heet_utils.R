# heet utils #
## Utilities for ``heet" package ##


##' a function to compute the risk function
##' @param t a vector of times
##' @param delta a vector of event indicators
##' @param ci indicator of whether the function should output 95% CI on risk function (default is F)
##' @returns dataframe with risk function and standard error at each observed event time
##' @export
##' @examples
##' est.riskfxn(times, indicators)
##' 
##' 
##'
est.riskfxn <- function(t, delta, ci = F){
  require(survival)
  obsrisk <- survfit(Surv(t, delta) ~ 1, data = data.frame(t = t, delta = delta))
  if(ci == F){
    r_obs <- data.frame(time = summary(obsrisk)$time,
                        risk = 1 - summary(obsrisk)$surv,
                        se = summary(obsrisk)$std.err)
  }
  
  if(ci == T){
    r_obs <- data.frame(time = summary(obsrisk)$time,
                        risk = 1 - summary(obsrisk)$surv,
                        se = summary(obsrisk)$std.err, 
                        lcl = 1 - summary(obsrisk)$upper, 
                        ucl = 1 - summary(obsrisk)$lower)
  }
  return(r_obs)
}



##' a function to output risk at specific timepoints
##' @param time a vector of times
##' @param risk a vector of risk estimates
##' @param taus a vector of timepoints of interest.
##' @param output_se an indicator of whether to output standard error (T/F)
##' @export
##' @examples
##' get.risk(obsrisk$time, obsrisk$risk, taus = 4*365.25)
##' 
##'
get.risk <- function(time, risk, taus, se = NULL, output_se = F){
  require(dplyr)
  
  riskattaus <- vector()
  riskfxn <- data.frame(time = time, risk = risk)
  for(r in 1:length(taus)){
    riskattaus[r]  <- riskfxn %>%
      slice(which.min(abs(taus[r]-time))) %>%
      pull(risk)
  }
  
  if(!is.null(se) & output_se == T){
    seattaus <- vector()
    riskfxn <- data.frame(time = time, se = se)
    for(r in 1:length(taus)){
      seattaus[r]  <- riskfxn %>%
        slice(which.min(abs(taus[r]-time))) %>%
        pull(se)
    }
    return(data.frame(time = taus, risk = riskattaus, se = seattaus))
  }
  else return(riskattaus)
}


##' a function to compute a(t) and b(t) nonparametrically in validation data
##' @param times unique observed event times
##' @param obs_times observed event times in validation data
##' @param true_times true event times in validation data
##' @param obs_events observed event indicator in validation data
##' @param true_events true event indicator in validation data
##' @returns a dataframe with times and values of a and b at each time
##' @export
##' @examples
##' est.np.ab(obsrisk$time, val$w, val$t, val$eta, val$delta)
##'
est.np.ab <- function(times, obs_times, true_times, obs_events, true_events){
  require(dplyr)
  jump_times <- times #times when the curves are allowed to jump (observed event times)
  if (length(obs_times) != length(true_times)) {
    stop("Error: observed time vector should have the same number of times as true time vector")
  }
  if (length(obs_events) != length(true_events)) {
    stop("Error: observed vector of event indicators should have the same number of times as true event indicator vector")
  }
  ## empirical a(t)
  wmat <- (outer(obs_times, jump_times, function(a, b) as.integer(b-a>0 | near(a, b)))) * obs_events
  tmat <- (outer(true_times, jump_times, function(a, b) as.integer(b-a>0 | near(a, b) ))) * true_events
  at <- colSums(wmat*tmat)/colSums(tmat)
  at <- ifelse(is.na(at), 1, at)
  
  ## empirical b(t)
  tbmat <- (outer(true_times, jump_times, function(a, b) as.integer(b-a<0 )))
  bt <- colSums(wmat*tbmat)/colSums(tbmat)
  bt <- ifelse(is.na(bt), 0, bt)
  
  return(data.frame(time = times, a = at, b = bt))
}


#' a function to estimate fp_rate, event detection rate, and theta from validation data
#' this version repeats for each observed event time to get curves correct
#' @param times vector of ordered unique observed event times (if computing risk function) or tau
#' @param obs_times obs event times in validation data
#' @param true_times true event times in validation data
#' @param obs_events obs event indicator in validation data
#' @param true_events true event indicator in validation data
#' @returns a matrix with the estimated false positive rate, detection rate, and theta at each time
#' @examples
#' est.mc.params(obsrisk$time, val$w, val$t, val$eta, val$delta)
#' 
##' @export
##'
est.mc.params <- function(times, obs_times, true_times, obs_events, true_events, suppress = F){
  
  if(suppress == F) message("This function recomputes lambda_fp, lambda_d, and theta at each time supplied in the `times` vector")
  
  ests <- matrix(nrow = length(times), ncol = 3)
  for(i in 1:length(times)){
    # redefine times, censoring at each event time
    t <-     ifelse(true_times > times[i], times[i], true_times)
    w <- ifelse(obs_times > times[i], times[i], obs_times)
    delta <-     ifelse(true_times <= times[i], true_events, 0)
    eta <- ifelse(obs_times <= times[i], obs_events, 0)
    
    # compute parametric ests
    fp <- ifelse((delta == 0  & eta == 1) | (w < t &  delta == 1 & eta == 1), 1, 0)
    pt <- pmin(t, w)
    est_fp_rate <- sum(fp)/sum(pt)
    
    tp <- ifelse(delta == 1 & eta == 1 & w > t, 1, 0)
    pt <- sum(pmax(t, w) - t)
    est_d_rate <-  ifelse(pt == 0, 0, sum(tp)/sum(pt))
    
    est_theta <- sum(delta * eta  * as.numeric(w == t))/sum(delta * as.numeric(w >= t))
    est_theta <- ifelse(is.na(est_theta), 1, est_theta)
    ests[i,] <- c(as.numeric(est_fp_rate), as.numeric(est_d_rate), as.numeric(est_theta))
    colnames(ests) <- c("lambda_fp", "lambda_d", "theta")
  }
  return(ests)
}



##' a function to implement nonparametric correction
##' @param times a vector of ordered unique observed event times
##' @param obsfxn a vector of risk estimates of each time
##' @param afxn a vector of values of a(t)
##' @param bfxn a vector of values of b(t)
##' @returns dataframe with times and corrected risks
##' @export
##' @examples
##' np.cor(obsrisk$time, obsrisk$risk, a, b)
##' 
np.cor <- function(times, obsfxn, afxn, bfxn){
  ll <- list(times, obsfxn, afxn, bfxn)
  lengths <- sapply(ll, length)
  if(length(unique(lengths)) == 1) {
    a_t <-  ifelse(is.na(afxn), 1, afxn)
    b_t <-  ifelse(is.na(bfxn), 0, bfxn)
    tmp <-  (obsfxn - b_t)/(a_t - b_t)
    tmp <-  ifelse(is.na(tmp) | tmp < 0 | !is.finite(tmp) | tmp > 1, obsfxn, tmp)
    risk <-  cummax(tmp)
    return(data.frame(time = times, risk = risk))
  }
  else {
    stop("Error: lengths of time vector, observed risk function, a(t) and b(t) must be the same")
  }
}


##' a function to implement parametric correction
##' @param times a vector of ordered unique observed event times
##' @param obsfxn a vector of risk estimates of each time
##' @param fp_rate estimated false positive rate (can be a vector of the same length as time or a scalar)
##' @param d_rate estimated detection rate (can be a vector of the same length as time or a scalar)
##' @param theta estiamted probability of observing an event on time (can be a vector of the same length as time or a scalar)
##' @examples
##' p.cor(obsrisk$time, obsrisk$risk, est_Fp_rate, est_d_rate, est_theta)
##' 
##' @export
p.cor <- function(times, obsfxn, fp_rate, d_rate, theta){
  if(length(times) != length(obsfxn)){
    stop("Error: lengths of time vector and observed risk function must be equal")
  }
  if((length(fp_rate) != length(times) ) & length(fp_rate) != 1){
    stop("Error: fp_rate must be a vector of the same length as times or a scalar")
  }
  if((length(d_rate) != length(times) ) & length(d_rate) != 1){
    stop("Error: d_rate must be a vector of the same length as times or a scalar")
  }
  if((length(theta) != length(times) ) & length(theta) != 1){
    stop("Error: theta must be a vector of the same length as times or a scalar")
  }
  b <- (1 - exp(-fp_rate * times))
  y <-  1 - exp(-fp_rate * (times/2))
  x <-  theta * (1 - y)
  z <- (1 - exp(-d_rate * (times/2))) * (1 - x) * (1 - y)
  a <- y + x + z
  tmp  <-  (obsfxn - b)/(a - b)
  tmp = ifelse(is.na(tmp) | tmp < 0 | !is.finite(tmp) | tmp > 1, obsfxn, tmp)
  risk = cummax(tmp)
  return(data.frame(time = times, risk = cummax(risk)))
}


##' a function to perform nonparametric correction
##' @param obstimes observed event times in main study
##' @param obseta observed event indicators in main study
##' @param valw error-prone event times in validation study
##' @param valt gold standard event times in validation study
##' @param valeta error-prone event indicators in validation study
##' @param valdelta gold standard event indicators in validation study
##' @param taus_ time points of interest
##' @returns estimated risk at each time in `taus_`
##' @examples
##' analysis.np(dat$w, dat$eta, val$w, val$t, val$eta, val$delta, 4*365.25)
##' 
##' @export
analysis.np <- function(obstimes, obseta, valw, valt, valeta, valdelta, taus_){
  if(length(obstimes) != length(obseta)){
    stop("Error: lengths of obstimes vector and obseta must be equal")
  }
  if(length(valw) != length(valt) | length(valw) != length(valeta) |
     length(valeta) != length(valdelta)){
    stop("Error: lengths of times and indicators in the validation data must be equal")
  }
  obsriskfxn <- est.riskfxn(obstimes, obseta)
  emp_ab <- est.np.ab(obsriskfxn$time, valw, valt, valeta, valdelta)
  est_np <- np.cor(obsriskfxn$time, obsriskfxn$risk, emp_ab$a, emp_ab$b) %>%
    select(time, risk) %>%
    mutate(method = "np")
  est_npt <- get.risk(est_np$time, est_np$risk, taus = taus_)
  est_npt
}


##' a function to perform parametric correction
##' @description
##' `analysis.p()` computes risk functions corrected for misclassification using the parametric estimator. 
##' This function constrains the misclassificaiton parameters lambda_fp, lambda_d, and theta to be constant over time. 
##' 
##' @param obstimes observed event times in main study
##' @param obseta observed event indicators in main study
##' @param valw error-prone event times in validation study
##' @param valt gold standard event times in validation study
##' @param valeta error-prone event indicators in validation study
##' @param valdelta gold standard event indicators in validation study
##' @param taus_ time points of interest
##' @returns estimated risk at each time in `taus_`
##' @examples
##' analysis.p(dat$w, dat$eta, val$w, val$t, val$eta, val$delta, 4*365.25)
##' 
##' @export
analysis.p <- function(obstimes, obseta, valw, valt, valeta, valdelta, taus_){
  if(length(obstimes) != length(obseta)){
    stop("Error: lengths of obstimes vector and obseta must be equal")
  }
  if(length(valw) != length(valt) | length(valw) != length(valeta) |
     length(valeta) != length(valdelta)){
    stop("Error: lengths of times and indicators in the validation data must be equal")
  }
  obsriskfxn <- est.riskfxn(obstimes, obseta)
  params <- as.data.frame(est.mc.params(max(taus_), valw, valt, valeta, valdelta, suppress = T))
  est_fp_rate <- ifelse(is.na(params[,1]), 0, params[,1])
  est_d_rate <- ifelse(is.na(params[,2]), 0, params[,2])
  est_theta <- params[,3]
  est_p <- p.cor(obsriskfxn$time, obsriskfxn$risk, est_fp_rate, est_d_rate, est_theta) %>%
    select(time, risk) %>%
    mutate(method = "p")
  est_pt <- get.risk(est_p$time, est_p$risk, taus = taus_)
  return(est_pt)
}


##' a function to bootstrap nonparametric and parametric estimators
##' @param B number of bootstrap iterations
##' @param bootdata main study data
##' @param valdata validation study data
##' @param taus_ timepoints of interest
##' @param risk_func choice of risk function. specify analysis_np, analysis_p, or a user-defined risk function with the same arguments
##' @returns a matrix with B rows containing risk at each time taus_ (in columns)
##' @export
##'
misclass.boot <- function(B, bootdata, valdata, taus_, risk_func){
  require(resample)
  if(B == 0) return(0)
  if(B>0){
    bsr <- matrix(NaN, nrow = B, ncol = length(taus_))
    datbi<-samp.bootstrap(nrow(bootdata), B)
    valbi<-samp.bootstrap(nrow(valdata), B)
    for(k in 1:B){
      dati <- bootdata[datbi[,k],]
      vali <- valdata[valbi[,k],]
      bsr[k,] <- risk_func(dati$wstar, dati$eta,
                           vali$wstar, vali$tstar,
                           vali$eta,
                           vali$delta, taus_)
    }
    return(bsr)
  }
}
