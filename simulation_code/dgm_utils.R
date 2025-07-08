## DGM functions ##

#' a function to generate survival data
gen <- function(i, n, lambda, alpha = 2, tau = 2, cenrate = 0.4){
  set.seed(i) 
  # confounder x and exposure a
  id <- c(1:n)
  
  # time to event
  h <- exp(log(lambda)) 
  ts<- rweibull(n, shape=alpha, scale=1/h)
  
  # time to censoring
  u <- runif(n, 0, 1)
  c <- runif(n=n, min=0, max=tau-0.001)*(u<cenrate) + 
    tau * (u>=cenrate)
  
  # observed event times and types
  t <- pmin(ts, c, tau)
  j <- ifelse(t == c, 0, 1)
  
  dat <- data.frame(id,  t, ts, j, c)
  return(dat)
}



misclassify <- function(data, se, fp_rate, tp_rate = 0){
  # potentially misclassified event types
  data <- data %>% 
    mutate( # apply se to get jstar
      idevent = ifelse(j == 1, rbinom(nrow(data), 1, se), 0),
      jstar = ifelse(j == 1 & idevent == 1, 1, 0), 
      # compute tstar based on t and jstar
      tstar1 = ifelse(j == 1 & jstar == 1, t , 
                      ifelse(j == 1 & jstar == 0, pmin(tau, c), t)),
      
      # add in late detections
      delay = rexp(nrow(data), rate = tp_rate),
      tstar2 = ifelse(idevent == 0 & j == 1, t + delay, tstar1),
      tstar2 = pmin(tau, c, tstar2, na.rm = T),
      jstar = ifelse(tstar2 < pmin(tau, c), 1, jstar),
      
      # add in false positives
      t_to_fp = rexp(nrow(data), 
                     rate = fp_rate),
      t_to_fp  = ifelse(is.na(t_to_fp), tau, t_to_fp),
      # update tstar based on false positives
      tstar3 = pmin(t_to_fp, c, tau), 
      
      # compute tstar and jstar
      tstar = pmin(tstar3, tstar2),
      # update jstar based on new tstar
      jstar = ifelse(tstar < t, 1, jstar))
}

#' #' a function to induce outcome measurement error
#' #' @returns a data frame with new variables jstar (mismeasured event indicator) and tstar (mismeasured event time)
#' misclassify <- function(data, se, fp_rate, tp_rate = 0){
#'   # potentially misclassified event types
#'   data <- data %>% 
#'     mutate( # apply se to get jstar
#'       idevent = ifelse(j == 1, rbinom(nrow(data), 1, se), 0),
#'       jstar = ifelse(j == 1 & idevent == 1, 1, 0), 
#'       # compute tstar based on t and jstar
#'       tstar1 = ifelse(j == 1 & jstar == 1, t , 
#'                       ifelse(j == 1 & jstar == 0, pmin(tau, c), t)),
#'       
#'       # add in late detections
#'       delay = rexp(nrow(data), rate = tp_rate),
#'       tstar2 = ifelse(idevent == 0 & j == 1, t + delay, tstar1),
#'       tstar2 = pmin(tau, c, tstar2, na.rm = T),
#'       jstar = ifelse(tstar2 < pmin(tau, c), 1, jstar),
#'       
#'       # add in false positives
#'       t_to_fp = rexp(nrow(data), 
#'                      rate = fp_rate),
#'       t_to_fp  = ifelse(is.na(t_to_fp), tau, t_to_fp),
#'       # update tstar based on false positives
#'       tstar3 = pmin(t_to_fp, c, tau), 
#'       
#'       # compute tstar and jstar
#'       tstar = pmin(tstar3, tstar2),
#'       # update jstar based on new tstar
#'       jstar = ifelse(tstar < t, 1, jstar))
#' }
