###########################################################
# function to assess simulation results
###########################################################
assess_results <- function(N=NULL,rho=NULL,scenario=NULL,results=NULL,estimand=NULL,
                           type = "comp")
{
  if(is.null(results))
  {
    load(paste0("C:/Users/ylong/OneDrive - LUMC/Projects/Longitudinal_win_odds/codes/results/results-",type,
                "-N",N,"-rho",rho,"-",scenario,".RData"))
  }
  # assess coverage of 95% CI for log win odds
  coverage <- rowMeans(sapply(results, function(i) 
  {
    i[[4]]
  }),na.rm=T)
  # power/type I error
  power_each_visit <- rowMeans(sapply(results, function(i) 
  {
    CI <- i[[3]]
    # CI does not contain 0
    CI[,1] > 0 | CI[,2] < 0
  }),na.rm=T)
  # assess bias of log win odds estimator
  est_log_wos <- sapply(results, function(i) return(i[[1]]))
  bias_log_wos <- rowMeans(est_log_wos,na.rm = TRUE) - log(estimand)
  # convert to odds scale
  est_wos <- exp(est_log_wos)
  bias_wos <- rowMeans(est_wos,na.rm = TRUE) - estimand
  # on the probability scale it should be unbiased
  bias_MWS <- rowMeans(est_wos/(1+est_wos),na.rm = TRUE) - estimand/(1+estimand)
  
  # assess estimated variance
  # "true" variance by directly calculating the variance of MC estimates
  var_log_wos_MC <- apply(est_log_wos, 1, function(x) var(x,na.rm = TRUE))
  # average estimated variance
  var_log_wos_avg <- rowMeans(sapply(results, function(i) return(i[[2]])),na.rm = TRUE)
  # global power
  power_global <- mean(sapply(results, function(i) i[[5]]<0.05))
  res_assessed <- data.frame(cbind(coverage,power_each_visit,estimand,bias_log_wos,bias_wos,bias_MWS,
                                   var_log_wos_MC,var_log_wos_avg))
  rownames(res_assessed) <- NULL
  out <- list("global power" = power_global,
              "table of performance metrics" = res_assessed)
  return(out)
}
