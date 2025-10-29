library(geepack)
# working correlation used for modeling
corstr_working <- "ar1"
# store the simulation results
results <- vector(mode = "list",length = nsim)

for(i in 1:nsim)
{

  dat <- gen_data(N,
                  baseprobs,
                  covs_effects,
                  time_effects,
                  time_trt_effects,
                  visits,
                  corMatrix)
  dat_wide <- dat$wide_format
  dat_long <- dat$long_format
  
  # create indicator for week 4 and week 8
  dat_long <- dat_long |>
    mutate(w4 = 1*(visit=="week4"),w8 = 1*(visit=="week8"))
  # model time categorically
    mod_lwo <- tryCatch(expr = {lwo(GBS_DS ~ trt + 
                                      trt:w4 + trt:w8+
                                      age + pre_diarrhea ,
                                    data = dat_long,
                                    id.var = "id",
                                    visit.var = "visit",
                                    time.vars = c("w4","w8"),
                                    corstr = corstr_working,
                                  #   std.err = "san.se")
                                    std.err = "san.se.modified")
      },
                        error = function(e)
                          {
                          return(NA)
                        })
    if(sum(is.na(mod_lwo)))
    {
      # estimated log win odds at each visit
      est_log_wos <- rep(NA,length(visits))
      # estimated variance for the log win odds
      var_log_wos <-  rep(NA,length(visits))
      # coverage on the log win odds scale
      CI_95 <- cbind( rep(NA,length(visits)), rep(NA,length(visits)))
      coverage_CI_95 <-  rep(NA,length(visits))
    } else
    {
      trans_matrix <- matrix(c(1,0,0,0,0,
                               1,0,0,1,0,
                               1,0,0,0,1),nrow = 3,byrow = T)
      L <- matrix(c(1,0,0,0,0,
                    0,0,0,1,0,
                    0,0,0,0,1),nrow = 3,byrow = T)
      # generalized Wald test for global null hypothesis
      W_2 <- t(L%*%mod_lwo$coefficients)%*%solve(L%*%mod_lwo$var%*%t(L))%*%(L%*%mod_lwo$coefficients)
      p_val_global_null <- pchisq(q=W_2,df=3,lower.tail = FALSE)
      # }
      # estimated log win odds at each visit
      est_log_wos <- trans_matrix%*%c(mod_lwo$coefficients)
      # estimated variance for the log win odds
      var_log_wos <- diag(trans_matrix%*%mod_lwo$var%*%t(trans_matrix))
      # coverage on the log win odds scale
      CI_95 <- cbind(est_log_wos - qnorm(0.975)*sqrt(var_log_wos),est_log_wos + qnorm(0.975)*sqrt(var_log_wos))
      coverage_CI_95 <- log(estimand) >= CI_95[,1] & log(estimand) <= CI_95[,2]
    }
  results[[i]] <- list(est_log_wos,var_log_wos,CI_95,coverage_CI_95,p_val_global_null)
  if(i%%10==0)
   {print(i)
   save(results,file=paste0("results/results-modified-N",N,"-rho",rho,"-",scenario,".RData"))
 }
}
