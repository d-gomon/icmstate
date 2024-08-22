test_that("msm_frydman", {
  skip_if(TRUE)
  testdat <- data.frame( Delta = c(0,1,0,1), delta = c(0,0,1,1), 
                         L = c(NA, NA, 2, 2.5), R = c(NA, NA, 3, 4),
                         time = c(4,5,6,7))
  msmtestdat <- msm_frydman(testdat)
  
  testdat2 <- data.frame( Delta = c(0,1,0,1), delta = c(0,0,1,1), 
                          L = c(7, 100, 2, 2.5), R = c(20, 1, 3, 4),
                          time = c(4,5,6,7), 
                          trunc = c(1,1,2,2))
  
  msm_frydman(testdat2)
  
  
  set.seed(1)
  Lt = rnorm(20, 10, 2)
  Rt <- Lt + rnorm(20, 6, 2)
  timet <- Rt + rnorm(20, 4, 0.5)
  largerdata <- data.frame(delta = rbinom(20, 1, prob = 0.5),
                           Delta = rbinom(20, 1, prob = 0.5),
                           L = Lt,
                           R = Rt,
                           time = timet)
  asd <- msm_frydman(largerdata)
  
  testdat3 <- data.frame(Delta = c(0,0,1), delta = c(0,1,1), 
                         L = c(NA,  2, 2.5), R = c(NA,  3, 4),
                         time = c(4,6,7))
  
  asd3 <- msm_frydman(testdat3)
  
  
  testdathudgens <- data.frame(delta = c(rep(1, 4), 0, 0),
                               Delta = c(rep(1,2), rep(0, 2), 0, 0),
                               L = c(1, 2, 4, 6, 20, 22),
                               R = c(3,5,7,8, 30, 33),
                               time = c(19, 20, 21, 22, 1, 1.2))
  
  msmhudgens <- msm_frydman(testdathudgens, tol = 1e-20)
  visualise_data(testdathudgens, msmhudgens)
  
  
  
  
  ###########################################################################
  ######################Truncated Data Testing###############################
  ###########################################################################
  
  testdat_trunc <- data.frame( Delta = c(0,1,0,1), delta = c(0,0,1,1), 
                               L = c(NA, NA, 2, 2.5), R = c(NA, NA, 3, 4),
                               time = c(4,5,6,7), trunc = c(2.2, 1, 0.5, 0))
  
  asd_trunc <- msm_frydman(testdat_trunc)
})









