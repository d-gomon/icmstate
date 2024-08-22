test_that("visualise_data", {
  skip_if(TRUE)
  #Just visualising data
  
  mockvisdata <- data.frame(delta = c(0,0,1,1), 
                            Delta = c(0,1,0,1), 
                            L = c(NA, NA, 2, 2.5), 
                            R = c(NA, NA, 3, 4),
                            time = c(4,5,6,7))
  rownames(mockvisdata) <- c("Group1", "Group2", "Group3", "Group4")
  
  visualise_data(mockvisdata)
  
  faildata <- data.frame(delta = c(0,0,1,1), 
                         Delta = c(0,1,0,1), 
                         L = c(NA, NA, 2, 2.5), 
                         R = c(NA, NA, 3, 4),
                         time = c(4,5,6,7))
  
  #Largerdata
  set.seed(4)
  Lt = rnorm(20, 10, 2)
  Rt <- Lt + rnorm(20, 6, 2)
  timet <- Rt + rnorm(20, 4, 0.5)
  largerdata <- data.frame(delta = rbinom(20, 1, prob = 0.5),
                           Delta = rbinom(20, 1, prob = 0.5),
                           L = Lt,
                           R = Rt,
                           time = timet)
  visualise_data(largerdata)
  
  
  
  #Visualising data + visualising support
  
  
  msmlarger <- msm_frydman(largerdata)
  visualise_data(largerdata, msmlarger)
  
  #When e_max >= s_max > R_max, F_{12} is not defined on [s_max, Inf]
  largerdata2 <- largerdata
  largerdata2$time[1] <- 26
  
  msmlarger2 <- msm_frydman(largerdata2)
  visualise_data(largerdata2, msmlarger2)
  #correct!
  
})















