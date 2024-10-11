#############################################################
################Compare with Turnbull########################
#############################################################

test_that("Turnbull vs MSM",
{
  ##############Data preparation#################
  set.seed(1)
  #Absorbing state 2
  qmatrix <- rbind(
    c(-0.1, 0.1),
    c(0, 0)
    )
  
  n <- 3
  gd <- NULL
  for (i in 1:n) {
    smsm <- msm::sim.msm(qmatrix, 14)
    # inspection times uniformly on 2 year intervals until 14 years
    itimes <- seq(0, 12, by=2) + runif(7, 0, 2)
    gdi <- data.frame(id=i, state=evalstep(time=smsm$times, stepf=smsm$states,
                                           newtime=c(0,itimes)), time=c(0,itimes))
    # throw away superfluous duplicate state 3 (absorbing)
    gdi <- gdi[!(gdi$state==3 & duplicated(gdi$state)),]
    gd <- rbind(gd, gdi)
  }
  
  #Which states can be reached from each state?
  tmat <- mstate::transMat(x = list( c(2), c() ))
  
  icdata <- NULL
  for(i in unique(gd$id)){
    gdi <- subset(gd, id == i)
    L_idx <- Position(isTRUE, gdi$state <= 1, right = TRUE)
    L <- gdi$time[L_idx] #Exit time from state 1
    if(L_idx < nrow(gdi)){
      R <- gdi$time[L_idx + 1]
    } else{
      R <- Inf
    }
    icdata <- rbind(icdata, c(L, R))
  }
  icdata <- as.data.frame(icdata)
  colnames(icdata) <- c("L", "R")
  
  ############Apply algorithms####################
  MSMres <- npmsm(gd = gd, tmat = tmat, tol = 1e-5)
  supportHein <- support_npmsm(MSMres, cutoff = 1e-5)
  #Visualise:
  #theinsupport <- visualise_msm(gd, tmat, MSMres, cutoff = 1e-5)
  
  Turnbull2 <- icenReg::ic_np(cbind(L, R) ~ 0, data = icdata)
  #Turnbull intervals
  #Turnbull2$T_bull_Intervals
  #Mass in those intervals
  #Turnbull2$p_hat
  
  t2res <- cbind(t(Turnbull2$T_bull_Intervals), matrix(Turnbull2$p_hat))
  colnames(t2res) <- c("L", "R", "phat")
  expect_true(all.equal(t2res[, 1:2], supportHein$`State 1 -> State 2`$support[, 1:2]))
  
})




test_that("Frydman (1995) vs MSM",{
  tmat <- mstate::transMat(x = list( c(2, 3), c(3), c() ))
  set.seed(1)
  #Absorbing states 4 and 6
  qmatrix <- rbind(
    c(-0.15, 0.1, 0.05),
    c(0, -0.1, 0.1),
    c(0, 0, 0)
  )
  
  n <- 10
  gd <- NULL
  for (i in 1:n) {
    smsm <- msm::sim.msm(qmatrix, 14)
    # inspection times uniformly on 2 year intervals until 14 years
    itimes <- seq(2, 24, by=3) + runif(8, 0, 2)
    gdi <- data.frame(id=i, state=evalstep(time=smsm$times, stepf=smsm$states,
                                           newtime=c(0,itimes)), time=c(0,itimes))
    # throw away superfluous duplicate state 3 (absorbing)
    gdi <- gdi[!(gdi$state==3 & duplicated(gdi$state)),]
    gd <- rbind(gd, gdi)
  }
  
  #Create Frydman data
  msmtoFrydman <- function(gd){
    #Create Frydman data
    gd_frydman <- NULL
    for(j in unique(gd$id)){
      tempdat <- subset(gd, id == j)
      tempstates <- unique(tempdat$state)
      if(length(tempstates) == 1){ #If we only observe the subject in 1 state, right censored in 1
        gdi_frydman <- data.frame(delta = 0, Delta = 0,
                                  L = NA,
                                  R = NA,
                                  time = tempdat$time[length(tempdat$time)])
      } else if(length(tempstates) == 2){ #If we only observe the subject in 2 states, either 1->2 or 1->3 has occured
        if(all(tempstates %in% c(1,2))){
          gdi_frydman <- data.frame(delta = 1, Delta = 0, 
                                    L = tempdat$time[which.min(tempdat$state == 1)-1],
                                    R = tempdat$time[which.min(tempdat$state == 1)],
                                    time = tempdat$time[length(tempdat$time)])
        } else if(all(tempstates %in% c(1,3))){
          gdi_frydman <- data.frame(delta = 0, Delta = 1, 
                                    L = NA,
                                    R = NA,
                                    time = tempdat$time[which.min(tempdat$state == 1)])
        }
      } else if(length(tempstates) == 3){ #If we observe 3 states, then 1->2->3 must have occured
        gdi_frydman <- data.frame(delta = 1, Delta = 1, 
                                  L = tempdat$time[which.min(tempdat$state == 1)-1],
                                  R = tempdat$time[which.min(tempdat$state == 1)],
                                  time = tempdat$time[length(tempdat$time)])
      }
      gd_frydman <- rbind(gd_frydman, gdi_frydman)
    }
    return(gd_frydman)
  }
  
  gd_frydman <- msmtoFrydman(gd)
  
  out_frydman <- msm_frydman(gd_frydman, tol = 1e-10)
  
  out_msm <- npmsm(gd, tmat, maxit = 300, tol = 1e-10, exact = c(3))
  out_msm_newmet <- npmsm(gd, tmat, maxit = 300, tol = 1e-10, exact = c(3), newmet = TRUE)
  suppMSM <- support_npmsm(out_msm,  cutoff = 1e-9)
  suppMSM_newmet <- support_npmsm(out_msm_newmet, cutoff = 1e-9)
  
  print("Frydman Support")
  out_frydman$supportMSM$Q_mat
  print("MSM Support")
  suppMSM$`State 1 -> State 2`$support
})





test_that("Multinomial vs Poisson", {
  set.seed(1)
  tmat <- mstate::trans.illdeath()
  #Function to generate evaluation times: at 0 and uniform inter-observation
  eval_times <- function(n_obs, stop_time){
    cumsum( c( 0,  runif( n_obs-1, 0, 2*(stop_time-4)/(n_obs-1) ) ) )
  }
  
  #Simulate illness-death model data with Weibull transitions
  sim_dat <- sim_id_weib(n = 20, n_obs = 6, stop_time = 15, eval_times = eval_times,
                         start_state = "stable", shape = c(0.5, 0.5, 2), scale = c(5, 10, 10/gamma(1.5)))
  
  mod_mult <- npmsm(gd = sim_dat, tmat = tmat, tol = 1e-2)
  mod_pois <- npmsm(gd = sim_dat, tmat = tmat, method = "poisson", tol = 1e-2)
  
  #The two models should not differ too much on the transition probabilities of the first and 
  #second transition
  mult_pt <- transprob(mod_mult, predt = 0)
  pois_pt <- transprob(mod_pois, predt = 0)
  
  #Check closeness of 1->1, 1->2 and 1->3 transitions
  expect_true(all(mult_pt[[1]] - pois_pt[[1]] < 0.2))
  #Check closeness of 2->2 transition
  expect_true(all(mult_pt[[2]][, 1:3] - pois_pt[[2]][, 1:3] < 0.2))
})



test_that("profvis", {
  #Code that should not run for testing.
  skip_if(TRUE)
  # Start with a difficult one right away
  set.seed(3)
  
  #Absorbing states 4 and 6
  qmatrix <- rbind(
    c(-0.6, 0.1, 0.5, 0, 0, 0),
    c(0.08, -0.205, 0, 0, 0.125, 0),
    c(0, 0, -0.25, 0.25, 0, 0),
    c(0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, -0.4, 0.4),
    c(0, 0, 0, 0, 0, 0)
  )
  
  n <- 100
  gd <- NULL
  for (i in 1:n) {
    smsm <- msm::sim.msm(qmatrix, 14)
    # inspection times uniformly on 2 year intervals until 14 years
    itimes <- seq(0, 12, by=2) + runif(7, 0, 2)
    gdi <- data.frame(id=i, state=evalstep(time=smsm$times, stepf=smsm$states,
                                           newtime=c(0,itimes)), time=c(0,itimes))
    # throw away superfluous duplicate 4's ad 6's (absorbing)
    gdi <- gdi[!(gdi$state==4 & duplicated(gdi$state)),]
    gdi <- gdi[!(gdi$state==6 & duplicated(gdi$state)),]
    gd <- rbind(gd, gdi)
  }
  head(gd, n=12)
  tail(gd, n=12)
  
  #Which states can be reached from each state?
  tmat <- mstate::transMat(x = list( c(2, 3), c(1, 5), c(4), c(), c(6), c() ))
  
  
  a <- npmsm(gd = gd, tmat = tmat)
  
  
  #Compare with Frydman, 3 state illness-death model
  set.seed(140703)
  
  #Absorbing states 4 and 6
  qmatrix <- rbind(
    c(-0.2, 0.1, 0.1),
    c(0, -0.1, 0.1),
    c(0, 0, 0)
  )
  
  n <- 100
  gd <- NULL
  for (i in 1:n) {
    smsm <- sim.msm(qmatrix, 14)
    # inspection times uniformly on 2 year intervals until 14 years
    itimes <- seq(0, 12, by=2) + runif(7, 0, 2)
    gdi <- data.frame(id=i, state=evalstep(time=smsm$times, stepf=smsm$states,
                                           newtime=c(0,itimes)), time=c(0,itimes))
    # throw away superfluous duplicate state 3 (absorbing)
    gdi <- gdi[!(gdi$state==3 & duplicated(gdi$state)),]
    gd <- rbind(gd, gdi)
  }
  head(gd, n=12)
  tail(gd, n=12)
  
  #Which states can be reached from each state?
  tmat <- mstate::transMat(x = list( c(2, 3), c(3), c() ))
  
  b <- npmsm(gd = gd, tmat = tmat)
  
  
  #Write code here to perform Frydman EM for above data.
  
  
  #Profiling npmsm
  #library(profvis)
  #Compare with Frydman, 3 state illness-death model
  set.seed(140703)
  
  #Absorbing states 4 and 6
  qmatrix <- rbind(
    c(-0.2, 0.1, 0.1),
    c(0, -0.1, 0.1),
    c(0, 0, 0)
  )
  
  n <- 500
  gd <- NULL
  for (i in 1:n) {
    smsm <- msm::sim.msm(qmatrix, 14)
    # inspection times uniformly on 2 year intervals until 14 years
    itimes <- seq(0, 12, by=2) + runif(7, 0, 2)
    gdi <- data.frame(id=i, state=evalstep(time=smsm$times, stepf=smsm$states,
                                           newtime=c(0,itimes)), time=c(0,itimes))
    # throw away superfluous duplicate state 3 (absorbing)
    gdi <- gdi[!(gdi$state==3 & duplicated(gdi$state)),]
    gd <- rbind(gd, gdi)
  }
  
  #Which states can be reached from each state?
  tmat <- mstate::transMat(x = list( c(2, 3), c(3), c() ))
  
  b <- npmsm(gd = gd, tmat = tmat, tol = 1e-2)
  
  #profvis::profvis({
  #  b <- npmsm(gd = gd, tmat = tmat, tol = 1e-2)
  #})
})



test_that("Turnbull vs MSM (not everyone observed at 0)", {
  
  
}
)





