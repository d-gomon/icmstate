test_that("profiling", {
  skip_if(TRUE)
  
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
  
  
  profvis::profvis({
    b_old <- npmsm(gd = gd, tmat = tmat, tol = 3e-3)
  })
  
  profvis::profvis({
    b_new <- npmsm(gd = gd, tmat = tmat, tol = 3e-3)
  })
  
  
  
  # tictoc::tic("old")
  # b_old <- npmsm(gd = gd, tmat = tmat, tol = 1e-2)
  # toc()
  # 
  # tictoc::tic("new")
  # b_new <- npmsm(gd = gd, tmat = tmat, tol = 1e-2)
  # toc()
  # 
  # tictoc::tic("new rm bins")
  # b_new_rm <- npmsm(gd = gd, tmat = tmat, tol = 1e-2, remove_bins = TRUE)
  # toc()
})

