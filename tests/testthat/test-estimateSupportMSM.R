test_that("estimate_support_msm", {
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
  
  
  tmat <- mstate::transMat(x = list( c(2, 3), c(1, 5), c(4), c(), c(6), c() ))
  tmat2 <- mstate::to.trans2(tmat)
})
