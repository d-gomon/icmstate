
test_that("support_npmsm", {
  skip_if(TRUE)
  
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
  
  gd_eid <- NULL
  for(j in unique(gd$id)){
    tempdat <- subset(gd, id == j)
    tempstates <- unique(tempdat$state)
    #If we observe 3 unique states, 1->2->3 must have happened
    if(length(tempstates) == 3){
      tempdat$state[which(tempdat$state == 3)] <- 4
    }
    gd_eid <- rbind(gd_eid, tempdat)
  }
  
  tmat_eid <- mstate::transMat(x = list( c(2, 3), c(4), c(), c() ))
  
  out_msm_eid <-  npmsm(gd_eid, tmat_eid, maxit = 200, tol = 1e-8, exact = c(3, 4))
})








