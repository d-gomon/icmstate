test_that("Right-censoring works as expected", {
  set.seed(1)
  #We want right-censoring to exclude any observations afterwards.
  #Simple 1->2->3 model
  tmat <- mstate::transMat(list(c(2), c(3), c()))
  
  #Set transition rate to state 2 reasonably high, and censoring rate to 3 very high.
  #We should not be able to see state 3 any more
  rc_dat <- sim1_weibmsm(obstimes = seq(0, 20, 1), tmat = tmat, tmat2 = mstate::to.trans2(tmat),
                         shape = c(1, 1), scale = c(3, 1000), censshape = c(1, 1),
                         censscale = c(1000, 0.01))
  
  #Should not contain the state 3
  expect(all(rc_dat[, "state"]!=3), failure_message = "Right-censoring is not working appropriately")
  
  #If we transition into state 2 very quickly, we should not even ever see state 2
  rc_dat2 <- sim1_weibmsm(obstimes = c(0, 5), tmat = tmat, tmat2 = mstate::to.trans2(tmat),
                          shape = c(1, 1), scale = c(0.1, 1000), censshape = c(1, 1),
                          censscale = c(1000, 0.001))
  expect(all(rc_dat[, "state"]!=2), failure_message = "Right-censoring is not working appropriately")
})
