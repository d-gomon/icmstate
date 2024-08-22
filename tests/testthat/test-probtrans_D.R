test_that("probtrans_D and probtrans yield same estimates", {
  testthat::skip_if_not_installed("mstate")
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
  
  
  MSMres <- npmsm(gd = gd, tmat = tmat, tol = 1e-3)
  
  
  mstate_probtrans <- probtrans(MSMres$A, predt = MSMres$A$Haz$time[8], direction = "forward", variance = FALSE)
  int_matrices <- get_intensity_matrices(MSMres$A)
  mstate_icmstate <- probtrans_D(int_matrices, predt = MSMres$A$Haz$time[8], direction = "forward", as.df = TRUE)
  
  
  expect_equal(mstate_probtrans[[1]], mstate_icmstate[[1]])
  expect_equal(mstate_probtrans[[2]], mstate_icmstate[[2]])
  expect_equal(mstate_probtrans[[3]], mstate_icmstate[[3]])
  
  
  mstate_probtrans <- probtrans(MSMres$A, predt = MSMres$A$Haz$time[8], direction = "fixedhorizon", variance = FALSE)
  int_matrices <- get_intensity_matrices(MSMres$A)
  mstate_icmstate <- probtrans_D(int_matrices, predt = MSMres$A$Haz$time[8], direction = "fixedhorizon", as.df = TRUE)
  
  
  expect_equal(mstate_probtrans[[1]], mstate_icmstate[[1]])
  expect_equal(mstate_probtrans[[2]], mstate_icmstate[[2]])
  expect_equal(mstate_probtrans[[3]], mstate_icmstate[[3]])
  
  
  mstate_probtrans <- probtrans(MSMres$A, predt = MSMres$A$Haz$time[8] + 0.01, direction = "fixedhorizon", variance = FALSE)
  int_matrices <- get_intensity_matrices(MSMres$A)
  mstate_icmstate <- probtrans_D(int_matrices, predt = MSMres$A$Haz$time[8] + 0.01, direction = "fixedhorizon", as.df = TRUE)
  
  
  expect_equal(mstate_probtrans[[1]], mstate_icmstate[[1]])
  expect_equal(mstate_probtrans[[2]], mstate_icmstate[[2]])
  expect_equal(mstate_probtrans[[3]], mstate_icmstate[[3]])
  
  # nodf is fastest, about 40 times faster than probtrans
  #library(microbenchmark)
  # microbenchmark(
  #   D_df <- probtrans_D(int_matrices, predt = MSMres$A$Haz$time[8] + 0.01, direction = "fixedhorizon", as.df = TRUE),
  #   mstate <- probtrans(MSMres$A, predt = MSMres$A$Haz$time[8] + 0.01, direction = "fixedhorizon", variance = FALSE),
  #   D_nodf <- probtrans_D(int_matrices, predt = MSMres$A$Haz$time[8] + 0.01, direction = "fixedhorizon", as.df = FALSE),
  #   times = 20
  # )
  
  #Test for single entry
  mstate_probtrans <- probtrans(MSMres$A, predt = MSMres$A$Haz$time[1], direction = "fixedhorizon", variance = FALSE)
  int_matrices <- get_intensity_matrices(MSMres$A)
  mstate_icmstate <- probtrans_D(int_matrices, predt = MSMres$A$Haz$time[1], direction = "fixedhorizon", as.df = TRUE)
  
  expect_equal(mstate_probtrans[[1]], mstate_icmstate[[1]])
  expect_equal(mstate_probtrans[[2]], mstate_icmstate[[2]])
  expect_equal(mstate_probtrans[[3]], mstate_icmstate[[3]])
  
  
})





test_that("probtrans_C is consistent with probtrans_D", {
  testthat::skip_if_not_installed("mstate")
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
  
  
  MSMres <- npmsm(gd = gd, tmat = tmat, tol = 1e-3)
  times <- unique(MSMres$A$Haz$time)
  
  
  int_matrices <- get_intensity_matrices(MSMres$A)
  
  #Forward with max range
  prob_D <- probtrans_D(int_matrices, predt = times[1], direction = "forward", as.df = FALSE)
  prob_C <- probtrans_C(int_matrices, predt = times[1], cutoff = max(times), direction = "forward", as.df = FALSE)
  
  expect_equal(prob_D[, , 1], prob_C[, , 1])
  expect_equal(prob_D[, , 2], prob_C[, , 2])
  expect_equal(prob_D[, , 3], prob_C[, , 3])
  
  #Forward with cutoff between two times
  prob_D <- probtrans_D(int_matrices, predt = times[1], direction = "forward", as.df = FALSE)
  prob_C <- probtrans_C(int_matrices, predt = times[1], cutoff = (times[14] + times[15])/2, direction = "forward", as.df = FALSE)
  
  expect_equal(prob_D[1:14, , 1], prob_C[, , 1])
  expect_equal(prob_D[1:14, , 2], prob_C[, , 2])
  expect_equal(prob_D[1:14, , 3], prob_C[, , 3])
  
  
  #Forward starting at 0
  prob_D <- probtrans_D(int_matrices, predt = 0, direction = "forward", as.df = FALSE)
  prob_C <- probtrans_C(int_matrices, predt = 0, cutoff = (times[14] + times[15])/2, direction = "forward", as.df = FALSE)
  
  expect_equal(prob_D[1:15, , 1], prob_C[, , 1])
  expect_equal(prob_D[1:15, , 2], prob_C[, , 2])
  expect_equal(prob_D[1:15, , 3], prob_C[, , 3])
  
  
  #Fixedhorizon with max range
  prob_D <- probtrans_D(int_matrices, predt = max(times), direction = "fixedhorizon", as.df = FALSE)
  prob_C <- probtrans_C(int_matrices, predt = max(times), cutoff = 0, direction = "fixedhorizon", as.df = FALSE)
  
  expect_equal(prob_D[, , 1], prob_C[, , 1])
  expect_equal(prob_D[, , 2], prob_C[, , 2])
  expect_equal(prob_D[, , 3], prob_C[, , 3])  
  
  #Fixedhorizon with cutoff at a time
  prob_D <- probtrans_D(int_matrices, predt = max(times), direction = "fixedhorizon", as.df = FALSE)
  prob_C <- probtrans_C(int_matrices, predt = max(times), cutoff = times[14], direction = "fixedhorizon", as.df = FALSE)
  
  expect_equal(prob_D[15:(length(times)+1), , 1], prob_C[, , 1])
  expect_equal(prob_D[15:(length(times)+1), , 2], prob_C[, , 2])
  expect_equal(prob_D[15:(length(times)+1), , 3], prob_C[, , 3])
  
  #Fixedhorizon with cutoff at 0
  prob_D <- probtrans_D(int_matrices, predt = max(times), direction = "fixedhorizon", as.df = FALSE)
  prob_C <- probtrans_C(int_matrices, predt = max(times), cutoff = 0, direction = "fixedhorizon", as.df = FALSE)
  
  expect_equal(prob_D[, , 1], prob_C[, , 1])
  expect_equal(prob_D[, , 2], prob_C[, , 2])
  expect_equal(prob_D[, , 3], prob_C[, , 3])
})

