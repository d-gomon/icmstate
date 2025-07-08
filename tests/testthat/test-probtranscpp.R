test_that("get_intensity_matrices_cpp same as get_intensity_matrices and probtrans_D_cpp same
          as probtrans_D", {
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
  gd <- sim_weibmsm(tmat = tmat, shape = c(1, 1, 1), scale = c(10, 20, 10),
                    n_subj = n, obs_pars = c(3, 0.5, 24))
  
  
  MSMres <- npmsm(gd = gd, tmat = tmat, tol = 1e-3)
  times <- unique(MSMres$A$Haz$time)
  
  
  int_matrices <- get_intensity_matrices(MSMres$A)
  int_matrices_cpp <- get_intensity_matrices_cpp(MSMres$A)
  
  expect_equal(int_matrices, int_matrices_cpp)
  
  #library(microbenchmark)
  #microbenchmark(int_matrices = get_intensity_matrices(MSMres$A),
  #               int_matrices_cpp = get_intensity_matrices_cpp(MSMres$A),
  #               times = 1000)
  #It's really much faster
  
  #Forward with max range
  prob_D <- probtrans_D(int_matrices, predt = times[1], direction = "forward", as.df = FALSE)
  prob_D_cpp <- probtrans_D_cpp(int_matrices, predt = times[1], direction = "forward", as_df = FALSE)
  
  expect_equal(prob_D, prob_D_cpp)
  
  #library(microbenchmark)
  #microbenchmark(prob_D = probtrans_D(int_matrices, predt = times[1], direction = "forward", as.df = FALSE),
  #              prob_D_cpp = probtrans_D_cpp(int_matrices, predt = times[1], direction = "forward", as_df = FALSE),
  #              times = 1000)
  #Not much faster, because we don't request many times
  
  
})


test_that("larger n, check if cpp faster", {
  testthat::skip_if(TRUE)
  tmat <- mstate::transMat(x = list( c(2, 3), c(3), c() ))
  set.seed(1)
  #Absorbing states 4 and 6
  qmatrix <- rbind(
    c(-0.15, 0.1, 0.05),
    c(0, -0.1, 0.1),
    c(0, 0, 0)
  )
  
  n <- 40
  gd <- sim_weibmsm(tmat = tmat, shape = c(1, 1, 1), scale = c(10, 20, 10),
                    n_subj = n, obs_pars = c(3, 0.5, 24))
  
  
  MSMres <- npmsm(gd = gd, tmat = tmat, tol = 1e-3)
  times <- unique(MSMres$A$Haz$time)
  
  
  int_matrices <- get_intensity_matrices(MSMres$A)
  int_matrices_cpp <- get_intensity_matrices_cpp(MSMres$A)
  
  expect_equal(int_matrices, int_matrices_cpp)
  
  #library(microbenchmark)
  #microbenchmark(int_matrices = get_intensity_matrices(MSMres$A),
  #               int_matrices_cpp = get_intensity_matrices_cpp(MSMres$A),
  #               times = 1000)
  #It's really much faster
  
  #Forward with max range
  prob_D <- probtrans_D(int_matrices, predt = times[1], direction = "forward", as.df = FALSE)
  prob_D_cpp <- probtrans_D_cpp(int_matrices, predt = times[1], direction = "forward", as_df = FALSE)
  
  expect_equal(prob_D, prob_D_cpp)
  
  #library(microbenchmark)
  #microbenchmark(prob_D = probtrans_D(int_matrices, predt = times[1], direction = "forward", as.df = FALSE),
  #             prob_D_cpp = probtrans_D_cpp(int_matrices, predt = times[1], direction = "forward", as_df = FALSE),
  #             times = 1000)
  #Still not much faster, so using C++ doesn't gain us anything here. Probably mostly due to the data transformations
  #back and forth.

  
})



