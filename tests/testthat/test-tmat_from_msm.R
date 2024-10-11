test_that("qmat_from_tmat", {
  #Simple survival
  tmat <- mstate::transMat(x = list( c(2), c() ))
  qmat <- qmat_from_tmat(tmat)
  e_res_simple <- rbind(c(0, 0.1), c(0, 0))
  expect_true(all(e_res_simple == qmat))
  
  #Illness-death
  tmat <- mstate::transMat(x = list(c(2, 3), c(3), c()))
  qmat <- qmat_from_tmat(tmat)
  e_res_id <- rbind(c(0, 0.1, 0.1), c(0, 0, 0.1),
                    c(0, 0, 0))
  expect_true(all(e_res_id == qmat))
})


test_that("tmat_from_msm", {
  #Extract tmat from ID msm
  eval_times <- function(n_obs, stop_time){
    cumsum( c( 0,  runif( n_obs-1, 0, 2*(stop_time-4)/(n_obs-1) ) ) )
  }
  set.seed(2)
  sim_dat <- sim_id_weib(n = 5, n_obs = 6, stop_time = 15, eval_times = eval_times,
                         start_state = "stable", shape = c(1, 1, 1), scale = c(5, 10, 1))
  
  tmat <- mstate::trans.illdeath()
  sim_dat <- remove_redundant_observations(sim_dat, tmat)
  
  msm_mod <- msm::msm(state ~ time, subject = id, data = sim_dat,
                      qmatrix = qmat_from_tmat(tmat), gen.inits = TRUE)
  tmat_msm <- tmat_from_msm(msm_mod)
  
  expect_true(all(tmat_msm == tmat, na.rm = TRUE))
})
