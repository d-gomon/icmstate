## Install development version of icmstate
#library(devtools)
# install_github("d-gomon/icmstate",
#               ref = "a0c61ae9da0a036a6d47e3a04ed75f4916e5597c",
#               auth_token = "DONTTELLANYONE:)")

## --------------------------------------------------------------------------------------------------------------
library(doParallel)
library(doRNG)
library(msm)
library(icmstate)

#-------------------------------------------------------------------------
#Read in arguments from command line
#-------------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)
#Arguments from command line are (in this order):
#scenario
#n
#n_obs
#N
#method (Poisson/multinomial/MSM)
#RNG (integer for RNG in doRNG pkg)
print(args)
scenario <- as.numeric(args[1])
n <- as.numeric(args[2])
n_obs <- as.numeric(args[3])
N <- as.numeric(args[4])
method <- args[5]
RNG_seed <- as.numeric(args[6])



#-------------------------------------------------------------------------
# General Simulation Parameters
#-------------------------------------------------------------------------
C_admin <- 15
params <- list()
params$method <- method

n.cores <- as.numeric(system("nproc", intern = TRUE))
stop_time <- 15



eval_times <- function(n_obs, stop_time){
  cumsum( c( 0,  runif( n_obs-1, 0, 2*(stop_time-4)/(n_obs-1) ) ) )
}


if(scenario != 4){
  generate_gd <- function(haz_type = c("hom", "weib"), start_type = c("stable", "equalprob"), qmatrix = NULL, ...){
    #Exponential hazard
    if(haz_type == "hom"){
      gd <- sim_id_weib(n = n, n_obs = n_obs, stop_time = stop_time, eval_times = eval_times,
                        start_state = start_type, shape = c(1, 1, 1), scale = 1/exp_rates)
    } else if(haz_type == "weib"){ #Weibull hazard
      gd <- sim_id_weib(n = n, n_obs = n_obs, stop_time = stop_time, eval_times = eval_times,
                        start_state = start_type, shape = w_shapes, scale = w_scales)
    }
    
    #Censor at appropriate time
    which_censored <- which(gd[, "time"] > C_admin)
    if(length(which_censored) > 0){
      gd <- gd[-which_censored,]
    }
    return(gd)
  }
} else{
  generate_gd <- function(haz_type = c("hom", "weib"), start_type = c("stable", "equalprob"), qmatrix = NULL, ...){
    if(start_type == "stable"){
      start_states <- rep(1, n)
    } else{
      start_states <- sample.int(2, size = n, replace = TRUE, prob = c(0.5, 0.5))
    }
    
    #We make eval_times so that the mean of the sum of n_obs - 1 uniform variables is equal to stop_time - 4 (time 11)
    eval_times <- function(n_obs, stop_time){
      cumsum( c( 0,  runif( n_obs-1, 0, 2*(stop_time-4)/(n_obs-1) ) ) )
    }
    
    sim.df <- data.frame(subject = rep(1:n, rep(n_obs, n)), 
                         time = c(replicate(n, eval_times(n_obs, stop_time = stop_time))))
    which_censored <- which(sim.df$time > C_admin)
    if(length(which_censored) > 0){
      sim.df <- sim.df[-which_censored,]
    }
    sim.df <- sim.df[order(sim.df$subject, sim.df$time),]
    gd <- msm::simmulti.msm(data = sim.df, qmatrix = qmatrix, start = start_states, death = c(3,4))
    gd <- gd[, c(1:3)]
    colnames(gd) <- c("id", "time", "state")
    
    return(gd)
  }
}









#-------------------------------------------------------------------------
# Scenario specific simulation parameters
#-------------------------------------------------------------------------

if( scenario == 1){
  
  qmatrix <- rbind(
    c(-0.15, 0.1, 0.05),
    c(0, -0.1, 0.1),
    c(0, 0, 0)
  )
  tmat <- mstate::transMat(list(c(2, 3), c(3), c()))
  params$tmat <- tmat
  params$method <- method

  haz_type <- "hom"
  start_type <- "stable"
  exp_rates <- c(0.1, 0.05, 0.1)
  
} else if( scenario == 2){
  
  qmatrix <- rbind(
    c(-0.15, 0.1, 0.05),
    c(0, -0.1, 0.1),
    c(0, 0, 0)
  )
  tmat <- mstate::transMat(list(c(2, 3), c(3), c()))
  params$tmat <- tmat
  params$method <- method

  haz_type <- "hom"
  start_type <- "equalprob"
  exp_rates <- c(0.1, 0.05, 0.1)
  
} else if( scenario == 3){
  
  #To-Do: generate Weibull hazards here (possible with MSM?)
  qmatrix <- rbind(
    c(-0.15, 0.1, 0.05),
    c(0, -0.1, 0.1),
    c(0, 0, 0)
  )
  tmat <- mstate::transMat(list(c(2, 3), c(3), c()))
  params$method <- method
  params$tmat <- tmat
  
  haz_type <- "weib"
  start_type <- "stable"
  w_shapes <- c(0.5, 0.5, 2)
  w_scales <- c(5, 10, 10/gamma(1.5))
  
} else if( scenario == 4){
  
  qmatrix <- rbind(
    c(-0.15, 0.1, 0.05, 0),
    c(0, -0.1, 0, 0.1),
    c(0, 0, 0, 0),
    c(0, 0, 0, 0)
  )
  tmat <- mstate::transMat(list(c(2, 3), c(4), c(), c()))
  params$tmat <- tmat
  params$method <- method
  params$exact <- c(3, 4)

  haz_type <- "hom"
  start_type <- "stable"
  exp_rates <- c(0.1, 0.05, 0.1)
  
}



#-------------------------------------------------------------------------
# Run models
#-------------------------------------------------------------------------

out.name <- paste0("sc", scenario, "n", n, "obs", n_obs, "N", N, method, "RNG", RNG_seed)

cl <- makeCluster(n.cores, outfile = paste0(out.name, ".txt"))
registerDoParallel(cl)
out <-  foreach(i = 1:N, .packages = c("icmstate", "msm"),
                .options.RNG = RNG_seed) %dorng% {
  print(i)
  print(Sys.time())
  # Generate Data (either using homogeneous or Weibull hazards)
  gd <- generate_gd(haz_type = haz_type, start_type = start_type, qmatrix = qmatrix)
  gd <- remove_redundant_observations(gd = gd, tmat = tmat)
  # Parameters for running the algorithm
  params$gd <- gd
  params$maxit <- 1200
  params$checkMLE_tol <- 1e-3
  params$conv_crit <- "haz"
  params$tol <- 1e-4
  
  if(params$method == "msm"){
    out <- tryCatch({msm::msm(state ~ time, subject = id, data = gd, qmatrix = qmatrix,
                              gen.inits = TRUE, deathexact = params$exact)},
                    error = function(cond){
                      as.character(cond)
                    })
    if(inherits(out, "msm")){
      out$data <- gd
    }
    #We can specify pci = c(times at which piecewise changes) (e.g. pci = c(1, 3, 5))
    #Gu et al. do this, but I see not reason to, as MSM is currently correctly defined.
  } else{
    out <- tryCatch({do.call(npmsm, params)},
                    error = function(cond){
                      as.character(cond)
                    })
  }
  
  out
}
stopCluster(cl) 

assign(out.name, out)

save(list = out.name, file = paste0(out.name, ".Rdata"))




