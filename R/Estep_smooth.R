
#' E-step of the EM algorithm for smooth estimation of transition intensities
#' 
#' 
#' 
#' @importFrom JOPS bbase
#' 
#' @keywords internal
#' 
#' 
#' 
#' 
#' 
#' 

#We want to calculate the at-risk matrix and the expect trans matrix
#for each bin for each transition for each subject (row, col, 3rd dim)

Estep_smooth <- function(fix_pars, subject_slices, EM_est, it_num){
  #Inputs:
  #tmat(2): to know the transition + numbers
  #model.matrix: to know the covariates of subjects
  #n_subjects: to know over how many subjects to iterrate
  #n_bins: to know how many bins we have
  #function (do we want P-splines? Constant hazards?)
  #gd: we need the transition times
  #subject_slices: to know what to slice over
  #bin_size: to determine expected number of transitions
  
  #Extract the required parameters
  
  n_subjects <- fix_pars[["n_subjects"]]
  n_states <- fix_pars[["n_states"]]
  n_transitions <- fix_pars[["n_transitions"]]
  n_bins <- fix_pars[["n_bins"]]
  n_covariates <- fix_pars[["n_covariates"]]
  n_splines <- fix_pars[["n_splines"]]
  tmat2 <- fix_pars[["tmat2"]]
  deg_splines <- fix_pars[["deg_splines"]]
  n_segments <- fix_pars[["n_segments"]]
  n_coefficients <- fix_pars[["n_coefficients"]]
  use_RA <- fix_pars[["use_RA"]]
  gd <- fix_pars[["gd"]]
  mod_matrix <- fix_pars[["mod_matrix"]]
  
  #Initial value for ode solver - maybe also push to fix_pars?
  ode_init <- diag(1, n_states, n_states)
  
  
  # Reset NumTrans and AtRisk to zero arrays.
  EM_est[["NumTrans"]][] <- 0 #array(0, dim=c(n_bins, n_transitions, n_subjects))
  EM_est[["AtRisk"]][] <- 0 #array(0, dim=c(n_bins, n_states, n_subjects))
  
  # Initialize log-likelihood
  loglik <- 0
  
  # Loop over subjects
  for (i in 1:n_subjects)
  {
    # deb(i, method="cat")
    #We use subject_slices (hashmap) to quickly slice gd 
    #(probably faster if we use data.table, but oh well)
    datai <- gd[subject_slices[[i]], ]
    nobs <- nrow(datai) - 1
    if (nobs > 0) {
      #Loop over observations of subject i
      for (iobs in 1:nobs) {
        # deb(iobs, method="cat")
        # Identify starting and ending times and states
        #browser()
        g <- datai[iobs, "state"] # starting state
        s <- datai[iobs, "time"] # starting time
        h <- datai[iobs + 1, "state"] # ending state
        t <- datai[iobs + 1, "time"] # ending time
        sbin <- ceiling(datai[iobs, "time"])
        tbin <- ceiling(datai[iobs + 1, "time"])
        #st <- t - s # interval length
        # deb(c(iobs, g, h, s, t, st), method="cat")
        # Within this interval reset time s to 0, time t to st
        # ntimes <- ceiling(st - 0.0001)
        
        #Nothing to do if both observations at the same time
        if (tbin == sbin) {
          break
        } else if (tbin == sbin + 1) {
          # In this case we either say (g=h) that no event has taken place,
          # (Y[tbin, , i]) = 0, which was already the case, so no need to explicitly code,
          # or (g != h), Y[tbin, g->h, i] = 1
          #
          # Assumption behind this is that bins are so close that at most one transition
          # could have happened between two bins
          if (g != h) {
            # Find out which transition
            transno <- which(tmat2$from == g & tmat2$to == h)
            EM_est[["NumTrans"]][tbin, transno, i] <- 1
          }
          # In both cases we have R[tbin, h, i] = 1 and all other R[tbin, , i] remain 0
          EM_est[["AtRisk"]][tbin, h, i] <- 1
        } else {
          # times <- c(0, (sbin+1):(tbin-1) - s, st)
          times <- sbin : tbin
          revtimes <- tbin : sbin
          ntimes <- tbin - sbin + 1L
          idx <- (sbin+1) : tbin # bins to be filled
          
          #
          # TODO: new structure, separate h=censor from otherwise
          #
          # TODO: CHECK in case of absorbing transition, can sum(Y)>1???
          #
          
          # deb(c(sbin, tbin, times), method="cat")
          
          
          # parms[[5]] <- sbin
          # Forward ODE
          # cat("Forward ODE ...\n")
          #browser()
          fwd <- suppressWarnings(ode(y = ode_init, times = times, func = ChapKolm_fwd_smooth,
                     parms = EM_est, method = fix_pars[["ode_solver"]], fix_pars = fix_pars, 
                     subject = i))
          fwd[fwd < 0] <- 0
          # fwd_array <- array(fwd[, -1], dim=c(ntimes1, n_states, n_states))
          fwd_array <- array(fwd[, -1], dim=c(ntimes, n_states, n_states))
          # Backward ODE
          # cat("Backward ODE ...\n")
          bwd <- suppressWarnings(ode(y = ode_init, times = revtimes, func = ChapKolm_bwd_smooth,
                     parms = EM_est, method = fix_pars[["ode_solver"]], fix_pars = fix_pars,
                     subject = i))
          bwd[bwd < 0] <- 0
          # bwd_array <- array(bwd[ntimes1:1, -1], dim=c(ntimes1, n_states, n_states)) # backwards in time
          
          bwd_array <- array(bwd[ntimes:1, -1], dim=c(ntimes, n_states, n_states)) # backwards in time
          # Combine
          fwdg <- fwd_array[, g, ]
          #if (h == censor)
          #  bwdh <- apply(bwd_array[, , censor.states], c(1, 2), sum)
          #else bwdh <- bwd_array[, , h]
          bwdh <- bwd_array[, , h]
          fwd_bwd <- fwdg * bwdh
          # Standardize
          #if (h == censor)
          #  denom <- sum(bwd_array[1, g, censor.states])
          #else denom <- bwd_array[1, g, h]
          denom <- bwd_array[1, g, h]
          #if (is.na(log(denom))) deb(c(i, iobs), method="cat")
          P_cond <- fwd_bwd / denom
          # AtRisk[sbin:tbin, , i] <- P_cond
          #Used to be P_cond[-ntimes, , drop = FALSE]
          EM_est[["AtRisk"]][idx, , i] <- P_cond[-1, , drop=FALSE]
          if(denom == 0){ #With our current estimates, observed path is impossible
            EM_est[["AtRisk"]][idx, , i] <- 0 #So subject is not at risk of any of the transitions
            warning(paste0("In E-step of iteration", it_num, "Observed path (", g, "->", h, ") for subject ", i, " between times ", 
                           s * fix_pars[["bin_length"]], " and ", t * fix_pars[["bin_length"]], 
                           " not possible with current estimates. Setting estimates to 0."))
          }
          
          for (transno in 1:n_transitions) {
            from <- tmat2$from[transno]
            to <- tmat2$to[transno]
            s2t <- (sbin : tbin)
            #Calculate the log-hazard
            loghaz <- c(JOPS::bbase(s2t, xl = 0, xr = n_bins, nseg = n_segments, 
                                      bdeg = deg_splines) %*% EM_est[["coeff_old"]][1:n_splines, transno])
            #Risk-adjust if we have to
            if(use_RA){
              loghaz_RA <- c(mod_matrix[i, ] %*% EM_est[["coeff_old"]][(n_splines+1):n_coefficients, transno])
            } else{
              loghaz_RA <- 0
            }
            loghaz <- loghaz + loghaz_RA
            # EM_est[["NumTrans"]][sbin:tbin, transno, i] <- fwdg[, from] * exp(loghaz) * bwdh[, to] / denom
            tmp <- fwdg[, from] * exp(loghaz) * bwdh[, to] / denom
            tmp2 <- (tmp[-1] + tmp[-ntimes])/2
            #If we have issues with the observed path using current estimates, set contribution to NumTrans to zero.
            if(any(is.na(tmp2))){
              tmp2[which(is.na(tmp2))] <- 0
            }
            EM_est[["NumTrans"]][idx, transno, i] <- tmp2
            #Check if this subject is not at risk at any times we are considering
            #Then simply set the expected Number of Transitions there to 0
            not_at_risk <- which(EM_est[["AtRisk"]][idx, from, i] == 0)
            if (length(not_at_risk) > 0) {
              if (any(EM_est[["NumTrans"]][(idx)[not_at_risk], transno, i] > 0)) {
                # deb(c(i, iobs, sbin, tbin, not_at_risk), method="cat")
                EM_est[["NumTrans"]][(idx)[not_at_risk], transno, i] <- 0
              }
            }
          }
          loglik <- loglik + log(denom)
        }
        #browser()
      }
    }
  }
  
  # What about this one: should never be allowed to be > 1
  
  EM_est[["NumTrans"]][EM_est[["NumTrans"]]>1] <- 1
  EM_est[["loglik_new"]] <- loglik
  return(EM_est)
}




