#' Helper function for npmsm()
#' 
#' @description For a general Markov chain multi-state model with interval censored 
#' transitions calculate the NPMLE using an EM algorithm with multinomial approach
#' 
#' 
#' @param gd A \code{data.frame} with the following named columns
#'\describe{
#'   \item{\code{id}:}{Subject idenitifier;}
#'   \item{\code{state}:}{State at which the subject is observed at \code{time};}
#'   \item{\code{time}:}{Time at which the subject is observed;}
#' } The true transition time between states is then interval censored between the times.
#' @param tmat A transition matrix as created by \code{transMat}
#' @param exact Numeric vector indicating to which states transitions are observed at exact times.
#' Must coincide with the column number in \code{tmat}.
#' @param maxit Maximum number of iterations.
#' @param tol Tolerance of the procedure.
#' @param manual Manually specify starting transition intensities?
#' @param newmet Should contributions after last observation time also be used 
#' in the likelihood? Default is FALSE.
#' @param include_inf Should an additional bin from the largest observed time to 
#' infinity be included in the algorithm? Default is FALSE.
#' @param checkMLE Should a check be performed whether the estimate has converged 
#' towards a true Maximum Likelihood Estimate? Default is TRUE.
#' @param checkMLE_tol Tolerance for checking whether the estimate has converged to MLE.
#' Whenever an estimated transition intensity is smaller than the tolerance, it is assumed 
#' to be zero.
#' @param ... Not used yet
#' 
#' @inheritParams npmsm
#' 
#' @importFrom mstate to.trans2 msfit probtrans
#' @importFrom msm crudeinits.msm
#' @importFrom stats runif rbeta
#' @keywords internal
#' 
#' @references Michael G. Hudgens, On Nonparametric Maximum Likelihood Estimation with 
#' Interval Censoring and Left Truncation, Journal of the Royal Statistical Society 
#' Series B: Statistical Methodology, Volume 67, Issue 4, September 2005, Pages 573-587,
#'  \doi{10.1111/j.1467-9868.2005.00516.x}
#' 
#' 

EM_multinomial <- function(gd, tmat, tmat2, inits, beta_params, support_manual, exact, maxit, tol, conv_crit, manual, 
                        verbose, newmet, include_inf, checkMLE, checkMLE_tol, prob_tol, remove_bins,
                        init_int = init_int, ...){
  # Remove notes from CRAN check
  to <- id <- time <- NULL
  
  #Data transformations for probtrans compatibility
  #We set the smallest time to 0, because probtrans always calculates from 0 onward.
  min_time <- min(gd[, "time"])
  gd[which(gd$time == min_time), "time"] <- 0
  
  
  # Indices + General Information -------------------------------------------
  
  H <- nrow(tmat) # no of states
  M <- nrow(tmat2) # no of transitions
  unique_id <- unique(gd[, "id"])
  n <- length(unique_id)
  gd_original <- gd
  final_step <- FALSE #We are not in our final iteration yet
  ll_storage <- c()
  KKT_violated <- 0
  
  
  

  
  
  ### Unique time points
  taus <- sort(unique(gd$time))
  taus <- taus[-1] # without tau=0
  K <- length(taus)
  if(include_inf){
    taus <- c(taus, Inf)  
  }
  #Check if we will have any numerical problems
  if(any(diff(taus) <= 10*.Machine$double.eps)){
    stop("Time differences between (some) transitions too small, please scale 
         time so that any(diff(gd[['time']]) > 10*.Machine$double.eps)")
  }
  
  
  
  # Exactly Observed States Pre-processing ----------------------------------
  
  if(!missing(exact)){
    exact <- sort(exact)
    H_exact <- length(exact)
    #From which states is the exactly observed state reachable?
    reachable_from <- vector(mode = "list", length = H_exact)
    #What is the associated transition number?
    reachable_from_transno <- vector(mode = "list", length = H_exact)
    for(q in 1:H_exact){
      reachable_from[[q]] <- subset(tmat2, to == exact[q])$from
      reachable_from_transno[[q]] <- subset(tmat2, to == exact[q])$transno
    }
    
    #If we want to use exact transitions, we need to add extra infinitesimal bins
    #where we calculate the probability of instantaneous transitions.
    #In the code below, we duplicate the rows in gd where an exact transition took place
    #and subtract 10 times the machine tolerance from the time
    #This creates 2 numerically different times: one for the initial state and one 
    #for the state the MSM transitions to.
    #This way we can use probtrans, as it requires times to be unique (for the machine)
    exact_times <- numeric(0)
    gd_exact <- NULL
    for(i in 1:n){
      #Consider only one subject
      gdi <- subset(gd, id == unique_id[i])
      
      #Which times does the transition happen at?
      trans_bool <- (c(0, diff(gdi$state)) != 0) & (gdi$state %in% exact)
      trans_idx <- which(trans_bool != 0)
      #Duplicate corresponding rows in df
      df_dupl_idx <- rep(1:nrow(gdi), times = trans_bool+1)
      gdi_exact <- gdi[df_dupl_idx, ]
      gdi_exact$transidx <- 0
      
      #Down the time of duplicated entries (creating infitesimal time intervals), if there are any:
      idx_duplicated <- which(duplicated(gdi_exact$time, fromLast = TRUE))
      if(length(idx_duplicated) != 0){
        gdi_exact[idx_duplicated,]$time <- gdi_exact[idx_duplicated, ]$time - 10*.Machine$double.eps
        gdi_exact[idx_duplicated,]$state <- gdi_exact[(idx_duplicated - 1), ]$state
        gdi_exact[idx_duplicated,]$transidx <- 1
      }
      #Bind together to create final df
      gd_exact <- rbind(gd_exact, gdi_exact)
    }
    gd <- gd_exact
    taus <- sort(unique(gd$time))
    taus <- taus[-1]
    K <- length(taus)
  }
  
  
  # Initialize Initial Cumulative Intensity (manual or automatic) -----------
  
  if(!missing(support_manual)){
    #If support has manually been determined, we would like to give initial mass
    #to transitions only on the support regions.
    warning("Estimated support sets used for transitions. Note that 
            resulting estimates may differ from estimateSupport = FALSE. Check 
            manually if in doubt.")
    Haz_manual <- c()
    for(m in 1:M){ #For each transition
      which_taus <- c()
      #Determine which taus are contained in the support region (we assign mass only here)
      for(s in 1:nrow(support_manual[[m]])){
        which_taus <- c(which_taus, which(taus > support_manual[[m]][s, 1] & taus <= support_manual[[m]][s, 2]))
      }
      #Give equal mass to every tau in the support set
      Haz_m <- numeric(length = K)
      Haz_m[which_taus] <- 1/K
      Haz_m <- cumsum(Haz_m)
      Haz_manual <- c(Haz_manual, Haz_m)
    }
    A <- list(Haz = data.frame(time=rep(taus, M), Haz=Haz_manual,
                               trans=rep(1:M, rep(K,M))),
              trans = tmat)
    attr(A, "class") <- "msfit"
  } else if(manual){
    print(paste0("Enter your initial estimates for the cumulative hazard as a vector of length ", M*K, " at the following times: ", paste(c(0, taus[-length(taus)]), taus, sep = "-")))
    Haz_manual <- scan(what = double(), nmax = M*K)
    print(Haz_manual)
    #Create list with initial hazard estimates and corresponding transition numbering
    A <- list(Haz = data.frame(time=rep(taus, M), Haz=Haz_manual,
                               trans=rep(1:M, rep(K,M))),
              trans = tmat)
    attr(A, "class") <- "msfit"
  } else if(all(init_int != c(0, 0))){
    percen_bins <- floor(K * init_int[2])
    mass_bins <- init_int[1]
    
    A <- list(Haz = data.frame(time=rep(taus, M), 
                               Haz=rep(c(seq(0, mass_bins, length.out = percen_bins), seq(mass_bins, 1, length.out = K - percen_bins)), M),
                               trans=rep(1:M, rep(K,M))),
              trans = tmat)
    attr(A, "class") <- "msfit"
  } else{ 
    #Randomly generate initial estimates from distribution given in inits.
    
    if(inits == "equalprob"){
      #Create list with initial hazard estimates and corresponding transition numbering
      A <- list(Haz = data.frame(time=rep(taus, M), Haz=rep((1:K)/K,M),
                                 trans=rep(1:M, rep(K,M))),
                trans = tmat)
    } else if(inits == "homogeneous"){
      crudeinits <- msm::crudeinits.msm(state ~ time, subject = id, qmatrix = qmat_from_tmat(tmat),
                          data = gd)
      #Estimates will now be time diff * estimate A(t) = \lambda * t
      cumhaz <- numeric(length = K * M)
      time_diffs <- c(taus[1], diff(taus))
      for(m in 1:M){
        #For each transition, calculate cumulative hazard values
        cumhaz[(m-1)*K + 1:K] <- cumsum(pmin(crudeinits[tmat2[m, c("from")], tmat2[m, c("to")]] * time_diffs, 0.99))
      }
      A <- list(Haz = data.frame(time=rep(taus, M), Haz=cumhaz,
                                 trans=rep(1:M, rep(K,M))),
                trans = tmat)
    } else if(inits == "unif"){
      A <- list(Haz = data.frame(time=rep(taus, M), Haz = c(replicate(M, cumsum(runif(K, 0, 1)))),
                                 trans = rep(1:M, rep(K, M))),
                trans = tmat)
    } else if(inits == "beta"){
      A <- list(Haz = data.frame(time = rep(taus, M), Haz = c(replicate(M, cumsum(rbeta(K, shape1 = beta_params[1], shape2 = beta_params[2])))),
                                 trans = rep(1:M, rep(K, M))),
                trans = tmat)
    }
    attr(A, "class") <- "msfit"
  }
  Ainit <- A
  

  
  # EM Algorithm - Start -----------------------------------------------------
  
  #Y_{gi}^{(k)} is defined on each of the:
  #states g (# = H)
  #persons i (# = n)
  #intervals k (# = K)
  Y <- array(0, dim=c(n, K, H))
  
  #d_{ghi}^{(k)} is defined on each of the:
  #persons i (# = n)
  #intervals k (# = K)
  #transitions m (# = M)
  d <- array(0, dim=c(n, K, M))
  
  # ---------Initiate EM values---------#
  llold <- -Inf
  d_old <- d
  Y_old <- Y
  
  
  # Main Loop (over iterations) -----------------------------------------------
  
  
  for (it in 1:(maxit+1)) { # Iterations
    nan_warning_issued <- FALSE
    if(it == (maxit+1)) {
      final_step <- TRUE
    }
    ll <- 0 #Start with 0 value for likelihood
    int_mats <- get_intensity_matrices(A)
    
    
    # Removing taus when they have 0 mass - IN PROGRESS IN PROGRESS
    #browser()
    if(remove_bins){
      can_be_removed <- which_identity_arr(int_mats$intensity_matrices, 3, H)
      which_can_be_removed <- which(can_be_removed)
      consecutive <- c(diff(which_can_be_removed) == 1, FALSE)
      remove <- which_can_be_removed[consecutive]
      if(length(remove) > 0){
        int_mats$intensity_matrices <- int_mats$intensity_matrices[, , -remove, drop = FALSE]
        int_mats$unique_times <- int_mats$unique_times[-remove]
        taus <- taus[-remove]
        K <- length(taus)
        Y <- array(0, dim=c(n, K, H))
        d <- array(0, dim=c(n, K, M))
      }  
    }
    #browser()
    
    #* Loop over subjects ------------------------------------------------------
    
    for (i in 1:n) { # all subjects
      gdi <- subset(gd_original, id==unique_id[i])
      ni <- nrow(gdi)-1
      
      #** Loop over subject times -------------------------------------------------
      
      for (j in 1:ni) {      # Go through all observation intervals of subject i
        li <- gdi$time[j]    # Left bound of current interval
        ri <- gdi$time[j+1]  # Right bound of current interval
        ai <- gdi$state[j]   # State at left bound of interval
        bi <- gdi$state[j+1] # State at right bound of interval
        
        
        #Removing bins if mass 0 - IN PROGRESS
        if(remove_bins){
          #If left end point was removed, set it to first one smaller
          if(!(li %in% taus) & !(li < taus[1])){
            li <- taus[binary_search_larger_equal(taus, li)]
          }
          #If right end point was removed, set it to first one larger
          if(!(ri %in% taus)){
            ri <- taus[binary_search_larger_equal(taus, ri) + 1]
          }  
        }
        
        
        
        wh <- which(taus<=ri & taus>li)  # 1{\tau \in (l_i, r_i]}
        
        
        
        #*** Case: Only Interval Censored --------------------------------------------
        
        if( missing(exact) || !(bi %in% exact) ) { 
          
          #**** Forward Probabilities ---------------------------------------------------
          ptf <- probtrans_C(int_mats, predt=li, cutoff = ri, direction="forward", as.df=FALSE) #Calculate P(li, t) - forward multiplication
          #Cut-off at ri, because we don't need probabilities after ri
          # ptf <- tryCatch(probtrans(A, predt=li, direction="forward", variance=FALSE),
          #                 warning = function(cond){
          #                   message(conditionMessage(cond))
          #                   print(tail(probtrans(A, predt=li, direction="forward", variance=FALSE)[[1]]))
          #                   probtrans(A, predt=li, direction="forward", variance=FALSE)
          #                 })
          ptf <- ptf[, , ai] #Extract only the probabilities from the state we are in
          #ptf is now P_{a_i *}(l_i, t) with t variable (rows) and * the states we can transition to (columns)
          denom <- ptf[nrow(ptf), bi+1] #Extract P_{ai, bi}(l_i, r_i)
          if(denom == 0 & !nan_warning_issued){
            warning(paste0("Transition ", ai, " to ", bi, " between times ", li, " and ", ri, " for subject ", i, " in iteration ", it, " is impossible with current estimates.
                        An estimated probability was 0, whereas it shouldn't have been. Try increasing prob_tol or using different initial intensity estimates.\n", 
                           "Ignore message if it disappears after a few iterations."))
            nan_warning_issued <- TRUE
          } else if(denom < 0 & !nan_warning_issued){
            warning(paste0("Transition ", ai, " to ", bi, " between times ", li, " and ", ri, " for subject ", i, " in iteration ", it, " impossible with current estimates.
                        A calculated probability is smaller than 0. Try increasing prob_tol or using different initial intensity estimates.\n",
                           "Ignore message if it disappears after a few iterations."))
            nan_warning_issued <- TRUE
          }          
          log_denom <- suppressWarnings(log(denom))
          ll <- ll + log_denom #Add this to likelihood
          #ptf = P_{ai, h}(l_i, t) - h = columns, t = rows
          #(P_{a_i, 1}(l_i, t_1)      ....        P_{a_i, H}(l_i, t_1))
          #(      ....                                                )
          #(      ....                                                )
          #(P_{a_i, 1}(l_i, t_n)      ....        P_{a_i, H}(l_i, t_n))
          
          #**** Backward Probabilities ---------------------------------------------------
          
          ptb <- probtrans_C(int_mats, predt=ri, cutoff = li, direction="fixedhorizon", as.df=FALSE) #Calculate P(t, ri) - backward multiplication
          #Cut-off at li, because we don't need probabilities before that time
          #Extract only probabilities to the state bi we arrive in at ri 
          ptbtemp <- matrix(NA, nrow = dim(ptb)[1], ncol = H+1)
          ptbtemp[,1] <- ptb[, 1, 1]
          for(l in 2:(H+1)) {
            ptbtemp[,l] <- ptb[, bi+1, l-1]
          }
          ptb <- ptbtemp
          #ptb <- as.data.frame(ptbtemp)
          #colnames(ptb) <- c("time", paste0("pstate", 1:H))
          #ptb = P_{g, bi}(t, r_i) - g = columns, t = rows
          #(P_{1, b_i}(t_1, r_i)      ....        P_{H, b_i}(t_1, r_i))
          #(      ....                                                )
          #(      ....                                                )
          #(P_{1, b_i}(t_n, r_i)      ....        P_{H, b_i}(t_n, r_i))
          #browser()
          tmp <- ptf[,-1] * ptb[,-1] / denom #Probabilities of transition via states/overall probability
          #In this way, each component of multiplication allows us to calculate Y_{gi}^{k}
          #Remove last row because we consider only taus \in (l_i, r_i],
          #and last row gives us entry if k = k_{r_i}+1
          
          #**** Update Y and d ----------------------------------------------------------
          
          for (h in 1:H) Y[i, wh, h] <- tmp[-nrow(tmp), h] 
          for (m in 1:M) {
            g <- tmat2$from[m]; h <- tmat2$to[m]
            agh <- int_mats$intensity_matrices[g, h, wh]
            #we used to have ptb[-nrow(ptb), h+1] but that was wrong, should be -1 instead.
            d[i, wh, m] <- ptf[-nrow(ptf),g+1]*agh*ptb[-1,h+1]/denom
          }
          #End of non-exactly observed case (updating of d and Y)
        } else { 
          
          #*** Case: Exactly Observed Transitions --------------------------------------
          
          
          #**** Forward Probabilities ---------------------------------------------------
          ptf <- probtrans_C(int_mats, predt=li, cutoff = ri, direction="forward", as.df=FALSE) #Calculate P(li, t) - forward multiplication
          # ptf <- tryCatch(probtrans(A, predt=li, direction="forward", variance=FALSE),
          #                 warning = function(cond){
          #                   message(conditionMessage(cond))
          #                   print(tail(probtrans(A, predt=li, direction="forward", variance=FALSE)[[1]]))
          #                   probtrans(A, predt=li, direction="forward", variance=FALSE)
          #                 })
          ptf <- ptf[, , ai] #Extract only the probabilities from the state we are in
          #ptf is now P_{a_i *}(l_i, t) with t variable (rows) and * the states we can transition to (columns)
          denom <- ptf[nrow(ptf), bi+1] #Extract P_{ai, bi}(l_i, r_i)
          if(denom == 0 & !nan_warning_issued){
            warning(paste0("Transition ", ai, " to ", bi, " between times ", li, " and ", ri, " for subject ", i, " in iteration ", it, " is impossible with current estimates.
                        An estimated probability was 0, whereas it shouldn't have been. Try increasing prob_tol or using different initial intensity estimates.\n", 
                           "Ignore message if it disappears after a few iterations."))
            nan_warning_issued <- TRUE
          } else if(denom < 0 & !nan_warning_issued){
            warning(paste0("Transition ", ai, " to ", bi, " between times ", li, " and ", ri, " for subject ", i, " in iteration ", it, " impossible with current estimates.
                        A calculated probability is smaller than 0. Try increasing prob_tol or using different initial intensity estimates.\n",
                           "Ignore message if it disappears after a few iterations"))
            nan_warning_issued <- TRUE
          }          
          #ptf = P_{ai, h}(l_i, t) - h = columns, t = rows
          #(P_{a_i, 1}(l_i, t_1)      ....        P_{a_i, H}(l_i, t_1))
          #(      ....                                                )
          #(      ....                                                )
          #(P_{a_i, 1}(l_i, t_n)      ....        P_{a_i, H}(l_i, t_n))
          
          #**** Extra variables for Exactly observed ------------------------------------
          
          #Match bi with entry in exact
          exact_idx <- which(exact == bi)
          #From which states is bi reachable?
          bi_reachable_from <- reachable_from[[exact_idx]]
          #From how many states can bi be reached?
          n_reachable_from <- length(bi_reachable_from)
          #Associated transno for these transitions?
          bi_reachable_from_transno <- reachable_from_transno[[exact_idx]]
          #Largest tau just before ri (we exploit the fact that C++ returns index - 1)
          ri_minus <- taus[binary_search_larger_equal(taus, ri)]
          #We need different transition probabilities in the exact case.
          
          #**** Backward Probabilities --------------------------------------------------
          
          #We no longer calculate probabilities until r_i, but instead until \tau_{r_i-1}
          ptb <- probtrans_C(int_mats, predt=ri_minus, cutoff = li, direction="fixedhorizon", as.df=FALSE) #Calculate P(0, tau_{k_{r_i}-1}) - backward multiplication
          #Extract only probabilities to states from which we can reach bi
          #Store this in 3D array: rows = time, columns = all states g in MSM
          #and 3D dimension = states m from which bi can be reached
          ptbtemp <- array(NA, dim = c(dim(ptb)[1], H+1, n_reachable_from))
          ptbtemp[, 1, ] <- ptb[, 1, 1]
          #Fill ptbtemp matrix
          for(r in seq_along(bi_reachable_from)) { #Iterate over states that bi can be reached from
            for(l in 2:(H+1)) { #Iterate from 2 because first column is time.
              ptbtemp[, l, r] <- ptb[, bi_reachable_from[r]+1 , l-1]  
            }
          }
          dimnames(ptbtemp)[[2]] <- c("time", paste0("pstate", 1:H))
          dimnames(ptbtemp)[[3]] <- bi_reachable_from
          #ptbtemp = P_{g, r}(t, r_i-) - g = columns, t = rows, r (m in overleaf) = 3rd dimension, r_i- fixed.
          #(P_{1, r}(t_1, r_i-)      ....        P_{H, r}(t_1, r_i-))
          #(      ....                                            )
          #(      ....                                            )
          #(P_{1, r}(t_n, r_i-)      ....        P_{H, r}(t_n, r_i-))
          #3rd dimension iterates over r from reachable_from
          
          #We only keep entries for time >= li
          ptbtemp <- ptbtemp[ptb[, 1, 1] >= li, , , drop = FALSE]
          
          
          #**** Calculate sum over reachable from exact state --------
          
          #------------------------------------------------------------------------#
          #Calculation of \sum_{m \in R_bi} \alpha_mbi^k P_gm(\tau_k-1, \tau_k_ri-1)
          #------------------------------------------------------------------------#
          #Create alpha_{mb_i^k}^{t} matrix to calculate sum over r in reachable_from
          arb <- vector(mode = "numeric", length = n_reachable_from)
          for(r in seq_along(bi_reachable_from_transno)) {
            #Get only cumulative hazard estimates for this transition
            Arb <- A$Haz[A$Haz$trans==bi_reachable_from_transno[r],]
            #Find index where A$Haz$time == r_i
            idx <- binary_search_larger_equal(Arb$time, ri) + 1
            #hazard = difference in cumulative hazard
            arb[r] <- Arb$Haz[idx] - Arb$Haz[idx-1]
            
            #Attempt at speed-up: didn't speed up somehow.
            #from <- tmat2$from[bi_reachable_from_transno[r]]
            #to <- tmat2$to[bi_reachable_from_transno[r]]
            #idx <- binary_search_larger_equal(taus, ri) + 1
            #arb[r] <- int_mats$intensity_matrices[from, to , idx]
          }
          
          #We obtain the required summation by matrix * vector multiplication
          ptb_exact <- matrix(NA, nrow = dim(ptbtemp)[1], ncol = dim(ptbtemp)[2]-1)
          for(s in 2:(dim(ptbtemp)[2])) {
            ptb_exact[,s-1] <- rowSums(t(t(as.matrix(ptbtemp[, s, ])) * arb))
          }
          #Below line can be commented away, programming help.
          rownames(ptb_exact) <- ptbtemp[, 1, 1]
          #------------------------------------------------------------------------#
          #ptb_exact = \sum_{m \in R} \alpha_{mb_i}^k P_{g, m}(t, \tau_{k_r_i-1}) - g = columns, t = rows
          #(\sum_{m \in R} alpha_{mbi}^{t_1} * P_{1m}(t_1, \tau_{k_{r_i-1})      ....        \sum_{m \in R} alpha_{mbi}^{t_1} * P_{Hm}(t_1, \tau_{k_{r_i-1}))
          #(      ....                                            )
          #(      ....                                            )
          #(\sum_{m \in R} alpha_{mbi}^{t_n} * P_{1m}(t_1, \tau_{k_{r_i-1})      ....        \sum_{m \in R} alpha_{mbi}^{t_n} * P_{Hm}(t_1, \tau_{k_{r_i-1}))
          #We obtain a matrix of dimensions (relevant times, #states)
          #------------------------------------------------------------------------#
          
          
          #The denominator in the exact case is also different:
          #\sum_{m \in R_bi} \alpha_mbi^k P_aim(l_i, \tau_k_ri-1)
          denom_exact <- ptb_exact[1,ai]
          log_denom_exact <- suppressWarnings(log(denom_exact))
          if(is.nan(log_denom_exact) & !nan_warning_issued){
            warning(paste0("Negative probabilities calculated during iteration ", it, ". \n", 
                           "This is probably due to the initial estimates chosen in 'inits'. \n",
                           "If warning doesn't disappear after a few iterations, try to change the 'inits' argument."))
            nan_warning_issued <- TRUE
          }
          ll <- ll + log_denom_exact
          
          
          #**** Update Y and d ----------------------------------------------------------
          
          
          #--------------------------Calculate Y-----------------------------------#
          #Remove last row from ptf (because \tau_{k-1} cannot be larger than \tau_{k_r_i-1})
          tmp <- ptf[-nrow(ptf),-1] * ptb_exact / denom_exact
          for (h in 1:H) Y[i, wh, h] <- tmp[, h] 
          
          #--------------------------Calculate d-----------------------------------#
          tmp_insert <- rep(0, ncol(ptb_exact))
          tmp_insert[bi] <- 1
          ptb_exact_temp <- rbind(ptb_exact, tmp_insert)
          for (m in 1:M) {
            g <- tmat2$from[m]; h <- tmat2$to[m]
            agh <- int_mats$intensity_matrices[g, h, wh]
            #we used to have ptb[-nrow(ptb), h+1] but that was wrong, should be -1 instead.
            d[i, wh, m] <- ptf[-nrow(ptf),g+1]*agh*ptb_exact_temp[-1,h]/denom_exact
          }
          #End of exactly observed case (updating of d and Y)
        } 
      } #End of for loop over j (times of observation of subject i)
      if(newmet == TRUE){
        
        
        #** Contribution after last observed time (newmet) --------------------------
        
        #UNRESOLVED QUESTION:
        #SHOULD THE CONTRIBUTIONS AFTER LARGEST OBSERVED TIME ALSO CONTIBUTE TO LIKELIHOOD VALUE?
        #EM AIMS TO MAXIMIZE THE MARGINAL LIKELIHOOD.
        #MARGINAL LIKELIHOOD IN OUR CASE IS LIKELIHOOD WHERE WE MISS EXACT TRANSITION TIMES
        
        #WE HAVE CONTRIBUTIONS AFTER LARGEST OBSERVED TIME DUE TO *CENSORING*
        #CENSORED OBSERVATIONS ALWAYS CONTRIBUTE TO LIKELIHOOD IN A DIFFERENT WAY
        #SO PROBABLY YES! WE SHOULD INCORPORATE INTO LIKELIHOOD THIS VALUE AS WELL
        #THIS IS SIMILAR TO CENSORED CONTRIBUTIONS IN NORMAL SURVIVAL.
        #NOTE THAT IF WE LAST SEE A SUBJECT IN AN ABSORBING STATE, THE CONTRIBUTION IS 0!
        
        #IF SO, WE MUST ADD A LINE:
        #ll <- ll + log(denom) similar to line 148 above.
        #UNRESOLVED QUESTION
        
        
        #Above we have calculated Y and d for taus in [t_{i0}, t_{in_i}], now we need to calculate
        #for taus in [t_{in_i}, max_i t_{in_i}]
        wh <- which(taus > max(gdi$time)) #1{\tau \in [t_{in_i}, max_i t_{in_i}]}
        if(!(length(wh) == 0)){
          j <- nrow(gdi)
          li <- gdi$time[j]    #Left bound of current interval
          ai <- gdi$state[j]   #State at left bound of interval
          if(!(length(wh) == 0)){ #If there exist such bins, we must update them as well
            ptf <- probtrans_D(int_mats, predt=li, direction="forward", as.df=FALSE) #Calculate P(li, t) - forward multiplication
            ptf <- ptf[, , ai] #Extract only the probabilities from the state we are in
          }
          for (h in 1:H) {
            Y[i, wh, h] <- ptf[-nrow(ptf), h+1]
          }
          for (m in 1:M) {
            g <- tmat2$from[m]; h <- tmat2$to[m]
            agh <- int_mats$intensity_matrices[g, h, wh]
            d[i, wh, m] <- ptf[-nrow(ptf),g+1]*agh
          }
        }
      } #End of if for newmet == TRUE (contributions after last observed time)
      
    } #End of loop over i (subjects)
    
    
    #* KKT conditions (current iteration) --------------------------------------
    
    #We need to check whether we would violate the KKT conditions if we would use 
    #The usual update rule alpha = d/Y.
    #If we would violate the conditions, we need to re-adjust the hazards in this step
    
    #Initiate mu_g matrix with rows representing states and columns representing times (taus)
    mu_g <- matrix(NA, nrow = H, ncol = K)
    #We need to have mu_g^k = sum_{h \leftarrow g} d_{gh}^k - Y_g^k
    for(g in 1:H){ #For each state g
      trans_t <- which(tmat2$from == g) #Which direct transitions start in g?
      #We calculate sum_{h \leftarrow g} d_{gh}^k resulting in a vector of length K
      sum_d <- apply(d[, ,trans_t], 2, sum)
      #Also calculate Y_g^k resulting in a vector of length K
      Y_g <- apply(Y[, , g], 2, sum)
      KKT_check <- any(sum_d > Y_g)
      Ylol <<- Y
      dlol <<- d
      if(is.na(KKT_check)){
        stop("An observed transition has become impossible with the current estimates. 
             This is likely due to a probability being set to 0. Try lowering 'prob_tol' to resolve the issue.")
      }
      if(KKT_check){
        KKT_violated <- KKT_violated + 1
        if(verbose){
          message(paste0("KKT condition (mu_g) violated for the intensities out of state ", g, " at time(s) ", taus[which(sum_d > Y_g)], " in iteration ", it))
          message("Estimates were adjusted to lie within the optimisation region.")  
        }
      }
      mu_g[g, ] <- pmax(sum_d - Y_g, 0) #When smaller than 0, we are not violating the constraint in next iteration.
    }
    
    if(final_step){
      break
    }
    
    
    #* Initiate + Assign new estimates -----------------------------------------
    
    #--------------------------------Initiate new estimates--------------------------#
    newA <- list(Haz = data.frame(time=rep(taus, M), Haz=rep(0,M),
                                  trans=rep(1:M, rep(K,M))),
                 trans = tmat)
    attr(newA, "class") <- "msfit"
    
    #-----------------------------Check whether anything can be set to zero----------#
    small_probs_d <- which(d < prob_tol)
    small_probs_Y <- which(Y < prob_tol)
    
    if(length(small_probs_d) != 0){
      d[small_probs_d] <- 0
    }
    if(length(small_probs_Y) != 0){
      Y[small_probs_Y] <- 0
    }
    
    
    
    #-------------------------------Assign new estimates-----------------------------#
    for (m in 1:M) {
      g <- tmat2$from[m]
      agh <- apply(d[,,m], 2, sum) / (apply(Y[,,g], 2, sum) + mu_g[g, ]) #mu_g = 0 when no violation
      agh[is.na(agh)] <- 0
      newA$Haz$Haz[newA$Haz$trans==m] <- cumsum(agh)
    }
    
    #browser()
    
    
    #--------------Update estimates and parameters for next iteration----------------#
    
    #Determine whether stopping criterion has been reached
    #stop_crit becomes TRUE when algorithm has converged.
    #stop_crit <- switch(conv_crit,
    #                    "haz" = all(abs(A$Haz$Haz - newA$Haz$Haz) < tol),
    #                    "prob" = all(abs(d_old - d) < tol) & all(abs(Y_old - Y) < tol),
    #                    "lik" = (ll - llold) < tol)  
    
    #Removing bins with bass 0 - IN PROGRESS
    #If we have removed bins, we cannot determine whether convergence has been reached
    if(remove_bins && length(remove > 0)){
      stop_crit <- FALSE
    } else{
      stop_crit <- switch(conv_crit,
                          "haz" = all(abs(A$Haz$Haz - newA$Haz$Haz) < tol),
                          "prob" = all(abs(d_old - d) < tol) & all(abs(Y_old - Y) < tol),
                          "lik" = (ll - llold) < tol)  
    }
    
    
    d_old <- d
    Y_old <- Y
    A <- newA
    #Change in likelihood value
    delta <- ll - llold
    ll_storage <- c(ll_storage, ll)
    llold <- ll
    if(verbose){
      message(paste0(it, ll))  
    }
    
    
    #Stopping criteria check
    if(stop_crit){
      if(checkMLE){
        final_step <- TRUE
      } else{
        break
      }
    }
    nan_warning_issued <- FALSE
    #browser()
  } #End of iteration over it (EM algorithm iterations)
  
  
  
  # Check Convergence of EM -------------------------------------------------
  
  if(checkMLE){
    #Reduced gradient is left side of KKT condition (1)
    reduced_gradient <- array(NA, dim = c(M, K)) #Store the reduced gradient in an array
    #We need Y_g^k for all g and k
    Ygk <- apply(Y, 2, colSums)
    #Ygk is now a (H, K) matrix with states h in rows and times k in columns
    #We need d_{gh}^k for all m(g->h transitions) and k
    Dmk <- apply(d, 2, colSums)
    #Dmk is now a (M, K) matrix with transitions g->h(m) in rows and times k in columns
    #We need current estimates of alpha
    aghmat <- matrix(NA, nrow = M, ncol = K)
    #Small fix for when we only have 1 transition
    if(M == 1){
      Dmk <- array(Dmk, dim = c(1, K))
    }
    for(m in 1:M){
      aghmat[m,] <- diff(c(0, A$Haz$Haz[A$Haz$trans==m]))
      g <- tmat2$from[m]; h <- tmat2$to[m]
      reduced_gradient[m, ] <- Ygk[g, ] - Dmk[m, ]/aghmat[m, ] + mu_g[g, ]
    }
    reduced_gradient[which(is.nan(reduced_gradient))] <- 0
    
    
    #----------Determine whether we have found the MLE-------------------------#
    if(any(reduced_gradient > checkMLE_tol)){ 
      #If any reduced_gradient is non-zero or any lagrange multiplier is negative
      #We have not arrived at the MLE (yet)
      #We need to check if the gradient is non-zero due to a numerical issue
      which_gradient_nonzero <- which(reduced_gradient > checkMLE_tol)
      if( any(aghmat[which_gradient_nonzero] > checkMLE_tol) ){
        isMLE <- FALSE
      } else{ #Numerical reason for non-MLE, the current estimate hasn't been set to zero yet (but is close)
        isMLE <- TRUE
      }
    } else { #Else we are at the MLE! :party:
      isMLE <- TRUE
    }
    rownames(reduced_gradient) <- 1:M
    colnames(reduced_gradient) <- round(taus, 2)
    names(attributes(reduced_gradient)$dimnames) <- c("transitions", "times")
    rownames(mu_g) <- 1:H
    colnames(mu_g) <- round(taus, 2)
    names(attributes(mu_g)$dimnames) <- c("states", "times")
  }
  
  
  
  # Output ------------------------------------------------------------------
  out <- list(A = A,
              Ainit = Ainit,
              gd = gd,
              ll = ll,
              delta = delta,
              it = it,
              taus = taus,
              tmat = tmat, 
              tmat2 = tmat2,
              ll_history = ll_storage,
              KKT_violated = KKT_violated,
              min_time = min_time)
  if(checkMLE){
    out$reduced_gradient <- reduced_gradient
    out$isMLE <- isMLE
    out$lagrangemultiplier <- mu_g
    out$aghmat <- aghmat
    out$Ygk <- Ygk
    out$Dmk <- Dmk
  }
  return(out)
}




