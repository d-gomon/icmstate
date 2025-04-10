#' @title
#' Calculate subject specific transition probabilities from a 
#' multi-state \code{\link[survival:coxph]{coxph}} model.
#' 
#' @description
#' Given a coxph model fit on multi-state data (prepared with \code{\link[mstate:msprep]{msprep}}),
#' determine transition probabilities for subjects in \code{newdata}.
#' 
#' @param object A \code{\link[survival:coxph]{coxph}} object fit on multi-state
#' data. Must contain a \code{strata(X)} term. Data used for the coxph() fit preferably prepared 
#' using \code{\link[mstate:msprep]{msprep}}.
#' @param predt A positive number indicating the prediction time. This is either 
#' the time at which the prediction is made (if \code{direction = "forward"}) or 
#' the time for which the prediction is to be made (if \code{direction = "backward"}).
#' @param direction One of \code{"forward"} (default) or \code{"fixedhorizon"},
#'  indicating whether prediction is forward or for a fixed horizon
#' @param newdata A \code{data.frame} containing a single row for each transition per
#' subject in the data. For a model with m possible transitions, and n subjects 
#' \code{newdata} must have m*n rows. It must contain the following named columns:
#' \describe{
#'   \item{\code{id}:}{Unique identifier of the subject, can be numeric or character;}
#'   \item{\code{from}:}{State from which the transition takes place;}
#'   \item{\code{to}:}{State to which the transition takes place;}
#'   \item{\code{trans}:}{Transition number in the \code{'transMat'} \code{trans} this transition relates to;}
#'   \item{\code{"variables"}:}{The variables and their values for the subject
#'   identified by "id" for the transition this entry relates to. Names must 
#'   match the names of the variables in coxph \code{object}.}
#' } Note that newdata must contain a column containing the variable which was 
#' used to determine the stratum of a transition in \code{object}. 
#' Usually the stratum is determined from one of the required columns. The 
#' "variables" columns can usually be obtained using \code{\link[mstate:expand.covs]{expand.covs}}.
#' @param trans A transition matrix as created by \code{\link[mstate:transMat]{transMat}}.
#' 
#' 
#' 
#' @return An object of class \code{"probtrans.subjects"}.
#' This is a list of length n (number of subjects in newdata), with each list element
#' an object of class \code{\link[mstate:probtrans]{probtrans}} for the associated 
#' subject. List elements can be accessed using \code{[[x]]}, with \code{x} 
#' ranging from 1 to n. Additionally, each list element 
#' has an element \code{$id}, representing the subject id and the output object 
#' also has an element \code{$subject_ids} representing the subject ids in order.
#' 
#' @details
#' When using this function for \code{newdata} with many subjects, consider running the 
#' function multiple times for parts of \code{newdata} to negate the risk of 
#' running our of memory.
#' 
#' 
#' @importFrom stats terms
#' 
#' @export
#' 
#' 
#' 
#' @examples
#' #Example from the mstate vignette
#' #We determine the subject specific transition probabilities for subjects
#' #in the ebmt3 data-set 
#' if(require("mstate")){
#'   data(ebmt3)
#'   n <- nrow(ebmt3)
#'   tmat <- transMat(x = list(c(2, 3), c(3), c()), names = c("Tx",
#'                                                            "PR", "RelDeath"))
#'   ebmt3$prtime <- ebmt3$prtime/365.25
#'   ebmt3$rfstime <- ebmt3$rfstime/365.25
#'   covs <- c("dissub", "age", "drmatch", "tcd", "prtime")
#'   msbmt <- msprep(time = c(NA, "prtime", "rfstime"), status = c(NA,
#'                   "prstat", "rfsstat"), data = ebmt3, trans = tmat, keep = covs)
#'   #Expand covariates so that we can have transition specific covariates
#'   msbmt <- expand.covs(msbmt, covs, append = TRUE, longnames = FALSE)
#'   
#'   #Create extra variable to allow gender mismatch to have the same effect 
#'   #for transitions 2 and 3.
#'   msbmt$drmatch.2.3 <- msbmt$drmatch.2 + msbmt$drmatch.3
#'   
#'   #Introduce pr covariate for proportionality assumption of transitions 2 and 3
#'   #(only used in models 2 and 4)
#'   msbmt$pr <- 0
#'   msbmt$pr[msbmt$trans == 3] <- 1
#'   
#'   
#'   
#'   #-------------Models---------------------#
#'   
#'   #Simple model, transition specific covariates, each transition own baseline hazard
#'   c1 <- coxph(Surv(Tstart, Tstop, status) ~ dissub1.1 + dissub2.1 +
#'                 age1.1 + age2.1 + drmatch.1 + tcd.1 + dissub1.2 + dissub2.2 +
#'                 age1.2 + age2.2 + drmatch.2 + tcd.2 + dissub1.3 + dissub2.3 +
#'                 age1.3 + age2.3 + drmatch.3 + tcd.3 + strata(trans), data = msbmt,
#'                 method = "breslow")
#'   
#'   #Model with same baseline hazards for transitions 2 (1->3) and 3(2->3)
#'   #pr then gives the ratio of the 2 hazards for these transitions
#'   c2 <- coxph(Surv(Tstart, Tstop, status) ~ dissub1.1 + dissub2.1 +
#'                 age1.1 + age2.1 + drmatch.1 + tcd.1 + dissub1.2 + dissub2.2 +
#'                 age1.2 + age2.2 + drmatch.2 + tcd.2 + dissub1.3 + dissub2.3 +
#'                 age1.3 + age2.3 + drmatch.3 + tcd.3 + pr + strata(to), data = msbmt,
#'                 method = "breslow")
#'   
#'   #Same as c2, but now Gender mismatch has the same effect for both 
#'   c4 <- coxph(Surv(Tstart, Tstop, status) ~ dissub1.1 + dissub2.1 +
#'                 age1.1 + age2.1 + drmatch.1 + tcd.1 + dissub1.2 + dissub2.2 +
#'                 age1.2 + age2.2 + drmatch.2.3 + tcd.2 + dissub1.3 + dissub2.3 +
#'                 age1.3 + age2.3  + tcd.3 + pr + strata(to), data = msbmt,
#'               method = "breslow")
#'               
#'   
#'   #We need to make a data.frame containing all subjects of interest
#'   ttmat <- to.trans2(tmat)[, c(2, 3, 1)]
#'   names(ttmat)[3] <- "trans"
#'   nd_n <- NULL
#'   for (j in 1:30) {
#'     # Select global covariates of subject j
#'     cllj <- ebmt3[j, covs]
#'     nd2 <- cbind(ttmat, rep(j, 3), rbind(cllj, cllj, cllj))
#'     colnames(nd2)[4] <- "id"
#'     # Make nd2 of class msdata to use expand.covs
#'     attr(nd2, "trans") <- tmat
#'     class(nd2) <- c("msdata", "data.frame")
#'     nd2 <- expand.covs(nd2, covs=covs, longnames = FALSE)
#'     nd2$drmatch.2.3 <- nd2$drmatch.2 + nd2$drmatch.3
#'     nd2$pr <- 0
#'     nd2$pr[nd2$trans==3] <- 1
#'     nd2$strata <- c(1, 2, 2)
#'     nd_n <- rbind(nd_n, nd2)
#'    }
#'    
#'    icmstate_pt <- probtrans_coxph(c2, predt = 0, direction = "forward", 
#'                                   newdata = nd_n, trans = tmat)
#'                                   
#'    #Now we can plot the transition probabilities for each subject separately:
#'    plot(icmstate_pt[[1]])
#'    #icmstate_pt has length number of subjects in newdata
#'    #And icmstate_pt[[i]] is an object of class "probtrans", so you can 
#'    #use all probtrans functions: summary, plot etc.
#' }
#' 
#' 



probtrans_coxph <- function(object, predt, 
                                       direction = c("forward", "fixedhorizon"),
                                       newdata, trans){
  #Given:
  #object: coxph model fitted for MSM data
  #newdata: data.frame containing the variables for new subjects to make predictions
  #on, must contain all variables in object
  #trans: a transition matrix 
  
  
  #---------------DATA CHECKS----------------#
  
  #We first copy some checks from msfit, which are actually from survfit
  #Here we make sure that the cox fit is appropriate, and extract the 
  #correct terms for each transition.
  
  #Some checks on coxmod
  if(!is.null((object$call)$weights) || !is.null(object$weights))
    stop("msfit cannot (yet) compute the result for a weighted model")
  Terms <- terms(object)
  strat <- attr(Terms, "specials")$strata
  if (is.null(strat)) stop("object has to have strata() term")
  cluster <- attr(Terms, "specials")$cluster
  if (length(cluster)) stop("cluster terms are not supported")
  if (!is.null(attr(object$terms, "specials")$tt))
    stop("msfit cannot yet process coxph models with a tt term")
  
  
  direction <- match.arg(direction, choices = c("forward", "fixedhorizon"))
  
  if(!(predt >= 0) | !(is.numeric(predt))){
    stop("Prediction time must be a positive numeric value.")
  }
  
  if(!inherits(trans, "matrix")){
    stop("trans must be a transition matrix (array).")
  }
  
  if(!is.data.frame(newdata)){
    stop("newdata must be a data frame containing columns id, from, to transno, strata")
  }
  
  
  #Keep track of subject ids.
  if("id" %in% colnames(newdata)){
    subject_ids <- unique(newdata[, "id"])
    n_subjects <- length(subject_ids)
  } else{
    subject_ids <- 1
    n_subjects <- 1
  }
  
  
  #----------Step 1: Baseline intensities-------------#
  #Recover baseline intensities from cox model
  #Output from get_intensity matrices
  #List with $intensity_matrices and $unique_times
  #print("Step 1")
  baseline_intensities <- baseline_intensities_from_coxmod(object, trans)

  
  #---------Step 2: Subject specific risks-------------#
  
  #####Step 2.1: Transform data we have into data which has all 
  #the necessary variables######
  
  #Next we want to transform newdata into a format which contains all 
  #the variables we have in the coxph object.
  
  #newdata <- expand_covariates_long_data(newdata)
  
  #Decided not to do this step, as this restricts the possible models which 
  #can be fit considerably.
  
  #####Step 2.2: Extract the subject specific risk from the transformed data#####
  #print("Step 2")
  subject_specific_risks <- trans_specific_risks(object = object, newdata = newdata,
                                                 trans = trans)
  
  
  #---------Step 3: Obtain subject specific predictions----------#
  
  #To obtain subject specific predictions we simply multiply the risk matrix
  #with the appropriate indices of the baseline intensity matrices
  #print("Step 3")
  subject_specific_intensity_matrices <- subject_specific_intensity_matrices(subject_specific_risks = subject_specific_risks,
                                                                             baseline_intensities = baseline_intensities,
                                                                             trans = trans)
  #The resulting output is a array of dimensions (S, S, n_times, n_subjects) 
  #with S the number of states in the model
  
  
  #--------Step 4: Calculate transition probabilities per subject------#
  
  #This step is relatively easy. We simply use probtrans for each subject separately.
  #print("Step 4")
  
  
  #We now determine n_subjects and subject_ids in the beginning of this function.
  #n_subjects <- dim(subject_specific_intensity_matrices$subject_intensity_matrices)[4]
  #Note that subject_ids will now have class "character", as it was transformed to a name before
  #subject_ids <- dimnames(subject_specific_intensity_matrices$subject_intensity_matrices)[[4]]
  res <- vector(mode = "list", length = n_subjects)
  for(i in 1:n_subjects){
    res[[i]] <- probtrans_D(list(intensity_matrices = subject_specific_intensity_matrices$subject_intensity_matrices[, , , i],
                                 unique_times = subject_specific_intensity_matrices$unique_times), 
                            predt = predt, direction = direction, as.df = TRUE)
    class(res[[i]]) <- "probtrans"
    res[[i]]$trans <- trans
    res[[i]]$method <- "aalen"
    res[[i]]$predt <- predt
    res[[i]]$direction <- direction
    res[[i]]$id <- subject_ids[i]
  }
  res$subject_ids <- subject_ids
  class(res) <- "probtrans.subjects"
  #-----------OUTPUT------------#
  return(res)
}


#' Second (hopefully faster) version of probtrans coxph
#' 
#' @description
#' Currently not faster than probtrans_coxph()
#' However, this function no longer stores the intensity matrices
#' for each subject in a 3D array, instead overwriting a single 2D matrix
#' with the values for a different subject iteratively. Unfortunately, due to
#' how R manages memory this does not lead to a speed/memory improvement.
#' 
#' 
#' @keywords internal
#' @noRd
#' 
#' 
#' 

probtrans_coxph2 <- function(object, predt, 
                            direction = c("forward", "fixedhorizon"),
                            newdata, trans){
  #Given:
  #object: coxph model fitted for MSM data
  #newdata: data.frame containing the variables for new subjects to make predictions
  #on, must contain all variables in object
  #trans: a transition matrix 
  
  
  #---------------DATA CHECKS----------------#
  
  #We first copy some checks from msfit, which are actually from survfit
  #Here we make sure that the cox fit is appropriate, and extract the 
  #correct terms for each transition.
  
  #Some checks on coxmod
  if(!is.null((object$call)$weights) || !is.null(object$weights))
    stop("msfit cannot (yet) compute the result for a weighted model")
  Terms <- terms(object)
  strat <- attr(Terms, "specials")$strata
  if (is.null(strat)) stop("object has to have strata() term")
  cluster <- attr(Terms, "specials")$cluster
  if (length(cluster)) stop("cluster terms are not supported")
  if (!is.null(attr(object$terms, "specials")$tt))
    stop("msfit cannot yet process coxph models with a tt term")
  
  
  direction <- match.arg(direction, choices = c("forward", "fixedhorizon"))
  
  if(!(predt >= 0) | !(is.numeric(predt))){
    stop("Prediction time must be a positive numeric value.")
  }
  
  if(!inherits(trans, "matrix")){
    stop("trans must be a transition matrix (array).")
  }
  
  if(!is.data.frame(newdata)){
    stop("newdata must be a data frame containing columns id, from, to transno, strata")
  }
  
  
  #----------Step 1: Baseline intensities-------------#
  #Recover baseline intensities from cox model
  #Output from get_intensity matrices
  #List with $intensity_matrices and $unique_times
  #print("Step 1")
  baseline_intensities <- baseline_intensities_from_coxmod(object, trans)
  
  
  #---------Step 2: Subject specific risks-------------#
  
  #####Step 2.1: Transform data we have into data which has all 
  #the necessary variables######
  
  #Next we want to transform newdata into a format which contains all 
  #the variables we have in the coxph object.
  
  #newdata <- expand_covariates_long_data(newdata)
  
  #TODOTODOTODO - for now newdata simply contains all necessary variables, 
  #and is in long format, with 1 line per possible transition for each subject!
  #If we don't have a line for each transition, we may miss some covariates
  #which were included for proportionality (see rel.3)
  
  #####Step 2.2: Extract the subject specific risk from the transformed data#####
  #print("Step 2")
  subject_specific_risks <- trans_specific_risks(object = object, newdata = newdata,
                                                 trans = trans)
  
  
  
  #--------Steps 3 & 4: Calculate intensity matrices and transition probabilities per subject------#
  
  #For each subject, determine their intensity matrices for all unique times.
  #Then simply use probtrans for each subject separately.
  
  n_subjects <- nrow(subject_specific_risks)
  n_transitions <- ncol(subject_specific_risks)
  n_times <- length(baseline_intensities$unique_times)
  
  #Calculate a baseline intensity matrix with zero diagonals
  #We can use this one to calculate the off-diagonal elements 
  #and then calculate the diagonal elements for that specific subject
  intensity_matrices_zero_diag <- baseline_intensities$intensity_matrices
  for(k in 1:n_times){
    diag(intensity_matrices_zero_diag[, , k]) <- 0
  }
  
  #pre-specify output (list of length n_subjects)
  res <- vector(mode = "list", length = n_subjects)
  
  for(i in 1:n_subjects){
    #subject_specific_intensity_matrix returns a 3D array containing the
    #intensity matrices for the current subject.
    subject_specific_intensity_matrix <- subject_specific_intensity_matrix(subject_specific_risk = subject_specific_risks[i,],
                                                                           intensity_matrices_zero_diag = intensity_matrices_zero_diag,
                                                                           trans = trans, n_transitions = n_transitions,
                                                                           n_times = n_times)
    res[[i]] <- probtrans_D(list(intensity_matrices = subject_specific_intensity_matrix,
                                 unique_times = baseline_intensities$unique_times), 
                            predt = predt, direction = direction, as.df = TRUE)
    class(res[[i]]) <- "probtrans"
    res[[i]]$trans <- trans
    res[[i]]$method <- "aalen"
    res[[i]]$predt <- predt
    res[[i]]$direction <- direction
  }
  
  
  #-----------OUTPUT------------#
  return(res)
}






#' Recover baseline intensities (in the form of intensity_matrices) from a 
#' coxph() fit on multi-state data.
#' 
#' @inherit probtrans_coxph params
#' @param tmat A transition matrix as created by \code{\link[mstate:transMat]{transMat}}.
#' 
#' @importFrom mstate msfit to.trans2
#' @import survival
#' @importFrom stats formula model.frame as.formula predict
#' 
#' @keywords internal


baseline_intensities_from_coxmod <- function(object, tmat){
  #Extract baseline intensities from a cox model
  #Before running this function, check whether the cox model has been fit 
  #for a MSM (contains strata & no interaction/tt terms)
  
  #object must be a coxph model with strata
  #tmat must be a transition matrix
  #Create a baseline subject, which has 0 variables everywhere!
  #We need strata here, as we want the baseline hazard for each strata separately
  ttmat <- to.trans2(tmat)[, c(2, 3, 1)]
  #Setting the name to "trans" is in line with msprep... This kind of requires people to use msprep.
  names(ttmat)[3] <- "trans"
  baseline_subject <- matrix(0, ncol = length(object$coefficients) + 3, 
                             nrow = nrow(ttmat), 
                             dimnames = list(NULL, c(colnames(ttmat), 
                                                     names(object$coefficients))))
  baseline_subject <- as.data.frame(baseline_subject)
  baseline_subject[1:ncol(ttmat), 1:nrow(ttmat)] <- ttmat
  #Add strata to the subject observations
  #Check how the strata were created
  Terms <- terms(object)
  stangle <- untangle.specials(Terms, 'strata')
  #stangle$vars is now the formula for strata creation
  strata_column <- model.frame(as.formula(paste0(stangle$vars, " ~ 1")), 
                               data = baseline_subject)
  baseline_subject <- cbind(baseline_subject, strata_column)
  
  #Use msfit once to obtain baseline hazards
  msfit_baseline <- mstate::msfit(object, newdata = baseline_subject, trans = tmat)
  
  #Now extract the intensities into matrix form
  baseline_intensities <- get_intensity_matrices(msfit_baseline)
  
  
  
  out <- baseline_intensities
  return(out)
}



#' Expand covariates for a data frame so that covariates can be transition 
#' specific.
#' 
#' @inherit probtrans_coxph params
#' 
#' @keywords internal


expand_covariates_long_data <- function(newdata){
  #Before using trans_specific risks, newdata must be put into the appropriate format.
  #For this, it must be expanded to include all covariates present in object.
  #The question on how to do this is still open!
  
  #Expand the data here to include all necessary covariates 
  #and make sure it's in long format, with a single row for each transition at least!
  #Must also include ttmat, so tmat2 with from, to, trans rows! Otherwise predict.coxph() cannot determine the 
  #strata
  #If we don't make sure there is a single row for each transition, we may
  #miss covariates which make sure we have a proportional hazard for some transitions
}



#' Calculate subject specific risks for subjects in newdata
#' 
#' @description Return a 2 dimensional array, with subjects in the first dimension and
#' transition (numbers) in the second. Can expand later to include time-dependent
#' covariates by introducing extra dimension. Entries of the matrix are exp(lp)_i^m,
#' with i denoting the subject and m the transition.
#' 
#' @inherit probtrans_coxph params
#' 
#' @import survival
#' 
#' @keywords internal
#' 

trans_specific_risks <- function(object, newdata, trans){
  #Given a coxph object (fitted on multi-state data)
  #We want to extract the subject-specific risk (exp(lp)_1, exp(lp)_2, ..., exp(lp)_k)
  #for each transition and each unique subject in the data.
  
  #object = coxph() fit
  #newdata = data.frame/matrix with covariates IN EXTENDED FORM. Must contain id 
  #trans = transition matrix
  #strata_column = a column vector indicating which transitions belong to which strata,
  #must match with the to.trans2(trans) ordering.
  
  #Required output:
  #2D array with:
  #rows: subjects i = 1, ..., n
  #columns: transitions m = 1, ..., M
  #3rd dimension (future): unique times t = t_1, ..., t_K
  #Each entry of the array should be exp(lp)_{i, m, t}
  #so the risk for person i in transition m (at time t (extension, consider later))
  
  
  
  tmat2 <- to.trans2(trans)
  n_transitions <- nrow(tmat2)
  
  #If no "id" column is supplied, assume all entries are from the same subject.
  #We then need to check whether there's exactly n_transitions lines in the data
  if("id" %in% colnames(newdata)){
    subject_ids <- unique(newdata[, "id"])
    n_subjects <- length(subject_ids)
  } else{
    if(nrow(newdata) == n_transitions){
      newdata <- cbind(newdata, data.frame(id = rep(1, nrow(newdata))))
      subject_ids <- 1
      n_subjects <- 1
    } else{
      stop("newdata must have a named column 'id', indicating subjects or newdata 
           must contain at most a single subject (n_transitions rows).")
    }
  }
  
  
  #Initialize output matrix
  subject_specific_risks_mat <- matrix(NA, nrow = n_subjects, ncol = n_transitions)
  rownames(subject_specific_risks_mat) <- subject_ids
  
  #Determine risk for each subject separately
  for(i in seq_along(subject_ids)){
    tempdat <- newdata[newdata[, "id"] == subject_ids[i],]
    #For each transition, we now obtain the subject risk in a vector of length n_transitions
    risk_subject <- predict(object = object, newdata = tempdat, type = "risk",
                            reference = "zero")
    subject_specific_risks_mat[i,] <- risk_subject
  }
  
  return(subject_specific_risks_mat)
}


#' Calculate the subject specific intensity matrices
#' 
#' @description
#' For each subject, calculate a 3D array containing the states (from, to) in the
#' first two dimension and the times in the third.
#' This is then stored in a 4D array, with the 4th dimension indicating each unique subject.
#' 
#' 
#' 
#' 
#' 
#' @keywords internal

subject_specific_intensity_matrices <- function(subject_specific_risks, baseline_intensities, trans){
  #Given subject specific risks, we want to obtain the intensity matrices for each subject
  
  
  tmat2 <- to.trans2(trans)
  n_subjects <- nrow(subject_specific_risks)
  n_transitions <- nrow(tmat2)
  n_times <- length(baseline_intensities$unique_times)
  
  #pre-specify the output
  #First 3 dimensions: from, to, time
  #Fourth dimension: subjects
  intensity_matrices_zero_diag <- baseline_intensities$intensity_matrices
  for(k in 1:n_times){
    diag(intensity_matrices_zero_diag[, , k]) <- 0
  }
  subject_intensity_matrices <- array(intensity_matrices_zero_diag, 
                                      dim = c(dim(intensity_matrices_zero_diag), n_subjects))
  dimnames(subject_intensity_matrices)[[4]] <- rownames(subject_specific_risks)
  
  for(i in 1:n_subjects){
    for(m in 1:n_transitions){
      from <- tmat2$from[m]
      to <- tmat2$to[m]
      subject_intensity_matrices[from, to, , i] <- subject_intensity_matrices[from, to, , i] * subject_specific_risks[i, m]
    }
    for(k in 1:n_times){
      diag(subject_intensity_matrices[, , k, i]) <- 1 - rowSums(subject_intensity_matrices[, , k, i])
    }
  }
  
  out <- list(subject_intensity_matrices = subject_intensity_matrices,
              unique_times = baseline_intensities$unique_times)
  
  return(out)
  
}


#' Determine subject specific intensity matrix for a single subject
#'
#' @keywords internal
#'
#' @importFrom mstate to.trans2
#'
#'


subject_specific_intensity_matrix <- function(subject_specific_risk, 
                                              intensity_matrices_zero_diag, trans,
                                              n_transitions, n_times){
  #Given subject specific risks, we want to obtain the intensity matrices for each subject
  
  
  tmat2 <- to.trans2(trans)
  
  #Initialize output to be the baseline intensities with zero diagonals
  subject_intensity_matrix <- intensity_matrices_zero_diag
  
  
  for(m in 1:n_transitions){
    from <- tmat2$from[m]
    to <- tmat2$to[m]
    subject_intensity_matrix[from, to, ] <- subject_intensity_matrix[from, to, ] * subject_specific_risk[m]
  }
  for(k in 1:n_times){
    diag(subject_intensity_matrix[, , k]) <- 1 - rowSums(subject_intensity_matrix[, , k])
  }
  
  out <- subject_intensity_matrix
  
  return(out)
  
}


