#' @title
#' Plot an object of class \code{"probtrans.subjects"}
#' 
#' @description
#' Plots the transition probabilities for a specific subject. Wrapper for 
#' \code{\link[mstate:plot.probtrans]{plot.probtrans}}
#' 
#' @param x An object of class \code{"probtrans.subjects"}
#' @param id Subject identifier
#' @param ... Further arguments to \code{\link[mstate:plot.probtrans]{plot.probtrans}}
#' 
#' @details
#' Note that 
#'
#' @author Hein Putter and Daniel Gomon
#' 
#' @import mstate
#' 
#' @method plot probtrans.subjects
#' @export
#' 
#' @examples
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
#'   #Simple model, transition specific covariates, each transition own baseline hazard
#'   c1 <- coxph(Surv(Tstart, Tstop, status) ~ dissub1.1 + dissub2.1 +
#'                 age1.1 + age2.1 + drmatch.1 + tcd.1 + dissub1.2 + dissub2.2 +
#'                 age1.2 + age2.2 + drmatch.2 + tcd.2 + dissub1.3 + dissub2.3 +
#'                 age1.3 + age2.3 + drmatch.3 + tcd.3 + strata(trans), data = msbmt,
#'                 method = "breslow")
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
#'     nd_n <- rbind(nd_n, nd2)
#'    }
#'    
#'    icmstate_pt <- probtrans_coxph(c1, predt = 0, direction = "forward", 
#'                                   newdata = nd_n, trans = tmat)
#' 
#'    #plot transition probabilities for subject 2
#'    plot(icmstate_pt, id = 2)
#' }
#' 


plot.probtrans.subjects <- function(x, id, ...){
  if(!inherits(x, "probtrans.subjects")){
    stop("Object to plot must be of class 'probtrans.subjects'.")
  }
  
  if(length(id) > 1){
    stop("Argument 'id' must be at most of length 1.")
  }
  
  which_id <- which(x[["subject_ids"]] == id)  

  if(length(which_id) == 1){
    plot(x[[which_id]], ...)
  } else{
    stop("Subject 'id' not found in 'probtrans.subjects' object. 
         Please provide valid 'id'.")
  }
}



#' Summary method for a probtrans.subjects object
#' 
#' Summary method for an object of class 'probtrans.subjects'. It prints a selection of
#' the estimated transition probabilities. Wrapper for 
#' \code{\link[mstate:summary.probtrans]{summary.probtrans}}.
#' 
#' @aliases summary.probtrans.subjects
#' @param object Object of class 'probtrans.subjects', containing estimated transition
#' probabilities from and to all states in a multi-state model
#' @param id Subject identifier
#' @param times Time points at which to evaluate the transition probabilites
#' @param from Specifies from which state the transition probabilities are to
#' be printed. Should be subset of 1:S, with S the number of states in the
#' multi-state model. Default is print from state 1 only. User can specify
#' from=0 to print transition probabilities from all states
#' @param to Specifies the transition probabilities to which state are to be
#' printed. User can specify to=0 to print transition probabilities to all
#' states. This is also the default
#' @param extend logical value: if \code{TRUE}, prints information for all
#' specified times, even if there are no subjects left at the end of the
#' specified times. This is only valid if the times argument is present
#' @param \dots Further arguments to \code{\link[mstate:summary.probtrans]{summary.probtrans}}
#' 
#' 
#' @import mstate
#' 
#' @return Function \code{summary.probtrans} returns an object of class
#' "summary.probtrans.subjects", which is a list (for each \code{from} state) of
#' transition probabilities at the specified (or all) time points. The
#' \code{print} method of a \code{summary.probtrans.subjects} doesn't return a value.
#' @author Hein Putter and Daniel Gomon
#' @seealso \code{\link[mstate:summary.probtrans]{summary.probtrans}}, 
#' \code{\link{predict_tp}}
#' @keywords print
#' @examples
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
#'   #Simple model, transition specific covariates, each transition own baseline hazard
#'   c1 <- coxph(Surv(Tstart, Tstop, status) ~ dissub1.1 + dissub2.1 +
#'                 age1.1 + age2.1 + drmatch.1 + tcd.1 + dissub1.2 + dissub2.2 +
#'                 age1.2 + age2.2 + drmatch.2 + tcd.2 + dissub1.3 + dissub2.3 +
#'                 age1.3 + age2.3 + drmatch.3 + tcd.3 + strata(trans), data = msbmt,
#'                 method = "breslow")
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
#'     nd_n <- rbind(nd_n, nd2)
#'    }
#'    
#'    icmstate_pt <- probtrans_coxph(c1, predt = 0, direction = "forward", 
#'                                   newdata = nd_n, trans = tmat)
#' 
#'    #Obtain summary of probtrans.subjects object
#'    plot(icmstate_pt, id = 2)
#' }
#' 
#' 
#' @export
summary.probtrans.subjects <- function(object, id, times, from=1, to=0,
                              extend=FALSE, ...)
{
  if(!inherits(object, "probtrans.subjects")){
    stop("Object to summarize must be of class 'probtrans.subjects'.")
  }
  
  if(length(id) > 1){
    stop("Argument 'id' must be at most of length 1.")
  }
  
  which_id <- which(object[["subject_ids"]] == id)  
  
  if(length(which_id) == 1){
    res <- summary(object[[which_id]], variance = FALSE, extend = extend, ...)
    attr(res, "id") <- id
    class(res) <- "summary.probtrans.subjects"
    return(res)
  } else{
    stop("Subject 'id' not found in 'probtrans.subjects' object. 
         Please provide valid 'id'.")
  }
}




#' Print method for a summary.probtrans.subjects object
#' 
#' @param x Object of class 'summary.probtrans.subjects', to be printed
#' @param complete Whether or not the complete estimated transition
#' probabilities should be printed (\code{TRUE}) or not (\code{FALSE}); default
#' is \code{FALSE}, in which case the estimated transition probilities will be
#' printed for the first and last 6 time points of each starting state or of
#' the selected times (or all when there are at most 12 of these time points
#' @param \dots Further arguments to print
#' 
#' @aliases print.summary.probtrans.subjects
#' 
#' @method print summary.probtrans.subjects
#' 
#' @import mstate
#' 
#' @examples 
#' 
#' \dontrun{
#' # If all time points should be printed, specify complete=TRUE in the print statement
#' print(x, complete=TRUE)
#' }
#' 
#' @export 
print.summary.probtrans.subjects <- function(x, complete=FALSE, ...)
{
  cat("\nPrediction for subject", attr(x, "id"), "\n")
  class(x) <- "summary.probtrans"
  print(x, ...)
  return(invisible())
}



#' Print method for a summary.probtrans object
#' 
#' @param x Object of class 'summary.probtrans', to be printed
#' @param complete Whether or not the complete estimated transition
#' probabilities should be printed (\code{TRUE}) or not (\code{FALSE}); default
#' is \code{FALSE}, in which case the estimated transition probilities will be
#' printed for the first and last 6 time points of each starting state or of
#' the selected times (or all when there are at most 12 of these time points
#' @param \dots Further arguments to print
#' 
#' @aliases print.summary.probtrans
#' 
#' @importFrom utils head
#' 
#' @method print summary.probtrans
#' 
#' @examples 
#' 
#' \dontrun{
#' # If all time points should be printed, specify complete=TRUE in the print statement
#' print(x, complete=TRUE)
#' }
#' 
#' @export 
print.summary.probtrans <- function(x, complete=FALSE, ...)
{
  if (!inherits(x, "summary.probtrans"))
    stop("'x' must be a 'summary.probtrans' object")
  from <- attr(x, "from")
  tt <- unique(x[[1]]$time) # the time points
  nt <- length(tt)
  if (nt<=12 | complete) {
    for (k in from) {
      cat("\nPrediction from state", k, ":\n")
      ptk <- x[[k]]
      print(ptk, ...)
    }
  } else {
    for (k in from) {
      cat("\nPrediction from state", k, "(head and tail):\n")
      ptk <- x[[k]]
      print(head(ptk), ...)
      cat("\n...\n")
      print(tail(ptk), ...)
    }
  }
  return(invisible())
}


