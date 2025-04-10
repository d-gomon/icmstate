#' @title
#' Plot an object of class \code{"probtrans.subjects"}
#' 
#' @description
#' Plots the transition probabilities for a specific subject. Wrapper for 
#' \code{\link[mstate:plot.probtrans]{plot.probtrans}}
#' 
#' @param x An object of class \code{"probtrans.subjects"}
#' @param id Subject identifier as in \code{newdata} of
#'  \code{\link[icmstate:probtrans_coxph]{probtrans_coxph}}
#' @param ... Further arguments to \code{\link[mstate:plot.probtrans]{plot.probtrans}}
#' 
#' @details
#' Note that 
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
    stop("Subject 'id'(s) not found in 'probtrans.subjects' object. 
         Please provide valid 'id'(s).")
  }
}




