% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_tp.R
\name{predict_tp}
\alias{predict_tp}
\title{Calculate subject specific transition probabilities from a 
multi-state proportional hazards model.}
\usage{
predict_tp(
  object,
  predt,
  direction = c("forward", "fixedhorizon"),
  newdata,
  trans
)
}
\arguments{
\item{object}{A \code{\link[survival:coxph]{coxph}} object fit on multi-state
data. Must contain a \code{strata(X)} term. Data used for the coxph() fit preferably prepared 
using \code{\link[mstate:msprep]{msprep}}.}

\item{predt}{A positive number indicating the prediction time. This is either 
the time at which the prediction is made (if \code{direction = "forward"}) or 
the time for which the prediction is to be made (if \code{direction = "backward"}).}

\item{direction}{One of \code{"forward"} (default) or \code{"fixedhorizon"},
indicating whether prediction is forward or for a fixed horizon}

\item{newdata}{A \code{data.frame} containing a single row per
subject in the data. It must contain the following named columns:
\describe{
  \item{\code{id}:}{(optional) Unique identifier of the subject, can be numeric or character;}
  \item{\code{"variables"}:}{The variables and their values for the subject(s)
  (optionally identified by "id"). After using \code{\link[mstate:expand.covs]{expand.covs}}, 
  the names of the new data must
  match the names of the variables in the coxph \code{object}.}
} Note that newdata must contain a column containing the variable which was 
used to determine the stratum of a transition in \code{object}. 
Usually the stratum is determined from a column that is generated automatically 
in the pre-processing steps of this function.}

\item{trans}{A transition matrix as created by \code{\link[mstate:transMat]{transMat}}.}
}
\value{
An object of class \code{"probtrans.subjects"}.
This is a list of length n (number of subjects in newdata), with each list element
an object of class \code{\link[mstate:probtrans]{probtrans}} for the associated 
subject. List elements can be accessed using \code{[[x]]}, with \code{x} 
ranging from 1 to n. Additionally, each list element 
has an element \code{$id}, representing the subject id and the output object 
also has an element \code{$subject_ids} representing the subject ids in order.
}
\description{
Given a coxph model fit on multi-state data (prepared with \code{\link[mstate:msprep]{msprep}}),
determine transition probabilities for subjects in \code{newdata}.
}
\details{
When using this function for \code{newdata} with many subjects, consider running the 
function multiple times for parts of \code{newdata} to negate the risk of 
running our of memory.
Using this function, it is only possible to consider models with transition specific covariates.
If you would like to have covariates shared over transitions or proportional 
hazards assumptions between transitions, see \code{\link{probtrans_coxph}}.
}
\examples{
#Example from the mstate vignette
#We determine the subject specific transition probabilities for subjects
#in the ebmt3 data-set 
if(require("mstate")){
  data(ebmt3)
  n <- nrow(ebmt3)
  tmat <- transMat(x = list(c(2, 3), c(3), c()), names = c("Tx",
                                                           "PR", "RelDeath"))
  #From days to years                                                          
  ebmt3$prtime <- ebmt3$prtime/365.25
  ebmt3$rfstime <- ebmt3$rfstime/365.25
  #Covariates we will use
  covs <- c("dissub", "age", "drmatch", "tcd", "prtime")
  msbmt <- msprep(time = c(NA, "prtime", "rfstime"), status = c(NA,
                  "prstat", "rfsstat"), data = ebmt3, trans = tmat, keep = covs)
  #Expand covariates so that we can have transition specific covariates
  msbmt2 <- expand.covs(msbmt, covs, append = TRUE, longnames = FALSE)
  
  #-------------Model---------------------#
  
  #Simple model, transition specific covariates, each transition own baseline hazard
  c1 <- coxph(Surv(Tstart, Tstop, status) ~ dissub1.1 + dissub2.1 +
                age1.1 + age2.1 + drmatch.1 + tcd.1 + dissub1.2 + dissub2.2 +
                age1.2 + age2.2 + drmatch.2 + tcd.2 + dissub1.3 + dissub2.3 +
                age1.3 + age2.3 + drmatch.3 + tcd.3 + strata(trans), data = msbmt2,
                method = "breslow")
  
   #Predict transition probabilities for first 30 subjects.
   tp_subjects <- predict_tp(c1, predt = 0, direction = "forward", 
                                  newdata = ebmt3[1:30,], trans = tmat)
                                  
   #Now we can plot the transition probabilities for each subject separately:
   plot(tp_subjects, id = 1)
   #tp_subjects has length number of subjects in newdata + 1
   #And tp_subjects[[i]] is an object of class "probtrans", so you can 
   #use all probtrans functions: summary, plot etc.
}


}
\seealso{
\code{\link{probtrans_coxph}}, \code{\link{plot.probtrans.subjects}},
\code{\link{summary.probtrans.subjects}}
}
