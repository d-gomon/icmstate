## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----IEDgeneration, warning=FALSE---------------------------------------------
library(icmstate)
library(mstate)
library(msm)
set.seed(1)
tmat_EID <- mstate::transMat(x = list( c(2, 3), c(4), c(), c() ))
qmatrix <- rbind(
    c(-0.15, 0.1, 0.05, 0),
    c(0, -0.1, 0, 0.1),
    c(0, 0, 0, 0),
    c(0, 0, 0, 0)
  )
n <- 30

#time = observation time, subject = subject identifier
simdat <- data.frame(time = c(replicate(n, c(0, seq(2, 12, by=2) + runif(6, 0, 2)))),
                    subject = rep(1:n, each = 7))
#Simulate interval-censored data. See help(simmulti.msm)
dat <- simmulti.msm(data = simdat, qmatrix = qmatrix, start = 1,
                    death = c(3,4))[, 1:3]
names(dat)[1] <- "id"


## ----visualiseIED-------------------------------------------------------------
visualise_msm(dat, tmat = tmat_EID)

## ----remove14trans------------------------------------------------------------
dat <- dat[!dat[, "id"] %in% c(5, 22),]

## ----msmtoFrydman-------------------------------------------------------------
#Create Frydman data
msmtoFrydman <- function(gd){
  #Create Frydman data
  gd <- remove_redundant_observations(gd, tmat = tmat_EID)
  gd_frydman <- NULL
  for(j in unique(gd$id)){
    tempdat <- subset(gd, id == j)
    tempstates <- unique(tempdat$state)
    if(length(tempstates) == 1){ #If we only observe the subject in 1 state, right censored in 1
      gdi_frydman <- data.frame(delta = 0, Delta = 0,
                                L = NA,
                                R = NA,
                                time = tempdat$time[length(tempdat$time)])
    } else if(length(tempstates) == 2){ #If we only observe the subject in 2 states, either 1->2 or 1->3 has occured
      if(all(tempstates %in% c(1,2))){
        gdi_frydman <- data.frame(delta = 1, Delta = 0, 
                                  L = tempdat$time[which.min(tempdat$state == 1)-1],
                                  R = tempdat$time[which.min(tempdat$state == 1)],
                                  time = tempdat$time[length(tempdat$time)])
      } else if(all(tempstates %in% c(1,3))){
        gdi_frydman <- data.frame(delta = 0, Delta = 1, 
                                  L = NA,
                                  R = NA,
                                  time = tempdat$time[which.min(tempdat$state == 1)])
      }
    } else if(length(tempstates) == 3){ #If we observe 3 states, then 1->2->3 must have occured
      gdi_frydman <- data.frame(delta = 1, Delta = 1, 
                                L = tempdat$time[which.min(tempdat$state == 1)-1],
                                R = tempdat$time[which.min(tempdat$state == 1)],
                                time = tempdat$time[length(tempdat$time)])
    }
    gd_frydman <- rbind(gd_frydman, gdi_frydman)
  }
  return(gd_frydman)
}

## ----visFrydmandat------------------------------------------------------------
dat_frydman <- msmtoFrydman(dat)
visualise_data(dat_frydman)

## ----fitmodels----------------------------------------------------------------
mod_npmsm <- npmsm(gd = dat, tmat = tmat_EID, maxit = 300, exact = c(3,4),
                   tol = 1e-6)
mod_frydman <- msm_frydman(data = dat_frydman, tol = 1e-6)

## ----printsupportmodels-------------------------------------------------------
supp_npmsm <- support_npmsm(mod_npmsm, cutoff = 1e-9)
supp_frydman <- mod_frydman$supportMSM$Q_mat
print("Frydman Support")
supp_frydman
print("icmstate Support")
supp_npmsm$`State 1 -> State 2`$support

## -----------------------------------------------------------------------------
#Right-endpoints of the 1->2 transition
RE12 <- supp_frydman[, 2]
#Right-endpoints of the 1->3 transition
RE13 <- mod_frydman$data_idx$e_k_star
#Right-endpoints of the 2->3 transition
RE23 <- mod_frydman$data_idx$t_n_star

## ----compsurvstate1-----------------------------------------------------------
#Transition probabilities from state 1 from time 0
P11 <- transprob(mod_npmsm, predt = 0)
#Extract times of interest
times1 <- P11[[1]][,1]
#1-F(s) for Frydman estimator:
#We take min(F_{12}(x)) as the cdf has only jumped at the right-endpoints
Frydman1minF <- sapply(times1, function(x) 1- (mod_frydman$cdf$F13(x) +
                                               min(mod_frydman$cdf$F12(x))))
#Comparison plot
plot(P11, main = "icmstate vs Frydman (dashed red line) 
     Comparison times (dotted blue lines)")
lines(times1, Frydman1minF, col = "red", type = "s", lwd = 2, lty = 2)
abline(v = unique(c(RE12, RE13)), col = "blue", lwd = 2, lty = 3)

## -----------------------------------------------------------------------------
FrydmanF13 <- sapply(times1, function(x) mod_frydman$cdf$F13(x))
#Comparison plot
plot(P11, main = "icmstate vs Frydman (dashed red line) 
     Comparison times (dotted blue lines)", ord = c(3, 1, 2, 4))
lines(times1, FrydmanF13, col = "red", type = "s", lwd = 2, lty = 2)
abline(v = unique(RE13), col = "blue", lwd = 2, lty = 3)

## -----------------------------------------------------------------------------
FrydmanF12 <- sapply(times1, function(x) min(mod_frydman$cdf$F12(x)))
#Comparison plot
plot(P11, main = "icmstate vs Frydman (dashed red line) 
     Comparison times (dotted blue lines)", ord = c(3, 1, 2, 4))
lines(times1, 1-FrydmanF12, col = "red", type = "s", lwd = 2, lty = 2)
abline(v = unique(RE12), col = "blue", lwd = 2, lty = 3)

## -----------------------------------------------------------------------------
#Calculate dA23
FrydmandA23 <- c(0, diff(sapply(times1, function(x) mod_frydman$cdf$Lambda23(x))))
#We calculate the product integral for the Frydman estimator
FrydmanG <- cumprod(1-FrydmandA23)
#Comparison plot
plot(P11, main = "icmstate vs Frydman (dashed red line) 
     Comparison times (dotted blue lines)", from = 2, ord = c(2, 4, 1, 3))
lines(times1, FrydmanG, col = "red", type = "s", lwd = 2, lty = 2)
abline(v = unique(RE23), col = "blue", lwd = 2, lty = 3)

