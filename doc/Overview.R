## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = 'figure-latex/', 
  cache = FALSE
)

## ----setup, echo=FALSE, warning = FALSE---------------------------------------
library(icmstate)
library(msm)
library(latex2exp)
set.seed(1)

## ----echo = FALSE, fig.keep='none', warning=FALSE-----------------------------
library(ggplot2)
library(latex2exp)
tmat <- mstate::trans.illdeath()

eval_times <- function(n_obs, stop_time){
  cumsum( c( 0,  runif( n_obs-1, 0, 2*(stop_time-4)/(n_obs-1) ) ) )
}
set.seed(1)
sim_dat <- sim_id_weib(n = 5, n_obs = 6, stop_time = 15, eval_times = eval_times,
                       start_state = "stable", shape = c(0.5, 0.5, 2), scale = c(5, 10, 10/gamma(1.5)))

ID_plot <- invisible(visualise_msm(remove_redundant_observations(sim_dat, tmat = tmat))$plot + xlab("Subject") + ylab("Time") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold")) +
  ggtitle(label = "Illness-Death Model Data"))


unique_times <- unique(sort(c(ID_plot$data$t1, ID_plot$data$t2)))
n_unique <- length(unique_times)
# Create the vector of LaTeX expressions
tau_labels <- sapply(1:n_unique, function(i) TeX(paste0("$\\tau_{", i, "}$")))

ID_plot3 <- ID_plot + 
  geom_segment(aes(x = 0, xend = ID_plot$data$ID, y = ID_plot$data$t1, yend = ID_plot$data$t1), linetype = "dashed", col = "red", lwd = 1.3) + 
  geom_segment(aes(x = 0, xend = ID_plot$data$ID, y = ID_plot$data$t2, yend = ID_plot$data$t2), linetype = "dashed", col = "red", lwd = 1.3) + 
  scale_y_continuous(breaks = unique(sort(c(ID_plot$data$t1, ID_plot$data$t2))), labels = tau_labels)

## ----IDplot, echo = FALSE, fig.cap = "Visualisation of 5 subjects following an illness-death model. Each line represents the trajectory of a single subject, with the numbers indicating the state the subject was observed in at a certain time.", fig.align='center', out.width='70%'----
ID_plot 

## ----IDplot3, echo = FALSE, fig.cap = "Unique observation times ordered chronologically and named.", fig.align='center', out.width='70%'----
ID_plot3 

## -----------------------------------------------------------------------------
library(mstate)
tmat <- transMat(x = list( c(2), c() ))

## -----------------------------------------------------------------------------
library(msm)
#Entry: constant transition rate from row state to column state
#Absorbing state 2
qmatrix <- rbind(
  c(-0.2, 0.2),
  c(0, 0)
)
#Number of subjects in simulated data
n <- 10
#time = observation time, subject = subject identifier
simdat <- data.frame(time = c(replicate(n, c(0, seq(2, 12, by=2) + runif(6, 0, 2)))),
                    subject = rep(1:n, each = 7))
#Simulate interval-censored data. See help(simmulti.msm)
dat <- simmulti.msm(data = simdat, qmatrix = qmatrix, start = 1)[, 1:3]
names(dat)[1] <- "id"

## -----------------------------------------------------------------------------
visualise_msm(dat)

## -----------------------------------------------------------------------------
icdata <- NULL
for(i in unique(dat$id)){
  gdi <- subset(dat, id == i)
  L_idx <- Position(isTRUE, gdi$state <= 1, right = TRUE)
  L <- gdi$time[L_idx] #Exit time from state 1
  if(L_idx < nrow(gdi)){
    R <- gdi$time[L_idx + 1]
  } else{
    R <- Inf
  }
  icdata <- rbind(icdata, c(L, R))
}
icdata <- as.data.frame(icdata)
colnames(icdata) <- c("L", "R")

## ----message=FALSE, warning=FALSE---------------------------------------------
library(icenReg)
SimpleMult <- npmsm(gd = dat, tmat = tmat, method = "multinomial")
Turnbull <- ic_np(cbind(L, R) ~ 0, data = icdata)

## ----warning = FALSE----------------------------------------------------------
#Extract Survival function from Turnbull estimate
surv_turnbull <- sapply(1:length(Turnbull$p_hat), function(i) 1 - sum(Turnbull$p_hat[1:i]))
#Extract Survival function using transition probabilities
surv_npmsm <- transprob(SimpleMult, predt = 0, direction = "forward")
#Plot the result
plot(surv_npmsm, main = "Comparison of survival functions: black (npmsm), red (Turnbull)")
lines(c(0, Turnbull$T_bull_Intervals[2,]), c(1, surv_turnbull), type = "s", col = "red",
      lwd = 2, lty = 2)

## -----------------------------------------------------------------------------
#Absorbing state 3, transition intensities constant = 0.1
qmatrix <- rbind(
  c(-0.2, 0.1, 0.1),
  c(0, -0.1, 0.1),
  c(0, 0, 0)
)
#Create a transition matrix
tmat_ID <- trans.illdeath()
#Number of subjects in simulated data
n <- 100
#Create data frame for simulation:
simdat_ID <- data.frame(time = c(replicate(n, c(0, seq(2, 12, by=2) + runif(6, 0, 2)))),
                    subject = rep(1:n, each = 7))
dat_ID <- simmulti.msm(data = simdat_ID, qmatrix = qmatrix, start = 1)[, 1:3]
names(dat_ID)[1] <- "id"

## ----visID, fig.cap="Visualisation of 20 subjects in an illness-death model"----
visualise_msm(subset(dat_ID, id < 21))

## -----------------------------------------------------------------------------
#Fit appropriate MSM using the msm package
ID_msm <- msm(state ~ time, data = dat_ID, subject = id, qmatrix = qmatrix)
#Fit NP MSM using the icmstate package
ID_npmsm <- npmsm(gd = dat_ID, tmat = tmat_ID)

## ----IDfits, fig.cap = "Cumulative intensity estimates using the msm package and icmstate package. True cumulative intensities are straight lines with slope $0.1$ for each transition."----
plot(ID_npmsm, main = "Cumulative intensity: icmstate (solid lines) vs msm (dashed lines)")
#To plot output from the 'msm' function we need to extract the estimates
cols <- c("black", "red", "green")
for(i in 1:length(ID_msm$estimates)){
  abline(a = 0, b = exp(ID_msm$estimates)[i], col = cols[i], lty = 2)
}


## -----------------------------------------------------------------------------
ID_npmsm

## -----------------------------------------------------------------------------
plot(transprob(object = ID_npmsm, predt = 0, direction = "forward", variance = FALSE),
     main = "Transition probabilities from npmsm() fit")

## -----------------------------------------------------------------------------
plot(transprob(ID_msm, predt = 0, times = seq(0, 14, 0.05)), 
     main = "Transition probabilities from msm() fit")

## ----gendata_hom_exact--------------------------------------------------------
#Absorbing state 3, transition intensities constant = 0.1
qmatrix <- rbind(
  c(-0.2, 0.1, 0.1),
  c(0, -0.1, 0.1),
  c(0, 0, 0)
)
#Create a transition matrix
tmat_ID <- trans.illdeath()
#Number of subjects in simulated data
n <- 100
#Create data frame for simulation:
simdat_ID_exact <- data.frame(time = c(replicate(n, c(0, seq(2, 12, by=2) + runif(6, 0, 2)))),
                    subject = rep(1:n, each = 7))
dat_ID_exact <- simmulti.msm(data = simdat_ID_exact, qmatrix = qmatrix, start = 1, 
                             death = 3)[, 1:3]
names(dat_ID_exact)[1] <- "id"

## ----fitexactmodels, warning = FALSE------------------------------------------
#Fit appropriate MSM using the msm package
ID_msm_exact <- msm(state ~ time, data = dat_ID_exact, subject = id, qmatrix = qmatrix, 
                    deathexact = 3)
#Fit NP MSM using the icmstate package
ID_npmsm_exact <- npmsm(gd = dat_ID_exact, tmat = tmat_ID, exact = 3)

## ----IDfitsexact, fig.cap = "Cumulative intensity estimates using the msm package and icmstate package. True cumulative intensities are straight lines with slope $0.1$ for each transition."----
plot(ID_npmsm_exact, main = "Cumulative intensity (exactly observed death):\n icmstate (solid lines) vs msm (dashed lines)")
cols <- c("black", "red", "green")
for(i in 1:length(ID_msm_exact$estimates)){
  abline(a = 0, b = exp(ID_msm_exact$estimates)[i], col = cols[i], lty = 2)
}

## ----plottransprobexactnpmsm, warning=FALSE-----------------------------------
plot(transprob(object = ID_npmsm_exact, predt = 0, direction = "forward", variance = FALSE),
     main = "Transition probabilities from npmsm() fit (exactly observed death)")

## ----plottransprobexactmsm, warning=FALSE-------------------------------------
plot(transprob(ID_msm_exact, predt = 0, times = seq(0, 14, 0.05)), 
     main = "Transition probabilities from msm() fit (exactly observed death)")

## ----simdatinhomogeneous------------------------------------------------------
#When do we observe the subjects? n_obs times with uniform inter observation times
eval_times <- function(n_obs, stop_time){
  cumsum( c( 0,  runif( n_obs-1, 0, 2*(stop_time-4)/(n_obs-1) ) ) )
}

#Simulate data with Weibull shapes and scales.
dat_inhom <- sim_id_weib(n = 50, n_obs = 6, stop_time = 15, eval_times = eval_times,
start_state = "stable", shape = c(0.5, 0.5, 2), scale = c(5, 10, 10/gamma(1.5)))

## -----------------------------------------------------------------------------
mult_fit <- npmsm(gd = dat_inhom, tmat = tmat_ID)
pois_fit <- npmsm(gd = dat_inhom, tmat = tmat_ID, method = "poisson")

## -----------------------------------------------------------------------------
plot(mult_fit, main = "Cumulative Intensities for Weibull illness-death model:  
     Multinomial (solid), Poisson (dotted), True (dashed)")
shapes <- c(0.5, 0.5, 2)
scales <- c(5, 10, 10/gamma(1.5))
for(i in 1:3){
  trans_dat <- subset(pois_fit$A$Haz, trans == i)
  #Add poisson estimates
  lines(trans_dat$time, trans_dat$Haz, lty = 3, col = cols[i])
  #Add true cumulative intensities
  lines(trans_dat$time, 
        -pweibull(trans_dat$time, shapes[i], scales[i], lower = FALSE, log = TRUE),
        lty = 2, col = cols[i])
}

## ----warning=FALSE------------------------------------------------------------
#Fit NP MSM using the icmstate package
dat_ID_small <- subset(dat_ID, id < 21)
ID_npmsm_equalprob <- npmsm(gd = dat_ID_small, tmat = tmat_ID, maxit = 300)
ID_npmsm_hom <- npmsm(gd = dat_ID_small, tmat = tmat_ID, inits = "homogeneous", maxit = 300)
ID_npmsm_unif <- npmsm(gd = dat_ID_small, tmat = tmat_ID, inits = "unif", maxit = 300)
ID_npmsm_beta <- npmsm(gd = dat_ID_small, tmat = tmat_ID, inits = "beta", beta_params = c(3, 6), maxit = 300)

## -----------------------------------------------------------------------------
fit_summary <- sapply(list(ID_npmsm_equalprob, ID_npmsm_hom, ID_npmsm_unif, ID_npmsm_beta), 
       function(x) c(x$it, x$ll))
attr(fit_summary, "dimnames") <- list("summary" = c("iterations", "likelihood"),
                                      "Initial" = c("EqProb", "Hom", "Unif", "Beta"))
fit_summary

