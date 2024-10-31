## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(icmstate)
set.seed(1)

## ----transmatID---------------------------------------------------------------
library(mstate)
tmat_ID <- transMat(list(c(2, 3), c(3), c()), names = c("Alive", "Illness",
                                                        "Death"))

## ----printransmat-------------------------------------------------------------
tmat_ID

## ----shapepars----------------------------------------------------------------
#The first entry corresponds to shape of first transition, and so on...
shape_ID <- c(1, 1.5, 0.5)

## ----scalepars----------------------------------------------------------------
scale_ID <- c(6, 10/gamma(1+1/1.5), 10/gamma(1+1/0.5))

## ----obsdata------------------------------------------------------------------
obs_data <- data.frame(time = c(rep(seq(0, 20, 2), 50), rep(seq(1, 19, 2), 50)),
                       id = c(rep(1:50, each = 11), rep(51:100, each = 10)))

## ----gendatamanual------------------------------------------------------------
data_manual <- sim_weibmsm(data = obs_data, tmat = tmat_ID,
                           shape = shape_ID, scale = scale_ID)

## ----vismanualdata, dpi=72----------------------------------------------------
visualise_msm(data_manual)

## ----autodatasim--------------------------------------------------------------
data_auto <- sim_weibmsm(tmat = tmat_ID, shape = shape_ID, scale = scale_ID,
                         n_subj = 100, obs_pars = c(2, 0.5, 20))

## ----visautodata, dpi = 72----------------------------------------------------
visualise_msm(data_auto)

## ----visstartprobs, dpi=72----------------------------------------------------
data_sprobs <- sim_weibmsm(tmat = tmat_ID, shape = shape_ID, scale = scale_ID,
                         n_subj = 40, obs_pars = c(2, 0.5, 20), 
                         startprobs = c(0.6, 0.3, 0.1))
visualise_msm(data_sprobs)

## ----visexact, dpi=72---------------------------------------------------------
data_exact <- sim_weibmsm(tmat = tmat_ID, shape = shape_ID, scale = scale_ID,
                         n_subj = 40, obs_pars = c(2, 0.5, 20), 
                         exact = c(3))
visualise_msm(data_exact)

## ----censoringdat, dpi=72-----------------------------------------------------
data_cens <- sim_weibmsm(tmat = tmat_ID, shape = shape_ID, scale = scale_ID,
                         n_subj = 40, obs_pars = c(2, 0.5, 20), 
                         censshape = c(1, NA, NA), censscale = c(3, NA, NA))
visualise_msm(data_cens)

## ----truetrajdata, dpi=72-----------------------------------------------------
data_true <- sim_weibmsm(tmat = tmat_ID, shape = shape_ID, scale = scale_ID,
                         n_subj = 40, obs_pars = c(2, 0.5, 20), 
                         true_trajec = TRUE)
visualise_msm(data_true$observed)
visualise_msm(data_true$true)

## ----simID--------------------------------------------------------------------
data_ID <- sim_weibmsm(tmat = tmat_ID, shape = shape_ID, scale = scale_ID,
                       n_subj = 2000, obs_pars = c(2, 0.5, 20), 
                       true_trajec = TRUE)

## ----survivalcheck------------------------------------------------------------
library(survival)
tmat2_ID <- to.trans2(tmat_ID)
dat_true <- data_ID$true

## ----truevssimulatedtrajectories----------------------------------------------
opar <- par(no.readonly = TRUE)
par(mfrow = c(2,2))
dat_surv <- reshape(dat_true, direction = "wide", idvar = "id", timevar = "state")
dat_surv$status <- rep(1, nrow(dat_surv))

#From 1 to 2
dat1 <- dat_surv
dat1[is.na(dat1[, "time.2"]), "status"] <- 0
dat1[is.na(dat1[, "time.2"]), "time.2"] <- dat1[is.na(dat1[, "time.2"]), "time.3"]
dat1 <- dat1[!dat1[, "time.2"] == 0,]
first_trans <- survfit(Surv(time.1, time.2, status) ~ 1, data = dat1)
plot(first_trans, fun = "cumhaz", xlim = c(0, 30), main = "First transition")
lines(first_trans$time, -pweibull(first_trans$time, shape_ID[1], scale_ID[1], lower = FALSE, log = TRUE), col = "red")

#From 1 to 3
dat2 <- dat_surv
dat2[!is.na(dat2[, "time.2"]), "status"] <- 0
dat2[!is.na(dat2[, "time.2"]), "time.3"] <- dat2[!is.na(dat2[, "time.2"]), "time.2"]
second_trans <- survival::survfit(Surv(time.3, status) ~ 1, data = dat2)
plot(second_trans, fun = "cumhaz", xlim = c(0, 30), main = "Second transition")
lines(second_trans$time, -pweibull(first_trans$time, shape_ID[2], scale_ID[2], lower = FALSE, log = TRUE), col = "red")

third_trans <- survfit(Surv(time.2, time.3, status) ~ 1, data = subset(dat_surv, !is.na(time.2)))
plot(third_trans, fun = "cumhaz", xlim = c(0, 30), main = "Third transition")
lines(third_trans$time, -pweibull(third_trans$time, shape_ID[3], scale_ID[3], lower = FALSE, log = TRUE), col = "red")
par(opar)

## ----npmsmicdata--------------------------------------------------------------
EM_fit <- npmsm(subset(data_ID$observed, id < 100), tmat = tmat_ID, tol = 1e-4, 
                maxit = 30)
plot(EM_fit)
tseq <- seq(0, 20, 0.01)
lines(tseq, (tseq/scale_ID[1])^shape_ID[1], col = "black", lty = 2)
lines(tseq, (tseq/scale_ID[2])^shape_ID[2], col = "red", lty = 2)
lines(tseq, (tseq/scale_ID[3])^shape_ID[3], col = "green", lty = 2)

## -----------------------------------------------------------------------------
par(mfrow = c(1, 2))
#Plot true transition probabilities
plot(probtrans_weib(tmat_ID, seq(0, 20, 0.01), shapes = shape_ID, 
                    scales = scale_ID))
#Plot estimated transition probabilities
plot(transprob(EM_fit, predt = 0))
par(opar)

