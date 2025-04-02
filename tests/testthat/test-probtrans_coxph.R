#-----To-Do corner------#
#in expand.covs(), we now assume that longnames = FALSE
#Could also account for longnames = TRUE, if we were to check the names for both expansions
#This might be a problem if the data set is quite large
#Maybe check the names first for a very small data set and then decide which one to use
#-----------------------#



#------------------Data generation for tests---------------------#
#First we create the data that we will use for the tests
#For this, we largely follow the vignette of mstate.
library(mstate)
data(ebmt3)

n <- nrow(ebmt3)

tmat <- transMat(x = list(c(2, 3), c(3), c()), names = c("Tx",
                                                         "PR", "RelDeath"))
ebmt3$prtime <- ebmt3$prtime/365.25
ebmt3$rfstime <- ebmt3$rfstime/365.25
covs <- c("dissub", "age", "drmatch", "tcd", "prtime")
msbmt <- msprep(time = c(NA, "prtime", "rfstime"), status = c(NA,
                "prstat", "rfsstat"), data = ebmt3, trans = tmat, keep = covs)

#For now, let's say that we always use longnames = FALSE
msbmt <- expand.covs(msbmt, covs, append = TRUE, longnames = FALSE)

#Create extra variable to allow gender mismatch to have the same effect for transitions 2 and 3.
msbmt$drmatch.2.3 <- msbmt$drmatch.2 + msbmt$drmatch.3

#Introduce pr covariate for proportionality assumption of transitions 2 and 3
msbmt$pr <- 0
msbmt$pr[msbmt$trans == 3] <- 1



#-------------Models---------------------#

#Simple model, transition specific covariates, each transition own baseline hazard
c1 <- coxph(Surv(Tstart, Tstop, status) ~ dissub1.1 + dissub2.1 +
              age1.1 + age2.1 + drmatch.1 + tcd.1 + dissub1.2 + dissub2.2 +
              age1.2 + age2.2 + drmatch.2 + tcd.2 + dissub1.3 + dissub2.3 +
              age1.3 + age2.3 + drmatch.3 + tcd.3 + strata(trans), data = msbmt,
              method = "breslow")

#Model with same baseline hazards for transitions 2 (1->3) and 3(2->3)
#pr then gives the ratio of the 2 hazards for these transitions
c2 <- coxph(Surv(Tstart, Tstop, status) ~ dissub1.1 + dissub2.1 +
              age1.1 + age2.1 + drmatch.1 + tcd.1 + dissub1.2 + dissub2.2 +
              age1.2 + age2.2 + drmatch.2 + tcd.2 + dissub1.3 + dissub2.3 +
              age1.3 + age2.3 + drmatch.3 + tcd.3 + pr + strata(to), data = msbmt,
              method = "breslow")

# #Model which evaluates whether time of arrival in 2 influences transition 3 (2->3)
# #This is achieved by including prtime.3, as prtime is the time of arrival in 2.
# c3 <- coxph(Surv(Tstart, Tstop, status) ~ dissub1.1 + dissub2.1 +
#               age1.1 + age2.1 + drmatch.1 + tcd.1 + dissub1.2 + dissub2.2 +
#               age1.2 + age2.2 + drmatch.2 + tcd.2 + dissub1.3 + dissub2.3 +
#               age1.3 + age2.3 + drmatch.3 + tcd.3 + pr + prtime.3 + strata(to),
#               data = msbmt, method = "breslow")

#Same as c2, but now Gender mismatch has the same effect for both 
c4 <- coxph(Surv(Tstart, Tstop, status) ~ dissub1.1 + dissub2.1 +
              age1.1 + age2.1 + drmatch.1 + tcd.1 + dissub1.2 + dissub2.2 +
              age1.2 + age2.2 + drmatch.2.3 + tcd.2 + dissub1.3 + dissub2.3 +
              age1.3 + age2.3  + tcd.3 + pr + strata(to), data = msbmt,
            method = "breslow")




test_that("Do we recover the same baseline hazards as msfit?", {
  #Create new data to recover baseline hazard
  newd <- data.frame(dissub = rep(0, 3), age = rep(0, 3), drmatch = rep(0,
                     + 3), tcd = rep(0, 3), trans = 1:3)
  newd$dissub <- factor(newd$dissub, levels = 0:2, labels = levels(ebmt3$dissub))
  newd$age <- factor(newd$age, levels = 0:2, labels = levels(ebmt3$age))
  newd$drmatch <- factor(newd$drmatch, levels = 0:1, labels = levels(ebmt3$drmatch))
  newd$tcd <- factor(newd$tcd, levels = 0:1, labels = levels(ebmt3$tcd))
  attr(newd, "trans") <- tmat
  class(newd) <- c("msdata", "data.frame")
  newd <- expand.covs(newd, covs[1:4], longnames = FALSE)
  newd$strata = 1:3
  
  #Check for first model
  msf1 <- msfit(c1, newdata = newd, trans = tmat)
  
  mstate_pt <- probtrans(msf1, predt = 0, direction = "forward")
  icmstate_pt <- probtrans_coxph(c1, predt = 0, direction = "forward", 
                                 newdata = newd, trans = tmat)
  
  expect_equal(mstate_pt[[1]][, 1:4], icmstate_pt[[1]][[1]])
  expect_equal(mstate_pt[[2]][, 1:4], icmstate_pt[[1]][[2]])
  expect_equal(mstate_pt[[3]][, 1:4], icmstate_pt[[1]][[3]])
  
  
  #Check for second model
  newd$strata = c(1, 2, 2)
  newd$pr <- c(0, 0, 1)
  msf2 <- msfit(c2, newdata = newd, trans = tmat)
  
  mstate_pt <- probtrans(msf2, predt = 0, direction = "forward")
  icmstate_pt <- probtrans_coxph(c2, predt = 0, direction = "forward", 
                                 newdata = newd, trans = tmat)
  
  expect_equal(mstate_pt[[1]][, 1:4], icmstate_pt[[1]][[1]])
  expect_equal(mstate_pt[[2]][, 1:4], icmstate_pt[[1]][[2]])
  expect_equal(mstate_pt[[3]][, 1:4], icmstate_pt[[1]][[3]])
  
  #Check for fourth model
  newd$drmatch.2.3 <- newd$drmatch.2 + newd$drmatch.3
  msf4 <- msfit(c4, newdata = newd, trans = tmat)
  
  mstate_pt <- probtrans(msf4, predt = 0, direction = "forward")
  icmstate_pt <- probtrans_coxph(c4, predt = 0, direction = "forward", 
                                 newdata = newd, trans = tmat)
  
  expect_equal(mstate_pt[[1]][, 1:4], icmstate_pt[[1]][[1]])
  expect_equal(mstate_pt[[2]][, 1:4], icmstate_pt[[1]][[2]])
  expect_equal(mstate_pt[[3]][, 1:4], icmstate_pt[[1]][[3]])
})







test_that("Do we recover the same transition probabilities as msfit + probtrans?", {
  #Let us store P_{12}(6) for the first 30 subjects in ebmt3
  #Using cox model 2 & 4
  n <- 30
  
  #Using msfit and probtrans from mstate pkg
  #We need to make predictions for each person separately.
  ttmat <- to.trans2(tmat)[, c(2, 3, 1)]
  names(ttmat)[3] <- "trans"
  p12_6 <- rep(NA, n)
  for (j in 1:n) {
    # Select global covariates of subject j
    cllj <- ebmt3[j, covs]
    nd <- cbind(ttmat, rbind(cllj, cllj, cllj))
    # Make nd of class msdata to use expand.covs
    attr(nd, "trans") <- tmat
    class(nd) <- c("msdata", "data.frame")
    nd <- expand.covs(nd, covs=covs, longnames = FALSE)
    nd$drmatch.2.3 <- nd$drmatch.2 + nd$drmatch.3
    nd$pr <- 0
    nd$pr[nd$trans==3] <- 1
    nd$strata <- c(1, 2, 2)
    msfj <- msfit(c2, newdata = nd, trans = tmat)
    ptj <- probtrans(msfj, predt = 0, variance = FALSE)
    sptj <- suppressWarnings(summary(ptj, times=6))
    p12_6[j] <- sptj[[1]]$pstate2
  }
  
  
  #Using probtrans_coxph()
  
  #We need to make a data.frame containing all subjects of interest
  nd_n <- NULL
  for (j in 1:n) {
    # Select global covariates of subject j
    cllj <- ebmt3[j, covs]
    nd2 <- cbind(ttmat, rep(j, 3), rbind(cllj, cllj, cllj))
    colnames(nd2)[4] <- "id"
    # Make nd of class msdata to use expand.covs
    attr(nd2, "trans") <- tmat
    class(nd2) <- c("msdata", "data.frame")
    nd2 <- expand.covs(nd2, covs=covs, longnames = FALSE)
    nd2$drmatch.2.3 <- nd$drmatch.2 + nd$drmatch.3
    nd2$pr <- 0
    nd2$pr[nd$trans==3] <- 1
    nd2$strata <- c(1, 2, 2)
    nd_n <- rbind(nd_n, nd2)
  }
  
  icmstate_pt <- probtrans_coxph(c2, predt = 0, direction = "forward", 
                                 newdata = nd_n, trans = tmat)
  
  p12_6_icmstate <- rep(NA, n)
  time6_idx <- which.min(icmstate_pt[[1]][[1]]$time < 6) - 1
  for(i in 1:n){
    p12_6_icmstate[i] <- icmstate_pt[[i]][[1]][time6_idx, 3]
  }
  
  expect_equal(p12_6, p12_6_icmstate)
})



test_that("profvis", {
  #Code that should not run for testing.
  skip_if(TRUE)
  library(profvis)
  
  #Let us store P_{12}(6) for the first 30 subjects in ebmt3
  #Using cox model 2 & 4
  n <- 100
  
  #Using msfit and probtrans from mstate pkg
  #We need to make predictions for each person separately.
  ttmat <- to.trans2(tmat)[, c(2, 3, 1)]
  names(ttmat)[3] <- "trans"
  p12_6 <- rep(NA, n)
  for (j in 1:n) {
    # Select global covariates of subject j
    cllj <- ebmt3[j, covs]
    nd <- cbind(ttmat, rbind(cllj, cllj, cllj))
    # Make nd of class msdata to use expand.covs
    attr(nd, "trans") <- tmat
    class(nd) <- c("msdata", "data.frame")
    nd <- expand.covs(nd, covs=covs, longnames = FALSE)
    nd$drmatch.2.3 <- nd$drmatch.2 + nd$drmatch.3
    nd$pr <- 0
    nd$pr[nd$trans==3] <- 1
    nd$strata <- c(1, 2, 2)
    msfj <- msfit(c2, newdata = nd, trans = tmat)
    ptj <- probtrans(msfj, predt = 0, variance = FALSE)
    sptj <- suppressWarnings(summary(ptj, times=6))
    p12_6[j] <- sptj[[1]]$pstate2
  }
  
  
  #Using probtrans_coxph()
  
  #We need to make a data.frame containing all subjects of interest
  nd_n <- NULL
  for (j in 1:300) {
    # Select global covariates of subject j
    cllj <- ebmt3[j, covs]
    nd2 <- cbind(ttmat, rep(j, 3), rbind(cllj, cllj, cllj))
    colnames(nd2)[4] <- "id"
    # Make nd of class msdata to use expand.covs
    attr(nd2, "trans") <- tmat
    class(nd2) <- c("msdata", "data.frame")
    nd2 <- expand.covs(nd2, covs=covs, longnames = FALSE)
    nd2$drmatch.2.3 <- nd$drmatch.2 + nd$drmatch.3
    nd2$pr <- 0
    nd2$pr[nd$trans==3] <- 1
    nd2$strata <- c(1, 2, 2)
    nd_n <- rbind(nd_n, nd2)
  }
  
  profvis_pt <- profvis({probtrans_coxph(c2, predt = 0, direction = "forward", 
                                          newdata = nd_n, trans = tmat)})
  profvis_pt2 <- profvis({probtrans_coxph2(c2, predt = 0, direction = "forward", 
                                          newdata = nd_n, trans = tmat)})
  
  icmstate_pt <- probtrans_coxph(c2, predt = 0, direction = "forward", 
                                         newdata = nd_n, trans = tmat)
  icmstate_pt2 <- probtrans_coxph2(c2, predt = 0, direction = "forward", 
                                           newdata = nd_n, trans = tmat)
  
})




