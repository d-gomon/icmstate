test_that("Correct hazard when generating", {
  skip_if(TRUE)
  eval_times <- function(n_obs, stop_time){
    cumsum( c( 0,  runif( n_obs-1, 0, 2*(stop_time-4)/(n_obs-1) ) ) )
  }
  
  #########################---Weibull stable initial conditions---###############
  asd <- sim_id_weib_exact(n = 10000, n_obs = 6, stop_time = 15, eval_times = eval_times, shape = c(0.5, 0.5, 2), scale = c(5, 10, 10/gamma(1.5)))
  qwe <- reshape(asd, direction = "wide", idvar = "id", timevar = "state")
  qwe$status <- rep(1, nrow(qwe))
  

  #Need to set transition times to 3 to censored so we can estimate correct cumulative hazard
  #From 1 to 2:
  dat1 <- qwe
  dat1[is.na(dat1[, "time.2"]), "status"] <- 0
  dat1[is.na(dat1[, "time.2"]), "time.2"] <- dat1[is.na(dat1[, "time.2"]), "time.3"]
  first_trans <- survival::survfit(Surv(time.2, status) ~ 1, data = dat1)
  
  #From 1 to 3
  dat2 <- qwe
  dat2[!is.na(dat2[, "time.2"]), "status"] <- 0
  dat2[!is.na(dat2[, "time.2"]), "time.3"] <- dat2[!is.na(dat2[, "time.2"]), "time.2"]
  second_trans <- survival::survfit(Surv(time.3, status) ~ 1, data = dat2)
  
  #From 2 to 3:
  third_trans <- survival::survfit(Surv(time.2, time.3, status) ~ 1, data = subset(qwe, !is.na(time.2)))
  
  
  #Make plots:
  #Transition 1
  plot(first_trans$time, first_trans$cumhaz, type = "s", main = "first transition", xlim = c(0, 15))
  lines(first_trans$time, -pweibull(first_trans$time, 0.5, 5, lower = FALSE, log = TRUE), col = "red")
  #Aligns well
  
  #Transition 2
  plot(second_trans$time, second_trans$cumhaz, type = "s", main = "second transition", xlim = c(0, 15))
  lines(second_trans$time, -pweibull(second_trans$time, 0.5, 10, lower = FALSE, log = TRUE), col = "red")
  #Aligns well
  
  #Transition 3
  plot(third_trans$time, third_trans$cumhaz, type = "s", main = "third transition", xlim = c(0, 15), ylim = c(0, 5))
  lines(third_trans$time, -pweibull(third_trans$time, 2, 10/gamma(1.5), lower = FALSE, log = TRUE), col = "red")
  #Aligns well
  
  
  #########################---Weibull equalprob initial conditions---###############
  asd <- sim_id_weib_exact(n = 10000, n_obs = 6, stop_time = 15, eval_times = eval_times,
                           start_state = "equalprob", shape = c(0.5, 0.5, 2), scale = c(5, 10, 10/gamma(1.5)))
  qwe <- reshape(asd, direction = "wide", idvar = "id", timevar = "state")
  qwe$status <- rep(1, nrow(qwe))
  
  #Need to set transition times to 3 to censored so we can estimate correct cumulative hazard
  #From 1 to 2:
  dat1 <- qwe
  dat1[is.na(dat1[, "time.2"]), "status"] <- 0
  dat1[is.na(dat1[, "time.2"]), "time.2"] <- dat1[is.na(dat1[, "time.2"]), "time.3"]
  dat1 <- dat1[!dat1[, "time.2"] == 0,]
  first_trans <- survival::survfit(Surv(time.2, status) ~ 1, data = dat1)
  
  #From 1 to 3
  dat2 <- qwe
  dat2[!is.na(dat2[, "time.2"]), "status"] <- 0
  dat2[!is.na(dat2[, "time.2"]), "time.3"] <- dat2[!is.na(dat2[, "time.2"]), "time.2"]
  second_trans <- survival::survfit(Surv(time.3, status) ~ 1, data = dat2)
  
  #From 2 to 3:
  third_trans <- survival::survfit(Surv(time.2, time.3, status) ~ 1, data = subset(qwe, !is.na(time.2)))
  
  
  #Make plots:
  #Transition 1
  plot(first_trans$time, first_trans$cumhaz, type = "s", main = "first transition", xlim = c(0, 15))
  lines(first_trans$time, -pweibull(first_trans$time, 0.5, 5, lower = FALSE, log = TRUE), col = "red")
  #Aligns well
  
  #Transition 2
  plot(second_trans$time, second_trans$cumhaz, type = "s", main = "second transition", xlim = c(0, 15))
  lines(second_trans$time, -pweibull(second_trans$time, 0.5, 10, lower = FALSE, log = TRUE), col = "red")
  #Aligns well
  
  #Transition 3
  plot(third_trans$time, third_trans$cumhaz, type = "s", main = "third transition", xlim = c(0, 15), ylim = c(0, 5))
  lines(third_trans$time, -pweibull(third_trans$time, 2, 10/gamma(1.5), lower = FALSE, log = TRUE), col = "red")
  #Aligns well
  
  
  
  #########################---Exponential equalprob initial conditions---###############
  asd <- sim_id_weib_exact(n = 10000, n_obs = 6, stop_time = 15, eval_times = eval_times,
                           start_state = "equalprob", shape = c(1, 1, 1), scale = c(10, 20, 10))
  qwe <- reshape(asd, direction = "wide", idvar = "id", timevar = "state")
  qwe$status <- rep(1, nrow(qwe))
  

  #Need to set transition times to 3 to censored so we can estimate correct cumulative hazard
  #From 1 to 2:
  dat1 <- qwe
  dat1[is.na(dat1[, "time.2"]), "status"] <- 0
  dat1[is.na(dat1[, "time.2"]), "time.2"] <- dat1[is.na(dat1[, "time.2"]), "time.3"]
  dat1 <- dat1[!dat1[, "time.2"] == 0,]
  first_trans <- survival::survfit(Surv(time.2, status) ~ 1, data = dat1)
  
  #From 1 to 3
  dat2 <- qwe
  dat2[!is.na(dat2[, "time.2"]), "status"] <- 0
  dat2[!is.na(dat2[, "time.2"]), "time.3"] <- dat2[!is.na(dat2[, "time.2"]), "time.2"]
  second_trans <- survival::survfit(Surv(time.3, status) ~ 1, data = dat2)
  
  #From 2 to 3:
  third_trans <- survival::survfit(Surv(time.2, time.3, status) ~ 1, data = subset(qwe, !is.na(time.2)))
  
  
  #Make plots:
  #Transition 1
  plot(first_trans$time, first_trans$cumhaz, type = "s", main = "first transition", xlim = c(0, 15))
  lines(first_trans$time, -pexp(first_trans$time, 0.1,  lower = FALSE, log = TRUE), col = "red")
  #Aligns well
  
  #Transition 2
  plot(second_trans$time, second_trans$cumhaz, type = "s", main = "second transition", xlim = c(0, 15))
  lines(second_trans$time, -pexp(second_trans$time, 0.05, lower = FALSE, log = TRUE), col = "red")
  #Aligns well
  
  #Transition 3
  plot(third_trans$time, third_trans$cumhaz, type = "s", main = "third transition", xlim = c(0, 15), ylim = c(0, 5))
  lines(third_trans$time, -pexp(third_trans$time, 0.1, lower = FALSE, log = TRUE), col = "red")
  #Aligns well
  
  
})
