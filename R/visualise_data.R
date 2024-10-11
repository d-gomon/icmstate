#' Visualise data for illness-death model, only applicable to Frydman(1995) setting.
#' 
#' @param msmFrydman from \code{\link{msm_frydman}}
#' @inheritParams msm_frydman
#' 
#' 
#' @import checkmate
#' @import ggplot2
#' 
#' @references Frydman, H. (1995). Nonparametric Estimation of a Markov 
#' 'Illness-Death' Process from Interval- Censored Observations, with 
#' Application to Diabetes Survival Data. Biometrika, 82(4), 773-789. 
#' \doi{10.2307/2337344}
#' 
#' @export


visualise_data <- function(data, msmFrydman){
  # Remove CRAN notes
  delta <- Delta <- ID <- group <- t1 <- t2 <- t3 <- time <- censored <- NULL
  
  
  arg_checks <- makeAssertCollection()
  assertDataFrame(data, min.cols = 5, max.cols = 6, add = arg_checks)
  assertNames(names(data), must.include = c("delta", "Delta", "L", "R", "time"), add = arg_checks)
  assertSubset(data[["delta"]], c(0,1,2), add = arg_checks)
  assertSubset(data[["Delta"]], c(0,1), add = arg_checks)
  assertNumeric(data[["L"]], lower = 0, add = arg_checks)
  assertNumeric(data[["R"]], lower = 0, add = arg_checks)
  if(!missing(msmFrydman)){
    assertClass(msmFrydman, "msmFrydman", add = arg_checks)  
  }
  assertNumeric(data[["time"]], lower = 0, any.missing = FALSE, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  
  #Assign group labels to data
  data$group <- factor(ifelse(data[["delta"]] == 0, ifelse(data[["Delta"]] == 0, 1L, 2L), 
                       ifelse(data[["Delta"]] == 0, 3L, 4L)))
  data <- subset(data, select = -c(delta, Delta))
  #Order data by group
  data <- data[order(data$group, data$time),]
  data$ID <- factor(rownames(data), levels = rownames(data))
  data$t1 <- ifelse(data$group %in% c(1,2), data$time, data$L)
  data$t2 <- ifelse(data$group %in% c(1,2), NA, data$R)
  data$t3 <- ifelse(data$group %in% c(1,2), NA, data$time)
  data$censored <- factor(ifelse(data$group %in% c(1,3), 0, 1))
  
  if(!missing(msmFrydman)){
    levels(data$group) <- c(levels(data$group), 5)
    supportdf <- as.data.frame(msmFrydman$supportMSM$Q_mat)
    supportdf$ID <- rep("support", nrow(supportdf))
    supportdf$group <- rep(5, nrow(supportdf))
    supportdf$t1 <- supportdf$L
    supportdf$t2 <- supportdf$R
    supportdf$t3 <- rep(NA, nrow(supportdf))
    supportdf$time <- rep(NA, nrow(supportdf))
    supportdf$censored <- rep(NA, nrow(supportdf))
    data <- rbind(data, supportdf)
  }
  
  
  out <- ggplot(data, aes(x = ID, colour = group)) +
    #First segment
    geom_linerange(aes(ymin = 0, ymax = t1)) +
    #Second segment
    geom_linerange(aes(ymin = t1, ymax = t2), lwd = 3, na.rm = TRUE) +
    #Third segment
    geom_linerange(aes(ymin = t2, ymax = t3), lwd = 1, na.rm = TRUE) +
    #Censoring/event indicator
    geom_point(data = data[!is.na(data$censored),], aes(y = time, shape = censored), size = 4) +
    #Scale shape for censoring/event
    scale_shape_manual(name = "Event", values = c("0" = 1, "1" = 4), 
                       labels = c("censored", "observed"), na.value = NA) +
    #Scale colour
    scale_colour_manual(name = "Group", values = c("1" = "green", "2" = "lightblue", "3" = "red", "4" = "purple", "5" = "orange"),
                        labels = c("RC", "Obs", "Interval + RC", "Interval + Obs", "support")) +
    #flip coordinates
    coord_flip() +
    #theme
    theme_bw()

  out
}




