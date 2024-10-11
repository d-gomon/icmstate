#' Visualise multi-state data
#' 
#' @description
#' Produce a plot with the y-axis representing subjects in the data and the 
#' x-axis representing the time at which states have been observed.
#' 
#' @param gd A \code{data.frame} containing the following named columns:
#'\describe{
#'   \item{\code{id}:}{Identifier of subject;}
#'   \item{\code{state}:}{state of subject at \code{time};}
#'   \item{\code{time}:}{time at which subject is observed;}
#' }
#' @param npmsm Output from \code{\link{npmsm}} function
#' @param tmat A transition matrix as created by \code{transMat}
#' @param neat Boolean indicating whether redundant observations should be 
#' removed in the plot. Default is TRUE
#' @param cutoff cutoff value for numerically determining the support using
#' \code{\link{support_npmsm}}
#' 
#' @return A plot will be produced in the plotting window.
#' 
#' @import checkmate
#' @import ggplot2
#' @importFrom prodlim row.match
#' 
#' @export
#' 
#' 
#' @examples 
#' #Write a function for evaluation times: observe at 0 and uniform inter-observation times.
#' eval_times <- function(n_obs, stop_time){
#'   cumsum( c( runif(1, 0, 0.5),  runif( n_obs-1, 0, 2*(stop_time-4)/(n_obs-1) ) ) )
#' }
#' #Use built_in function to simulate illness-death data
#' #from Weibull distributions for each transition
#' sim_dat <- sim_id_weib(n = 50, n_obs = 6, stop_time = 15, eval_times = eval_times,
#'                       start_state = "stable", shape = c(0.5, 0.5, 2), 
#'                       scale = c(5, 10, 10/gamma(1.5)))
#'
#' #Visualise the data
#' visualise_msm(sim_dat)                       
#' 

visualise_msm <- function(gd, npmsm, tmat, neat = TRUE, cutoff){
  
  # Remove CRAN notes
  id <- ID <- t1 <- t2 <- color <- state1 <- state2 <- NULL

# Argument checks ---------------------------------------------------------

  arg_checks <- makeAssertCollection()
  assertMultiClass(gd, c("matrix", "data.frame"), add = arg_checks)
  if(inherits(gd, "data.frame")){
    assertDataFrame(gd, min.cols = 3, max.cols = 3, add = arg_checks)
    
  } else if(inherits(gd, "matrix")){
    assertMatrix(gd, min.cols = 3, max.cols = 3, add = arg_checks)
  }
  assertNames(colnames(gd), must.include = c("id", "state", "time"), add = arg_checks)
  assertNumeric(gd[, "time"], lower = 0, any.missing = FALSE, add = arg_checks)
  assertIntegerish(gd[, "state"])
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  

# Data Transformations ----------------------------------------------------
  if(isTRUE(neat)){
    if(!missing(npmsm)){
      gd <- remove_redundant_observations(gd = gd, tmat = npmsm$tmat) 
    } else if(!missing(tmat)){
      gd <- remove_redundant_observations(gd = gd, tmat = tmat)
    } else{
      gd <- remove_redundant_observations(gd = gd)
    }
  }
  
  
  #Order gd according to ID and then time
  gd <- gd[order(gd[, "id"], gd[, "time"]),]
  

# Create Visualisation Data Frame -----------------------------------------

  plot_frame <- NULL
  for(i in unique(gd[, "id"])){
    gdi <- gd[gd[, "id"] == i,]
    gdi <- gdi[order(gdi[, "time"]),]
    n <- nrow(gdi)
    ploti <- matrix(NA, nrow = 0, ncol = 3)
    diffi <- which(diff(gdi[, "state"]) != 0)
    #Determine minimal observation time & state
    mintime <- min(gdi[, "time"])
    minstate <- gdi[1, "state"]
    #Display transitions:
    if(length(diffi) != 0){
      for(j in diffi){
        plot_frame <- rbind(plot_frame, c(i, gdi[j, "time"], gdi[j+1, "time"], gdi[j, "state"], gdi[j+1, "state"], mintime, minstate))
      }  
    } else{
      plot_frame <- rbind(plot_frame, c(i, gdi[nrow(gdi), "time"], gdi[nrow(gdi), "time"], gdi[nrow(gdi), "state"], gdi[nrow(gdi), "state"], mintime, minstate))
    }
    if(length(diffi) == 0){
      diffi <- -Inf
    }
    if((max(diffi)+1) != n){
      plot_frame <- rbind(plot_frame, c(i, gdi[n, "time"], gdi[n, "time"], gdi[n, "state"], gdi[n, "state"], mintime, minstate))
    }
  }
  plot_frame <- as.data.frame(plot_frame)
  plot_frame[, 1] <- as.factor(plot_frame[, 1])
  colnames(plot_frame) <- c("ID", "t1", "t2", "state1", "state2", "mintime", "minstate")
  plot_frame$color <- as.factor(row.match(plot_frame[, c("state1", "state2")], unique(plot_frame[, c("state1", "state2")])))


# Support Plotting (if possible) ------------------------------------------

  if(!missing(npmsm)){
    tmat <- npmsm$tmat
    if(!missing(cutoff)){
      support <- support_npmsm(npmsm, cutoff)  
    } else{
      support <- support_npmsm(npmsm)
    }
    for(i in 1:length(support)){
      plot_framei <- as.data.frame(support[[i]]$support)
      if(nrow(plot_framei) == 0){
        next
      }
      colnames(plot_framei) <- c("t1", "t2")
      plot_framei$ID <- names(support)[i]
      plot_framei$state1 <- support[[i]]$transition[1]
      plot_framei$state2 <- support[[i]]$transition[2]
      plot_framei$color <- "dummy"
      plot_framei$mintime <- 0
      plot_framei$minstate <- "|"
      plot_framei <- plot_framei[, c("ID", "t1", "t2", "state1", "state2", "mintime", "minstate", "color")]
      plot_frame <- rbind(plot_frame, plot_framei)
    }
    plot_frame$color <- as.factor(row.match(plot_frame[, c("state1", "state2")], unique(plot_frame[, c("state1", "state2")])))
  }
  
  

# Create ggplot2 plot -----------------------------------------------------

  out <- ggplot(plot_frame, aes(x = ID)) +
    geom_linerange(aes(ymin = mintime, ymax = t1)) +
    #First segment
    geom_linerange(aes(ymin = t1, ymax = t2, lwd = 3, col = color), alpha = .4) +
    geom_text(aes(x = ID, y = mintime, label = minstate)) +
    #Add state text
    geom_text(aes(x = ID, y = t1, label = state1)) +
    geom_text(aes(y = t2, label = state2)) +
    #flip coordinates
    coord_flip() +
    #theme
    theme_bw() +
    theme(legend.position = "none") +
    ylab("Time") +
    xlab("Subject")
  ret <- list(plot = out)
  if(!missing(npmsm)){
    ret$support <- support  
  }
  print(out)
  return(invisible(ret))
}