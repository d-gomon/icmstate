#' Visualise data 
#' 
#' @param gd A \code{data.frame} containing the following named columns:
#'\describe{
#'   \item{\code{id}:}{Identifier of subject;}
#'   \item{\code{state}:}{state of subject at \code{time};}
#'   \item{\code{time}:}{time at which subject is observed;}
#' }
#' @param npmsm Output from \code{\link{npmsm}} function
#' @param neat Boolean indicating whether redundant observations should be 
#' removed in the plot. Default is TRUE
#' @param cutoff cutoff value for numerically determining the support using
#' \code{\link{support_npmsm}}
#' @param ... Further arguments to plot
#' 
#' @import checkmate
#' @import ggplot2
#' @importFrom prodlim row.match
#' 
#' @export

visualise_msm <- function(gd, npmsm, neat = TRUE, cutoff, ...){
  
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
    } else{
      gd <- remove_redundant_observations(gd = gd)
    }
  }
  

# Create Visualisation Data Frame -----------------------------------------

  plot_frame <- NULL
  for(i in unique(gd[, "id"])){
    gdi <- gd[gd[, "id"] == i,]
    gdi <- gdi[order(gdi[, "time"]),]
    n <- nrow(gdi)
    ploti <- matrix(NA, nrow = 0, ncol = 3)
    diffi <- which(diff(gdi[, "state"]) != 0)
    if(length(diffi) != 0){
      for(j in diffi){
        plot_frame <- rbind(plot_frame, c(i, gdi[j, "time"], gdi[j+1, "time"], gdi[j, "state"], gdi[j+1, "state"]))
      }  
    } else{
      plot_frame <- rbind(plot_frame, c(i, gdi[nrow(gdi), "time"], gdi[nrow(gdi), "time"], gdi[nrow(gdi), "state"], gdi[nrow(gdi), "state"]))
    }
    if(length(diffi) == 0){
      diffi <- -Inf
    }
    if((max(diffi)+1) != n){
      plot_frame <- rbind(plot_frame, c(i, gdi[n, "time"], gdi[n, "time"], gdi[n, "state"], gdi[n, "state"]))
    }
  }
  plot_frame <- as.data.frame(plot_frame)
  plot_frame[, 1] <- as.factor(plot_frame[, 1])
  colnames(plot_frame) <- c("ID", "t1", "t2", "state1", "state2")
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
      plot_framei <- plot_framei[, c("ID", "t1", "t2", "state1", "state2", "color")]
      plot_frame <- rbind(plot_frame, plot_framei)
    }
    plot_frame$color <- as.factor(row.match(plot_frame[, c("state1", "state2")], unique(plot_frame[, c("state1", "state2")])))
  }
  
  

# Create ggplot2 plot -----------------------------------------------------

  out <- ggplot(plot_frame, aes(x = ID)) +
    geom_linerange(aes(ymin = 0, ymax = t1)) +
    #First segment
    geom_linerange(aes(ymin = t1, ymax = t2, lwd = 3, col = color), alpha = .4) +
    #Add state text
    geom_text(aes(x = ID, y = t1, label = state1)) +
    geom_text(aes(y = t2, label = state2)) +
    #flip coordinates
    coord_flip() +
    #theme
    theme_bw() +
    theme(legend.position = "none")
  ret <- list(plot = out)
  if(!missing(npmsm)){
    ret$support <- support  
  }
  print(out)
  return(invisible(ret))
}