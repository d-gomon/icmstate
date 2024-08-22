#' Plot the transition specific survival probabilities for a fitted \code{npmsm} model
#' 
#' @description For a fitted \code{npmsm} model plot the transition specific 
#' survival probabilities. These are given by the product integral of the hazard 
#' increments estimated for a single transition. This is equivalent to a Kaplan-Meier 
#' estimator ignoring the existence of all other transitions.
#' 
#' 
#' @param npmsmlist An "npmsm" object or a list containing multiple "npmsm" objects
#' @param landmark A landmark time indicating from which time on survival should be determined.
#' If missing, the smallest time in the first "npmsm" object will be used.
#' @param support Should the support regions be displayed as rectangles?
#' @param sup_cutoff Cutoff to be used for determining the support intervals.
#' 
#' @importFrom mstate to.trans2 msfit probtrans
#' @import checkmate
#' @export
#' 
#' 
#' 
#' 



plot_surv <- function(npmsmlist, landmark, support = FALSE, sup_cutoff = 1e-8){
  
  
  # Remove CRAN check notes
  trans <- time <- surv <- obj_name <- NULL
  
  #----------Data checks----------
  #TODOTODOTODO
  #Need to ensure that each tmat of each element in the list is the same.
  #Need to check that landmark is numeric or missing. 
  #Also, landmark should be smaller than largest time in npmsm$A$Haz$time
  
  #------Data-preprocessing------------------
  if(inherits(npmsmlist, "npmsm")){
    dat <- list(npmsmlist)  
  } else {
    dat <- npmsmlist
  }
  
  if(missing(landmark)){
    landmark = min(dat[[1]]$A$Haz$time, 0)
  }
  
  #------Extract summaries-----------
  tmat2 <- to.trans2(dat[[1]]$tmat)
  M <- nrow(tmat2)

  
  
  #------Function to calculate survival-----------
  prodint_surv <- function(haz){ #Returns product integral (survival estimate) of intensities
    cumprod((1-haz))
  }
  
  #---------Fill plotting data frame-------------
  plot_df <- data.frame(time = numeric(0), surv = numeric(0), trans = numeric(0), obj_name = factor())
  for(i in 1:length(dat)){
    for(m in 1:M){
      Haz <- subset(dat[[i]]$A$Haz, trans == m) #Only this transition
      Haz$haz <- diff(c(0, Haz$Haz)) #Calculate intensity from cumulative intensity
      Haz <- subset(Haz, time > landmark) #Retain only intensities after landmark time
      Haz <- rbind(c(landmark, 0, m, 0), Haz)
      plot_df <- rbind(plot_df, data.frame(time = Haz$time, 
                                           surv = prodint_surv(Haz$haz), 
                                           trans = rep(m, length(Haz$time)), 
                                           obj_name = as.factor(rep(paste0("df", i), length(Haz$time)) )))
    }
  }
  trans_names <- tmat2$transname
  names(trans_names) <- 1:M
  
  
  g_plot <- ggplot(data = plot_df) + 
    aes(x = time, y = surv, col = obj_name, group = obj_name, fill = obj_name) + 
    geom_step(aes(linetype = obj_name), lwd = 1.3, alpha = 0.8) + 
    facet_grid(trans~., labeller = as_labeller(trans_names)) + 
    geom_vline(xintercept = landmark, linetype = "longdash") +
    xlim(min(dat[[1]]$A$Haz$time, 0), NA) +
    ylim(0, 1) +
    labs(colour='Model') 
  
  #---------Create support data_frame-----------
  supp_df <- data.frame(L = numeric(0), R = numeric(0), trans = numeric(0), obj_name = factor())
  if(support){
    for(i in 1:length(dat)){
      suppi <- support_npmsm(dat[[i]], cutoff = sup_cutoff)
      for(m in 1:M){
        L <- suppi[[m]]$support[, 1]
        R <- suppi[[m]]$support[, 2]
        L <- L[which(R > landmark)]
        R <- R[which(R > landmark)]
        landmark_within_support <- which(L < landmark & R > landmark)
        if(length(landmark_within_support) > 0){
          warning(paste0("Landmark time chosen within support interval for transition: ", trans_names[m], ". Survival estimates may not be accurate. 
                         Choose landmark time smaller than ", L[landmark_within_support], " or larger than ", R[landmark_within_support]))
        }
        supp_df <- rbind(supp_df, data.frame(L = L,
                                             R = R,
                                             trans = m,
                                             obj_name = as.factor( rep(paste0("df", i), length(L) ) )))
      }
    }
    g_plot <- g_plot + geom_rect(data = supp_df, inherit.aes = FALSE, mapping = aes(xmin = L, xmax = R, ymin = 0, ymax = 1, fill = obj_name), alpha = 0.3)
  }
  
  
  
  
  print(g_plot)
  out <- list(plot_df = plot_df,
              g_plot = g_plot)
  return(invisible(out))
}