#' Plot the transition probabilities for a fitted \code{npmsm} model
#' 
#' @description For a fitted \code{npmsm} model plot the transition probabilities
#' from a certain state for all possible (direct and indirect) transitions.
#' 
#' 
#' @param npmsmlist An "npmsm" object or a list containing multiple "npmsm" objects
#' @param from A numeric value indicating the state from which we consider the 
#' transition probabilities. Default is NULL, meaning we consider transition probabilities 
#' from all states from which a direct transition is possible.
#' @param to A numeric vector indicating to which states we consider the transition 
#' probabilities. Only states that can be reached from the \code{"from"} state are considered.
#' @param transitions A numeric vector indicating which transitions to consider (plot). 
#' Can only be used if \code{"from"} is not specified, as it only works for direct transitions.
#' @param landmark A landmark time indicating from which time on survival should be determined.
#' If missing, the smallest between the time in the first "npmsm" object or 0 will be used.
#' @param interpolate Should the cumulative hazard be linearly interpolated before 
#' determining transition probabilities? Default is TRUE.
#' @param facet Should the resulting plot be faceted (one panel per transition)? 
#' Default is TRUE.
#' @param times_interpol At which times should the cumulative hazard be interpolated?
#' Only necessary to specify if interpolate = TRUE.
#' @param c.legend Should legend be displayed for colour (different entries in 
#' \code{npmsmlist})? Default is TRUE.
#' @param c.names A character vector indicating the names to display in the legend.
#' These names should represent the entries in \code{npmsmlist.}. Default = \code{NULL}.
#' 
#' @importFrom mstate to.trans2 msfit probtrans
#' @import checkmate
#' @export
#' 
#' 
#' 
#' @returns A plot will be produced in the plotting window. When assigning  
#' the output to an object, the underlying data frame used for plotting 
#' and a \code{'ggplot'} object will be returned in a list.
#' 
#' @examples 
#' require(mstate)
#' require(ggplot2)
#' #Generate from an illness-death model with exponential transitions with 
#' #rates 1/2, 1/10 and 1 for 10 subjects over a time grid.
#' gd <- sim_weibmsm(tmat = trans.illdeath(), shape = c(1,1,1),
#'                   scale = c(2, 10, 1), n_subj = 10, obs_pars = c(2, 0.5, 20), 
#'                   startprobs = c(0.9, 0.1, 0))
#' #Fit 2 models: 1 with at most 4 iterations and 1 with at most 20
#' mod1 <- npmsm(gd, trans.illdeath(), maxit = 4)
#' mod2 <- npmsm(gd, trans.illdeath(), maxit = 20)
#' 
#' #Plot the transition probabilities from state 1, without interpolating 
#' #the cumulative hazard for the npmsm runs with max 4 and 20 iterations.
#' plot_probtrans(list(mod1, mod2), from = 1, interpolate = FALSE,
#'                c.names = c("4 iterations", "20 iterations"))
#' 
#' 



plot_probtrans <- function(npmsmlist, from = NULL, to = NULL, transitions = NULL,
                           landmark, interpolate = TRUE, facet = TRUE,
                           times_interpol = NULL, c.legend = TRUE, c.names = NULL){
  # Remove CRAN check notes
  trans <- time <- surv <- obj_name <- NULL
  
  #----------Data checks----------
  #TODOTODOTODO
  #Need to ensure that each tmat of each element in the list is the same.
  #Need to check that landmark is numeric or missing. 
  #Also, landmark should be smaller than largest time in npmsm$A$Haz$time
  
  #------Data-preprocessing------------------
  if(inherits(npmsmlist, "npmsm") | inherits(npmsmlist, "msm")){
    dat <- list(npmsmlist)  
  } else {
    dat <- npmsmlist
  }
  
  if(inherits(dat[[1]], "msm") & missing(times_interpol)){
    stop("Please specify times_interpol for list of class 'msm' to determine transition probabilities.")
  }
  
  if(missing(landmark)){
    if(inherits(dat[[1]], "npmsm")){
      landmark = min(dat[[1]]$A$Haz$time, 0)  
    } else{
      landmark = min(times_interpol)
    }
  }
  
  if(isTRUE(interpolate)){
    if(is.null(times_interpol)){
      stop("Please specify at which times interpolation should occur, or set interpolate to FALSE.")
    }
  }
  
  if(is.null(c.names)){
    c.names <- paste0("dat", 1:length(npmsmlist))
  }
  
  

  
  #------Extract summaries-----------
  if(inherits(dat[[1]], "npmsm")){
    tmat <- dat[[1]]$tmat
    tmat2 <- to.trans2(tmat)
    M <- nrow(tmat2)
    n_states <- nrow(tmat)  
  } else{
    tmat <- tmat_from_msm(dat[[1]])
    tmat2 <- to.trans2(tmat)
    M <- nrow(tmat2)
    n_states <- dat[[1]]$qmodel$nstates
  }
  
  
  if(!is.null(from)){ 
    #We want to extract the transition probabilities only for transitions from state "from"
    trans_names <- paste0(colnames(tmat)[from], "->", colnames(tmat))
    names(trans_names) <- 1:n_states
    remove_trans <- NULL
    #Determine which states we can transition to
    if(!is.null(to)){
      to <- intersect(to, 1:n_states)
    } else{
      to <- 1:n_states
    }
  } else{ #Consider all direct transitions only
    if(is.null(transitions)){
      transitions <- tmat2[, "transno"]
    }
    which_transitions <- intersect(tmat2[, "transno"], transitions)
    trans_names <- tmat2[which_transitions, "transname"]
    names(trans_names) <- which_transitions
  }
  
  
  
  
  
  #---------Fill plotting data frame-------------
  plot_df <- data.frame(time = numeric(0), surv = numeric(0), trans = numeric(0), obj_name = factor())
  for(i in 1:length(dat)){
    
    
    ###########################NPMSM PART#############################
    if(inherits(dat[[1]], "npmsm")){
      #Determine interpolated transition probabilities
      if(isTRUE(interpolate)){
        #Can be quite slow for many interpolation time points.
        probtrans_orig <- suppressWarnings({
          probtrans(interpol_msfit(dat[[i]]$A, times = times_interpol), predt = landmark, direction = "forward", variance = FALSE)
        }) 
      } else{
        probtrans_orig <- suppressWarnings({
          probtrans(dat[[i]]$A, predt = landmark, direction = "forward", variance = FALSE)
        })
      }
      if(is.null(from)){ #If we did not specify a state from which we are interested in transitions
        for(m in seq_along(which_transitions)){
          from_state <- tmat2[tmat2$transno == which_transitions[m], "from"]
          which_states_to <- tmat2[tmat2$transno == which_transitions[m], "to"]
          states_to_transno <- tmat2[tmat2$transno == which_transitions[m], "transno"]
          for(j in seq_along(which_states_to)){
            plot_df <- rbind(plot_df, data.frame(time = probtrans_orig[[from_state]][, 1], 
                                                 surv = probtrans_orig[[from_state]][, which_states_to[j]+1],
                                                 trans = rep(states_to_transno[j], length(probtrans_orig[[from_state]][, 1])),
                                                 obj_name = as.factor(rep(c.names[i], length(probtrans_orig[[from_state]][, 1]))) ))
          }
        }
      } else{ #If we are interested in all possibly indirect transitions from a state
        probtrans_orig <- probtrans_orig[[from]]
        for(j in to+1){
          if(all(probtrans_orig[, j] == 0)){
            remove_trans <- c(remove_trans, j-1)
            next
          }
          plot_df <- rbind(plot_df, data.frame(time = probtrans_orig[, 1], 
                                               surv = probtrans_orig[, j],
                                               trans = rep(j-1, length(probtrans_orig[, 1])),
                                               obj_name = as.factor(rep(c.names[i], length(probtrans_orig[, 1]))) ))
        }
      }
    } else if(inherits(dat[[1]], "msm")){
      ###################################MSM PART########################
      probtrans_orig <- suppressWarnings({
        transprob(dat[[i]], times = times_interpol)
      }) 
      
      if(is.null(from)){ #If we did not specify a state from which we are interested in transitions
        for(m in seq_along(which_transitions)){
          from_state <- tmat2[tmat2$transno == which_transitions[m], "from"]
          which_states_to <- tmat2[tmat2$transno == which_transitions[m], "to"]
          states_to_transno <- tmat2[tmat2$transno == which_transitions[m], "transno"]
          for(j in seq_along(which_states_to)){
            plot_df <- rbind(plot_df, data.frame(time = probtrans_orig[[from_state]][, 1], 
                                                 surv = probtrans_orig[[from_state]][, which_states_to[j]+1],
                                                 trans = rep(states_to_transno[j], length(probtrans_orig[[from_state]][, 1])),
                                                 obj_name = as.factor(rep(c.names[i], length(probtrans_orig[[from_state]][, 1]))) ))
          }
        }
      } else{ #If we are interested in all possibly indirect transitions from a state
        probtrans_orig <- probtrans_orig[[from]]
        for(j in to+1){
          if(all(probtrans_orig[, j] == 0)){
            remove_trans <- c(remove_trans, j-1)
            next
          }
          plot_df <- rbind(plot_df, data.frame(time = probtrans_orig[, 1], 
                                               surv = probtrans_orig[, j],
                                               trans = rep(j-1, length(probtrans_orig[, 1])),
                                               obj_name = as.factor(rep(c.names[i], length(probtrans_orig[, 1]))) ))
        }
      }
    }  #End of MSM part
    
    
    
  }

  
  if(isTRUE(facet)){
    g_plot <- ggplot(data = plot_df) + 
      aes(x = time, y = surv, col = obj_name, group = obj_name, fill = obj_name) + 
      geom_step(linetype = "longdash", lwd = 1.3, alpha = 0.8) + 
      facet_grid(trans~., labeller = as_labeller(trans_names)) + 
      geom_vline(xintercept = landmark, linetype = "longdash") +
      xlim(min(dat[[1]]$A$Haz$time, 0), NA) +
      ylim(0, 1) +
      xlab("Time") + 
      ylab("Probability") + 
      labs(colour='Model') +
      ggtitle("Transition probabilities") + 
      scale_color_brewer(palette="Dark2")
  } else{
    g_plot <- ggplot(data = plot_df) + 
      aes(x = time, y = surv, colour = obj_name, group = interaction(obj_name, trans), linetype = as.factor(trans)) + 
      geom_step(lwd = 1.3, alpha = 0.8) + 
      geom_vline(xintercept = landmark, linetype = "longdash") +
      xlim(min(dat[[1]]$A$Haz$time, 0), NA) +
      ylim(0, 1) +
      xlab("Time") + 
      ylab("Probability") + 
      labs(colour='Model', linetype = 'Transition') +
      ggtitle("Transition probabilities") +
      scale_color_brewer(palette="Dark2") + 
      scale_linetype(labels = trans_names)
  }
  
  if(isFALSE(c.legend)){
    g_plot <- g_plot + guides(colour = "none")
  }
  
  
  print(g_plot)
  out <- list(plot_df = plot_df,
              g_plot = g_plot)
  return(invisible(out))
}