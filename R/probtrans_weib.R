#' 
#' Determine transition probabilities for a multi-state model with Weibull
#' hazards for the transitions.
#' 
#' 
#' @param transMat A transition matrix as generated by \code{mstate::transMat}
#' describing the possible transitions for the multi-state model.
#' @param times The times at which the transition probabilities should be 
#' determined. Will always determine the probabilities forward in time starting 
#' from \code{min(times)}.
#' @param shapes The Weibull shapes corresponding to the numbered transitions 
#' in \code{transMat}. See \code{?pweibull} for more info.
#' @param scales The Weibull scales corresponding to the numbered transitions 
#' in \code{transMat}. See \code{?pweibull} for more info.
#' @param type Should the transition probabilities be determined using 
#' product integration \code{"prodint"} or by solving the Kolmogorov forward ordinary differential 
#' equation \code{"ODE"}.
#' 
#' 
#' @return An object containing the "true" transition probabilities for the 
#' specified Weibull hazards.
#' 
#' @export
#' @importFrom mstate probtrans
#' @importFrom deSolve ode
#' @importFrom stats pweibull
#' 
#' @examples 
#' #Illness-death model
#' tmat <- mstate::trans.illdeath()
#' IDM <- probtrans_weib(tmat, seq(0, 15, 0.01), shapes = c(0.5, 0.5, 2), 
#'                       scales = c(5, 10, 10/gamma(1.5)), type = "prodint")
#' IDM2 <- probtrans_weib(tmat, seq(0, 15, 0.01), shapes = c(0.5, 0.5, 2), 
#'                        scales = c(5, 10, 10/gamma(1.5)), type = "ODE")
#' plot(IDM)
#' plot(IDM2)
#' 
#' #Extended illness-death model
#' tmat <- mstate::transMat(list(c(2, 3), c(4), c(), c()))
#' IDM <- probtrans_weib(tmat, seq(0, 15, 0.01), shapes = c(0.5, 0.5, 2), 
#'                       scales = c(5, 10, 10/gamma(1.5)), type = "prodint")
#' IDM2 <- probtrans_weib(tmat, seq(0, 15, 0.01), shapes = c(0.5, 0.5, 2), 
#'                        scales = c(5, 10, 10/gamma(1.5)), type = "ODE")
#' plot(IDM)
#' plot(IDM2)
#' 
#' 



probtrans_weib <- function(transMat, times, shapes, scales, type = c("prodint", "ODE")){
  type <- match.arg(type)
  #Number of transitions
  M <- sum(!is.na(transMat))
  n_times <- length(times)
  n_states <- nrow(transMat)
  if(length(shapes)!=length(scales)){
    stop("Please provide an equal amount of shape and scale parameters")
  } else if(length(shapes)!=M){
    stop("Please provide shapes and scales for each possible transition")
  }
  
  
  if(type == "prodint"){
    Haz_mat <- matrix(NA, nrow = M*length(times), ncol = 3)
    for(i in 1:M){
      Haz_mat[(1+(i-1)*n_times):(i*n_times), ] <- matrix(c(times, 
                                      -pweibull(times, shapes[i], scales[i], lower.tail = FALSE, log.p = TRUE),
                                      rep(i, length(times))), ncol = 3)
    }
    colnames(Haz_mat) <- c("time", "Haz", "trans")
    
    msfit_temp <- list(Haz = as.data.frame(Haz_mat),
                       trans = transMat)
    class(msfit_temp) <- "msfit"
    prob_out <- mstate::probtrans(msfit_temp, predt = min(times), direction = "forward", variance = FALSE)
  } else if(type == "ODE"){
    zero_adjusted <- FALSE
    if(any(shapes < 1) & (0 %in% times)){
      #Hazard is infinite for shapes < 1 at time = 0
      times <- times[times != 0]
      zero_adjusted <- TRUE
    }
    init_state <- diag(n_states)
    pars_ode <- list(shapes = shapes,
                     scales = scales, 
                     tmat = transMat,
                     M = M,
                     n_states = n_states)
    
    out_mat <- deSolve::ode(y = init_state, times = times, func = ChapKolm_fwd_mat,
                   parms = pars_ode,
                   method = "lsodes")
    
    if(zero_adjusted){
      extra_row <- rep(0, ncol(out_mat) -1)
      extra_row[seq(1, n_states^2, by = n_states+1)] <- 1
      out_mat <- rbind(c(0, extra_row), out_mat)
      times <- c(0, times)
    }
    
    #Put out_mat results into probtrans like list
    prob_out <- vector(mode = "list", length = n_states)
    for(i in 1:n_states){
      prob_out[[i]] <- data.frame(times, out_mat[, seq(i, n_states^2, by = n_states)+1])
      colnames(prob_out[[i]]) <- c("times", paste0("pstate", 1:n_states))
    }
    prob_out$trans <- transMat
    class(prob_out) <- "probtrans"
    
  }
  
  
  return(prob_out)
  
}




#' Function to use 
#' 
#' 
#' @keywords internal
#' 

ChapKolm_fwd_mat <- function(t, state, parameters, transMat) {
  #pars_ode <- 
  #list(shapes = shapes,
  #scales = scales, 
  #tmat = transMat,
  #M = M,
  #n_states = n_states)
  #Weibull hazard function
  #dweibull(t, shape, scale) / pweibull(t, shape, scale, lower.tail = FALSE)
  with(as.list(c(state, parameters)), {
    P <- matrix(state, parameters$n_states, parameters$n_states)
    A <- matrix(0, nrow = parameters$n_states, ncol = parameters$n_states)
    for(i in 1:parameters$M){
      idx <- which(parameters$tmat == i)
      haz <- dweibull(t, parameters$shapes[i], parameters$scales[i]) / pweibull(t, parameters$shapes[i], parameters$scales[i], lower.tail = FALSE)
      A[idx] <- ifelse(is.infinite(haz), .Machine$double.xmax, haz)
    }
    diag(A) <- -apply(A, 1, sum)
    # rate of change
    dP <- P %*% A
    # return the rate of change
    list(c(dP))
  }) # end with(as.list ...
}


hweibull <- function(t, shape, scale){
  return((shape/scale) * (t/scale)^(shape-1))
}


