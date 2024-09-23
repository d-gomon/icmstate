#' @title EM solver for extended illness-death model (Frydman 1995)
#' 
#' @description Solves for the cdf and transition intensities using the EM 
#' algorithm described in Frydman (1995). 
#' 
#' @param data_idx List containing data, outputted from \code{\link{msm_frydman}}
#' @param supportMSM List containing data on the support of the 1->2 transition,
#' output from supportMSM()
#' @param z Initial values for \eqn{F_{12}}{F_(12)} and \eqn{F_{13}}{F_(13)}, used to initiate the EM alg.
#' @param lambda Initial values for \eqn{\Lambda_{23}}{Lambda_(23)}, used to initiate the EM alg.
#' @param tol Tolerance of the EM algorithm. When the change in sum(abs(z)) and 
#' sum(abs(lambda)) no longer exceeds tol, the algorithm stops.
#' 
#' @return beta: Indicator whether Q subset A
#' mu: Value used in the EM algorithm, see Frydman (1995) and Notes.
#' I_Q_in_Lm_tn_star: Indicator whether Q is in [L_m, t_n^*]
#' gamma: Value used in the EM algorithm, see Frydman (1995) and Notes.
#' alpha: Indicator whether Q subset [s_j, Infinity)
#' mu_overline: Value used in the EM algorithm, see Frydman (1995) and Notes.
#' lambda: Intensity for the 2->3 transition
#' z: Mass assigned to the 1->2 and 1->3 transitions
#' iter: Number of iterations required for convergence
#'
#'
#' @references Frydman, H. (1995). Nonparametric Estimation of a Markov 
#' `Illness-Death’ Process from Interval- Censored Observations, with 
#' Application to Diabetes Survival Data. Biometrika, 82(4), 773–789. 
#' \doi{10.2307/2337344} 
#' @keywords internal
#'
#'


EM_solver <- function(data_idx, supportMSM, z, lambda, tol = 1e-8){
  
  #Initiate values for algorithm
  #Lambda should be vector of size N
  lambda_current <- rep(Inf, data_idx$N)
  #z should be vector of size I' = I + K
  z_current <- rep(Inf, supportMSM$I + data_idx$K)
  
  
  #############Pre-calculate some quantities used for EM algorithm#############
  #----------\beta_{im} = I(Q_i \subset A_m) (Equation (13))------------#
  #We obtain an M X I matrix (M = size of group 3 and 4, I = number of Q intervals in support)
  beta <- AsubsetB(supportMSM$Q_mat, data_idx$matdata_g34[, c(3,4)])
  
  #----------I(Q_i \in [L_m, t_n_star]) (Equation 13.5)------------------------#
  #We create a 3D array to store this: dim(M, N, I) representing L_m, t_n_star, Q_i respectively
  I_Q_in_Lm_tnstar <- array(0, dim = c(data_idx$M, data_idx$N, supportMSM$I))
  for(m in 1:data_idx$M){
    for(n in 1:data_idx$N){
      #Check if t_n_star < L. If so, skip to next (indicator will be 0)
      if(data_idx$t_n_star[n] < data_idx$matdata_g34[m, 3]){
        next
      }
      #for(i in 1:supportMSM$I){
      I_Q_in_Lm_tnstar[m, n, ] <- ifelse(supportMSM$Q_mat[,2] < data_idx$t_n_star[n] & supportMSM$Q_mat[, 1] >= data_idx$matdata_g34[m, 3], 1, 0)
      #}
    }
  }
  dimnames(I_Q_in_Lm_tnstar)[[1]] <- rownames(data_idx$matdata_g34)
  dimnames(I_Q_in_Lm_tnstar)[[2]] <- 1:data_idx$N
  dimnames(I_Q_in_Lm_tnstar)[[3]] <- 1:supportMSM$I
  
  #-----alpha_{ij} = 1{Q_i \in [s_j, Inf)} (Equation (14))---------------------#
  #for i = 1, ..., I, ..., I+K, j = 1, ..., J
  Q_I_plus_K <- c(supportMSM$Q_mat[,1], data_idx$e_k_star)
  alpha <- ageqb(Q_I_plus_K, data_idx$matdata_g1[,5])
  rownames(alpha) <- 1:(supportMSM$I + data_idx$K)
  
  #-------------------I(t_m \geq t_n_star) (Equation (12))---------------------#
  #Checks if t_m >= t_n_star
  t_m <- data_idx$matdata_g34[, 5]
  I_t_m_geq_t_n_star <- ageqb(t_m, data_idx$t_n_star)
  
  iter <- 0
  while(sum(abs(z_current - z)) > tol & sum(abs(lambda_current - lambda)) > tol){
    #Iterations of the EM algorithm
    z_current <- z
    lambda_current <- lambda
    
    #--------Calculate \mu_{mi} (Equation 13 in Frydman (1995))------------#
    #For each 1 \leq m \leq M we calculate:
    #--prod_{lambda}{(r_i, R_m]} = prod_N (1-lambda_n 1(lambda_n in (r_i, R_m])))--------#
    prod_lambda_G <- function(lambda, G, t_n_star, N){
      #Function to calculate Prod_\lambda {G}
      #lambda is input vector of length N (current estimate of lambda)
      #z is input vector of length I (current estimate of z)
      #G is a half-open interval (r_i, R_m] (vector of length 2)
      #t_n_star vector of length N associated with lambda
      #N = length(lambda) = length(t_n_star)
      out <- 1
      for(j in 1:N){
        if(t_n_star[j] > G[1] && t_n_star[j] <= G[2]){
          out <- out * (1-lambda[j])
        }
      }
      return(out)
    }
    mu <- matrix(NA, nrow = data_idx$M, ncol = supportMSM$I)
    rownames(mu) <- rownames(data_idx$matdata_g34)
    colnames(mu) <- 1:supportMSM$I
    for(m in 1:data_idx$M){
      for(i in 1:supportMSM$I){
        if(beta[i,m] != 0){
          #If Q_i is contained in A_m, then we calculate this product
          #Prod_lambda {(r_i, R_m]}
          mu[m, i] <- beta[i, m] * z_current[i] * prod_lambda_G(lambda = lambda_current,
                                                                G = c(supportMSM$Q_mat[i, 2], 
                                                                      data_idx$matdata_g34[m, 4]),
                                                                t_n_star = data_idx$t_n_star, N = data_idx$N)
        }else{
          mu[m, i] <- 0
        }
      }
      mu[m, ] <- mu[m, ]/sum(mu[m,])
    }
    
    #--------Calculate \gamma_{mn} (Equation 13.5 in Frydman (1995))-----------#
    gamma <- matrix(NA, nrow = data_idx$M, ncol = data_idx$N)
    for(m in 1:data_idx$M){
      for(n in 1:data_idx$N){
        gamma[m,n] <- sum(mu[m,] * I_Q_in_Lm_tnstar[m, n, ])
      }
    }
    rownames(gamma) <- rownames(data_idx$matdata_g34)
    colnames(gamma) <- 1:data_idx$N
    
    #--------Calculate \overline{\mu}_{ji} (Equation 14 in Frydman (1995))-----#
    mu_unscaled <- sweep(t(alpha), 2, z_current, '*')
    mu_overline <- mu_unscaled * 1/rowSums(mu_unscaled)
    
    #-------------------Update values of Z and Lambda--------------------------#
    #Pre-requisites
    mu_colsums <- colSums(mu)
    mu_overline_colsums <- colSums(mu_overline)
    #-------Update value of z--------------#
    #Equation (11)
    z[1:supportMSM$I] <- (mu_colsums + mu_overline_colsums[1:supportMSM$I])/data_idx$N_star
    #Equation (12)
    if(data_idx$K != 0){
      z[(supportMSM$I + data_idx$K)-((data_idx$K-1):0)] <- (data_idx$c_k + mu_overline_colsums[(supportMSM$I+1):(supportMSM$I + data_idx$K)])/data_idx$N_star
    }
    
    #-------Update value of lambda---------#
    lambda <- data_idx$d_n/(colSums(gamma * I_t_m_geq_t_n_star))
    #The line below (old) just cost me 3 months of researching nothing...
    #WTF..............!!!!!!!!!!!!!!!!!!!
    #lambda <- data_idx$d_n/(colSums(t(gamma) %*% I_t_m_geq_t_n_star))
    
    #Update iteration index
    iter <- iter+1
  }
  
  
  
  
  
  return(list(beta = beta, 
              mu = mu, 
              I_Q_in_Lm_tnstar = I_Q_in_Lm_tnstar, 
              gamma = gamma,
              alpha = alpha,
              mu_overline = mu_overline,
              lambda = lambda,
              z = z,
              iter = iter))
}


#' EM Solver for truncated data
#' 
#' @description See EM_solver
#' 
#' @keywords internal
#' @noRd
#' 

EM_solver_trunc <- function(data_idx, supportMSM, z, lambda, tol = 1e-8){
  
  #Initiate values for algorithm
  #Lambda should be vector of size N
  lambda_current <- rep(Inf, data_idx$N)
  #z should be vector of size I' = I + K
  z_current <- rep(Inf, supportMSM$I + data_idx$K)
  
  
  #############Pre-calculate some quantities used for EM algorithm#############
  #----------\beta_{im} = I(Q_i \subset A_m) (Equation (13))------------#
  #We obtain an M X I matrix (M = size of group 3 and 4, I = number of Q intervals in support)
  beta <- AsubsetB(supportMSM$Q_mat, data_idx$matdata_g34[, c(3,4)])
  
  #----------I(Q_i \in [L_m, t_n_star]) (Equation 13.5)------------------------#
  #We create a 3D array to store this: dim(M, N, I) representing L_m, t_n_star, Q_i respectively
  I_Q_in_Lm_tnstar <- array(0, dim = c(data_idx$M, data_idx$N, supportMSM$I))
  for(m in 1:data_idx$M){
    for(n in 1:data_idx$N){
      #Check if t_n_star < L. If so, skip to next (indicator will be 0)
      if(data_idx$t_n_star[n] < data_idx$matdata_g34[m, 3]){
        next
      }
      #for(i in 1:supportMSM$I){
      I_Q_in_Lm_tnstar[m, n, ] <- ifelse(supportMSM$Q_mat[,2] < data_idx$t_n_star[n] & supportMSM$Q_mat[, 1] >= data_idx$matdata_g34[m, 3], 1, 0)
      #}
    }
  }
  dimnames(I_Q_in_Lm_tnstar)[[1]] <- rownames(data_idx$matdata_g34)
  dimnames(I_Q_in_Lm_tnstar)[[2]] <- 1:data_idx$N
  dimnames(I_Q_in_Lm_tnstar)[[3]] <- 1:supportMSM$I
  
  #-----alpha_{ij} = 1{Q_i \in [s_j, Inf)} (Equation (14))---------------------#
  #for i = 1, ..., I, ..., I+K, j = 1, ..., J
  Q_I_plus_K <- c(supportMSM$Q_mat[,1], data_idx$e_k_star)
  alpha <- ageqb(Q_I_plus_K, data_idx$matdata_g1[,5])
  rownames(alpha) <- 1:(supportMSM$I + data_idx$K)
  
  #-------------------I(t_m \geq t_n_star) (Equation (12))---------------------#
  #Checks if t_m >= t_n_star
  t_m <- data_idx$matdata_g34[, 5]
  I_t_m_geq_t_n_star <- ageqb(t_m, data_idx$t_n_star)
  
  #--------delta_{ir} = I{Q_i \subset (v_r, Inf)} (Equation (15))--------------#
  #Check if support sets contained in truncation sets
  delta <- agreaterb(Q_I_plus_K, data_idx$v_r)

  
  iter <- 0
  while(sum(abs(z_current - z)) > tol & sum(abs(lambda_current - lambda)) > tol){
    #Iterations of the EM algorithm
    z_current <- z
    lambda_current <- lambda
    
    #--------Calculate \mu_{mi} (Equation 13 in Frydman (1995))------------#
    #For each 1 \leq m \leq M we calculate:
    #--prod_{lambda}{(r_i, R_m]} = prod_N (1-lambda_n 1(lambda_n in (r_i, R_m])))--------#
    prod_lambda_G <- function(lambda, G, t_n_star, N){
      #Function to calculate Prod_\lambda {G}
      #lambda is input vector of length N (current estimate of lambda)
      #z is input vector of length I (current estimate of z)
      #G is a half-open interval (r_i, R_m] (vector of length 2)
      #t_n_star vector of length N associated with lambda
      #N = length(lambda) = length(t_n_star)
      out <- 1
      for(j in 1:N){
        if(t_n_star[j] > G[1] && t_n_star[j] <= G[2]){
          out <- out * (1-lambda[j])
        }
      }
      return(out)
    }
    mu <- matrix(NA, nrow = data_idx$M, ncol = supportMSM$I)
    rownames(mu) <- rownames(data_idx$matdata_g34)
    colnames(mu) <- 1:supportMSM$I
    for(m in 1:data_idx$M){
      for(i in 1:supportMSM$I){
        if(beta[i,m] != 0){
          #If Q_i is contained in A_m, then we calculate this product
          #Prod_lambda {(r_i, R_m]}
          mu[m, i] <- beta[i, m] * z_current[i] * prod_lambda_G(lambda = lambda_current,
                                                                G = c(supportMSM$Q_mat[i, 2], 
                                                                      data_idx$matdata_g34[m, 4]),
                                                                t_n_star = data_idx$t_n_star, N = data_idx$N)
        }else{
          mu[m, i] <- 0
        }
      }
      mu[m, ] <- mu[m, ]/sum(mu[m,])
    }
    
    #--------Calculate \gamma_{mn} (Equation 13.5 in Frydman (1995))-----------#
    gamma <- matrix(NA, nrow = data_idx$M + data_idx$Z, ncol = data_idx$N)
    for(m in 1:data_idx$M){
      for(n in 1:data_idx$N){
        gamma[m,n] <- sum(mu[m,] * I_Q_in_Lm_tnstar[m, n, ])
      }
    }
    for(m in (data_idx$M + 1):(data_idx$M + data_idx$Z)){
      for(n in 1:data_idx$N){
        gamma[m,n] <- ifelse()
      }
    }
    rownames(gamma) <- rownames(data_idx$matdata_g34)
    colnames(gamma) <- 1:data_idx$N
    
    #--------Calculate \overline{\mu}_{ji} (Equation 14 in Frydman (1995))-----#
    mu_unscaled <- sweep(t(alpha), 2, z_current, '*')
    mu_overline <- mu_unscaled * 1/rowSums(mu_unscaled)
    
    #--------Calculate E[Y_{ri}] (Equation 16+ in Frydman (1995))-----#
    EY_ri <- sweep(t(1-delta), 2, z_current, '*')
    EY_ri <- EY_ri * as.numeric(1/(t(delta) %*% z_current))
    
    
    #-------------------Update values of Z and Lambda--------------------------#
    #Pre-requisites
    mu_colsums <- colSums(mu)
    mu_overline_colsums <- colSums(mu_overline)
    #-------Update value of z--------------#
    #Equation (11)
    z[1:supportMSM$I] <- (mu_colsums + mu_overline_colsums[1:supportMSM$I])/data_idx$N_star
    #Equation (12)
    if(data_idx$K != 0){
      z[(supportMSM$I + data_idx$K)-((data_idx$K-1):0)] <- (data_idx$c_k + mu_overline_colsums[(supportMSM$I+1):(supportMSM$I + data_idx$K)])/data_idx$N_star
    }
    
    #-------Update value of lambda---------#
    lambda <- data_idx$d_n/(colSums(gamma * I_t_m_geq_t_n_star))
    
    #Update iteration index
    iter <- iter+1
  }
  
  
  
  
  
  return(list(beta = beta, 
              mu = mu, 
              I_Q_in_Lm_tnstar = I_Q_in_Lm_tnstar, 
              gamma = gamma,
              alpha = alpha,
              mu_overline = mu_overline,
              lambda = lambda,
              z = z,
              iter = iter))
}


