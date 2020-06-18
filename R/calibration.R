#' @title Performs Gibbs sampling for calibration
#' @description Takes in estimated causes (or cause probabilities) for both a
#' representative set of deaths without labels, and an unrepresentative set
#' of deaths with labels, and estimates the calibrated CSMF
#'
#' @param A_U When using cause of death predictions from a single algorithm, this
#' will be a matrix, where each row gives the predicted cause of death probabilities for
#' an individual death, for individuals without cause of death labels. Each column represents
#' a cause. If using the top cause, one entry in each row should be 1, while the rest should be 0.
#' When using predictions from multiple algorithms for the ensemble approach, this should
#' be a list of matrices with algorithm predictions for the same individuals, where each
#' entry in the list are predictions from a given algorithm. See examples for more information
#' @param A_L A matrix or list in the same format as A_U, but for individuals with labeled
#' causes of death. If there are no individuals with labeled causes, leave as NULL
#' @param G_L A matrix where each row represents either the true cause for an individual
#' with a labeled cause of death (i.e. if the label for individual i is cause j, then
#' G_L[i,j] will be 1, and the other entries of that row will be 0), or the probabilities
#' that each individual died of a certain cause. The rows of \code{G_L} should correspond
#' to the rows of \code{A_L} (or the rows of each element of \code{A_L} if it is a list)
#' @param causes A character vector with the names of the causes. These should correspond to
#' the columns of \code{A_U}, \code{A_L}, and \code{G_L}
#' @param method One of either "mshrink" (default) for M-shrinkage or "pshrink" for p-shrinkage
#' @param nchains The number of chains. Default is 3
#' @param ndraws Number of draws in each chain. Default is 10,000
#' @param burnin Number of burnin samples. Default is 1,000
#' @param thin Thinning parameter. Default is no thinning
#' @param pseudo_samplesize The number of pseudo samples (T) used for the 
#' Gibbs Sampler using rounding and coarsening. Default is 100.
#' @param alpha A numeric value for the alpha in the prior of gamma when using M-shrinkage.
#' Higher values (relative to beta) leads to more shrinkage. Default is 5. If using
#' the ensemble model, a vector of length K can be used (where K is the number of algorithms).
#' @param beta A numeric value for the beta in the prior of gamma when using M-shrinkage.
#' Default is .5. 
#' @param lambda A numeric value for the lambda in the prior of p for p-shrinkage.
#' Higher values leads to more shrinkage. Default is 1.
#' #' @param delta A numeric value for the delta in the prior of p. Only used for
#' M-shrinkage sampling.
#' @param epsilon A numeric value for the epsilon in the prior of M. Default is .001.
#' @param tau.vec A numeric vector for the log standard deviation for the sampling distributions
#' of the gammas. Only used for M-shrinkage sampling.
#' @param init.seed The initial seed for L’Ecuyer-CMRG RNG stream generation. Default is 123.
#'
#' @return A list with the following components.
#' \describe{
#'   \item{samples}{A \code{\link[coda]{mcmc.list}} object containing the posterior samples
#'   for p, M, and gamma (if using M-shrinkage)}
#'   \item{A_U}{The value of \code{A_U} using for the posterior samples}
#'   \item{A_L}{The value of \code{A_L} using for the posterior samples}
#'   \item{G_L}{The value of \code{G_L} using for the posterior samples}
#' }
#' @export
#' 
#' @importFrom coda mcmc mcmc.list
#' @importFrom gtools rdirichlet
#' @importFrom future.apply future_lapply
calibratedva <- function(A_U, A_L = NULL, G_L = NULL, causes,
                         method = c("mshrink", "pshrink"), nchains = 3,
                         ndraws = 10000, burnin = 1000, thin = 1,
                         pseudo_samplesize = 100,
                         alpha = 5, beta = .5,
                         lambda = 1,
                         delta = 1,
                         epsilon = .001,
                         tau = .5, 
                         init.seed = 123) {
    ### power corresponds to 1 / pseudo_samplesize
    power <- 1 / pseudo_samplesize
    if(!is.list(A_U) & !is.matrix(A_U)) {
       stop("A_U must be either a list or a matrix") 
    }
    ### See whether this is ensemble or not
    is_ensemble <- is.list(A_U) 
    ### Make sure rows of A_U sum to 1
    if(!is_ensemble) {
        A_U_row_sums <- rowSums(A_U)
    } else {
        K <- length(A_U)
        A_U_row_sums <- do.call(c, lapply(1:K, function(k) rowSums(A_U[[k]])))
    }
    if(max(abs(A_U_row_sums - 1)) > 1e-8) {
        stop("Rows of A_U must sum to 1")
    }   
    if(!is.null(A_L)) {
        ### Make sure rows of A_L sum to 1
        if(!is_ensemble) {
            A_L_row_sums <- rowSums(A_L)
        } else {
            A_L_row_sums <- do.call(c, lapply(1:K, function(k) rowSums(A_L[[k]])))
        }
        if(max(abs(A_L_row_sums - 1)) > 1e-8) {
            stop("Rows of A_L must sum to 1")
        }
        ### Make sure dimensions of A_L and G_L are the same
        if(!is_ensemble) {
            if(!identical(dim(A_L), dim(G_L))) {
                stop("G_L and A_L do not have same number of rows and columns")
            }
        } else {
            ### Make sure dimensions of A_L and G_L are the same
            for(k in 1:K) {
                if(!identical(dim(A_L[[k]]), dim(G_L))) {
                    stop(paste0("G_L and the ", k, "th algorithm for A_L do not have same number of rows and columns"))
                }
            }
        }
    }
    ### Format A_U and A_L into arrays (only if they are in list format i.e. ensemble)
    if(is_ensemble) {
       C <- length(causes)
       N_U <-  nrow(A_U[[1]])
       N_L <- nrow(A_L[[1]])
       A_U_array <- array(NA, dim = c(N_U, C, K))
       A_L_array <- array(NA, dim = c(N_L, C, K))
       for(k in 1:K) {
           A_U_array[,,k] <- A_U[[k]]
           A_L_array[,,k] <- A_L[[k]]
       }
       A_U <- A_U_array
       A_L <- A_L_array
    }
    
    ### make sure rows of G_L sum to 1
    if(!is.null(G_L)) {
        G_L_row_sums <- rowSums(G_L)
        if(max(abs(G_L_row_sums - 1)) > 1e-8) {
            stop("Rows of G_L must sum to 1")
        }
    }
    ### Set method (either mshrink or pshrink)
    method <- method[1]
    if(!(method %in% c("mshrink", "pshrink"))) {
        stop("Method must be either mshrink or pshrink")
    } else if (method == "mshrink") {
        samples <- calibratedva_mshrink(A_U = A_U, A_L = A_L, G_L = G_L,
                                        causes = causes, nchains = nchains,
                                        ndraws = ndraws, burnin = burnin, thin = thin,
                                        power = power,
                                        alpha = alpha, beta = beta,
                                        delta = delta,
                                        epsilon = epsilon,
                                        tau = tau)
    } else {
        samples <- calibratedva_pshrink(A_U = A_U, A_L = A_L, G_L = G_L,
                                        causes = causes, nchains = nchains,
                                        power = power, ndraws = ndraws, burnin = burnin, thin = thin,
                                        lambda = lambda, epsilon = epsilon)
    }
    ### output will be a list
    ### return the list of posterior samples, the inputs, and the settings used
    ### (either ensemble vs single algorithm, and single-cause vs multi-cause)
    output <- list(samples = samples, A_U = A_U, A_L = A_L, G_L = G_L)
    return(output)
}

#' @title Chooses parameter values for calibration using WAIC
#' @description Uses a grid search to choose the shrinkage parameter for calibration
#' which gives the lowest WAIC values.
#'
#' @param A_U A matrix or list of causes for individuals without labeled
#' causes of death. See \code{\link{calibratedva}}.
#' @param A_L A matrix or list of causes for individuals with labeled
#' causes of death. See \code{\link{calibratedva}}.
#' @param G_L A matrix where each row represents either the true cause for an individual
#' with a labeled cause of death. See \code{\link{calibratedva}}.
#' @param causes A character vector with the names of the causes. These should correspond to
#' the columns of \code{A_U}, \code{A_L}, and \code{G_L}
#' @param method One of either "mshrink" (default) for M-shrinkage or "pshrink" for p-shrinkage
#' @param alpha_vec If using M-shrinkage vector of values for alpha for which the function will
#' evaluate the WAIC. If not provided, the function will auto-generate a vector.
#' @param lambda_vec If using p-shrinkage vector of values for lambda for which the function will
#' evaluate the WAIC. If not provided, the function will auto-generate a vector.
#' @param which.multimodal A character specifying whether both p and M (which.multimodal = "all")
#' should be evaluated for multimodality, or just p (which.multimodal = "p")
#' @param which.rhat A character specifying whether both p and M (which.rhat = "all")
#' should be evaluated for convergence, or just p (which.rhat = "p")
#' @param ... Additional arguments passed into \code{\link{calibratedva}}.
#'
#' @return A list with the following components.
#' \describe{
#'   \item{final_model}{A list with output specified in
#'   \code{\link{calibratedva}} object containing the posterior samples
#'   for the best value of the parameter (either alpha or lambda)}
#'   \item{alpha_final}{The chosen value of alpha for the posterior samples.
#'   Null if using p-shrinkage.}
#'   \item{lambda_final}{The chosen value of lambda for the posterior samples.
#'   Null if using M-shrinkage.}
#'   \item{waic_final}{The WAIC for the posterior samples, using the chosen
#'   value of alpha or lambda}
#'   \item{waic_uncalib}{The WAIC corresponding the uncalibrated model where we
#'   assume that each COD has a sensitivity of .95}
#'   \item{waic_df}{A data frame containing a row for each parameter evaluated,
#'   which gives the WAIC, whether or not the posteriors were multimodal, and the
#'   Rhat of the posterior samples.}
#' }
#' @export
#' 
#' @import dplyr
#' @importFrom LaplacesDemon Modes
#' @importFrom ggmcmc ggs ggs_Rhat
#' @importFrom loo waic
#' @importFrom purrr quietly
tune_calibratedva <- function(A_U, A_L = NULL, G_L = NULL, causes,
                              method = c("mshrink", "pshrink"),
                              alpha_vec = NULL, lambda_vec = NULL,
                              samples_list = NULL,
                              which.multimodal = "all",
                              which.rhat = "all",
                              ...) {
    log10vec <- c(-3, -2, seq(-1, 2, length = 23))
    if(is.null(alpha_vec)) {
        alpha_vec <- 10^log10vec
    }
    if(is.null(lambda_vec)) {
        lambda_vec=10^(log10vec)
    }
    if(is.null(samples_list)) {
        if(method == "mshrink") {
            samples_list <- lapply(alpha_vec, function(alpha) {
                calibrateva_out <- calibratedva(A_U, A_L, G_L, causes, method = "mshrink",
                                                alpha = alpha, ...)
                return(calibrateva_out)
            })
        } else {
            samples_list <- lapply(lambda_vec, function(lambda) {
                calibrateva_out <- calibratedva(A_U, A_L, G_L, causes, method = "pshrink",
                                                lambda = lambda, ...)
                return(calibrateva_out)
            })
        } 
    }
    
    ### now check whether we're using ensemble or not
    is_ensemble <- is.list(A_U)
    K <- ifelse(is_ensemble, length(A_U), 1)
    ### Format A_U and A_L into arrays (only if they are in list format i.e. ensemble)
    if(is_ensemble) {
        C <- length(causes)
        N_U <-  nrow(A_U[[1]])
        N_L <-  nrow(A_L[[1]])
        A_U_array <- array(NA, dim = c(N_U, C, K))
        A_L_array <- array(NA, dim = c(N_L, C, K))
        for(k in 1:K) {
            A_U_array[,,k] <- A_U[[k]]
            A_L_array[,,k] <- A_L[[k]]
        }
        A_U <- A_U_array
        A_L <- A_L_array
    } 
    waic_df <- do.call(rbind, lapply(seq_along(samples_list), function(i) {
        calibration <- samples_list[[i]]$samples
        multimodal <- ifelse(which.multimodal == "all",
                             is_multimodal(calibration, C = C, K = K),
                             is_multimodal(calibration, C = C, params = "p"))
        if(is_ensemble) {
            waic <- quietly_get_waic(calibration,
                             A_U = A_U,
                             A_L = A_L,
                             G_L = G_L,
                             method = "ensemble")$result
        } else {
            waic <- quietly_get_waic(calibration,
                                     A_U = A_U,
                                     A_L = A_L,
                                     G_L = G_L,
                                     method = "single_alg")$result
        }
        rhat_max <- ifelse(which.rhat == "all",
                           max_r_hat(calibration, C = C, K = K),
                           max_r_hat(calibration, C = C, K = K, params = "p"))
        param <- ifelse(method == "mshrink", alpha_vec[i], lambda_vec[i])
        df <- data.frame(param = param, 
                         waic_calib = waic,
                         multimodal = multimodal, 
                         rhat_max = rhat_max)
        return(df)
    }))
    ### pick final model
    waic_df <- arrange(waic_df, param)
    if(method == "mshrink") {
        my_param <- pick_param(waic_df, param_vec = alpha_vec)
    } else {
        my_param <- pick_param(waic_df, param_vec = lambda_vec)
    }
    
    alpha_final <- switch(method == "mshrink", my_param, NULL)
    lambda_final <- switch(method == "pshrink", my_param, NULL)
    param_index <- ifelse(method == "mshrink",
                          which(alpha_vec == my_param),
                          which(lambda_vec == my_param))
    df_index <- which(waic_df$param == my_param)
    waic_final <- waic_df$waic_calib[df_index]
    final_samples <- samples_list[[param_index]]
    ### Finally, get uncalibrated WAIC
    if(is_ensemble) {
        waic_uncalib <- quietly_get_waic_uncalib(final_samples$samples,
                                                 A_U = A_U,
                                                 A_L = A_L,
                                                 G_L = G_L,
                                                 method = "ensemble")$result
    } else {
        waic_uncalib <- quietly_get_waic_uncalib(final_samples$samples,
                                                 A_U = A_U,
                                                 A_L = A_L,
                                                 G_L = G_L,
                                                 method = "single_alg")$result
    }
    return(list(final_model = final_samples,
                alpha_final = alpha_final,
                lambda_final = lambda_final,
                waic_final = waic_final,
                waic_uncalib = waic_uncalib,
                waic_df = waic_df))
}