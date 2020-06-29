### params says which parameters to get multi-modal information for
### either all, or just p
is_multimodal <- function(calibration, C, K = NA, cutoff = .05, params = "all") {
    if(params == "all") {
        param_vec <- 1:(C + K * C^2)
    }
    if(params == "p") {
        param_vec <- 1:C
    }
    multimodalnodes <- sapply(param_vec, function(c) {
        laplace_modes <- Modes(unlist(calibration[,c]))$modes
        if(length(laplace_modes) == 1) {
            return(FALSE)
        } else {
            min.mode <- min(laplace_modes)  
            max.mode <- max(laplace_modes)
            if(abs(max.mode - min.mode) <= cutoff) {
                return(FALSE)
            } else {
                return(TRUE)
            }
        }
        #return(is.multimodal(unlist(calibration[,c]), min.size = .4))
    })
    multimodal  <- sum(multimodalnodes) >= 1
    return(multimodal)
}

max_r_hat <- function(calibration, C, K = NA, params = "all") {
    if(params == "all") {
        param_vec <- 1:(C + K * C^2)
    }
    if(params == "p") {
        param_vec <- 1:C
    }
    nchains <- length(calibration)
    ndraws <- nrow(calibration[[1]])
    rhats <- sapply(param_vec, function(i) {
        param_mat <- matrix(NA, nrow = ndraws, ncol = nchains)
        for(chain in 1:nchains) {
            param_mat[,chain] <- calibration[[chain]][,i]
        } 
        return(rstan::Rhat(param_mat))
    })
    rhats[is.na(rhats)] <- 1
    rhat_max <- max(rhats)
    return(rhat_max)
}

pick_param <- function(waic_df, param_vec) {
    for(i in 1:length(param_vec)) {
        ### Filter df to param greater than current value
        waic_df_filtered <- filter(waic_df, param >= param_vec[i])
        ### Number of remaining param that are eligible
        neligible <- sum(waic_df_filtered$rhat_max <= 1.05 & !waic_df_filtered$multimodal)
        ### Once we get down to all eligible, get the one with the smallest WAIC
        if(neligible == 0) {
            return(param_vec[length(param_vec)])
        }
        if(neligible == nrow(waic_df_filtered)){
            best_param <- waic_df_filtered$param[which.min(waic_df_filtered$waic_calib)]
            return(best_param)
        } 
    }
}


pick_lambda_p_multimodal <- function(waic_df, lambda_vec) {
    for(i in 1:length(lambda_vec)) {
        ### Filter df to lambda greater than current value
        waic_df_filtered <- filter(waic_df, lambda >= lambda_vec[i])
        ### Number of remaining lambda that are eligible
        neligible <- sum(waic_df_filtered$rhat_max <= 1.05 & !waic_df_filtered$multimodal_p)
        ### Once we get down to all eligible, get the one with the smallest WAIC
        if(neligible == 0) {
            return(lambda_vec[length(lambda_vec)])
        }
        if(neligible == nrow(waic_df_filtered)){
            best_lambda <- waic_df_filtered$lambda[which.min(waic_df_filtered$waic_calib)]
            return(best_lambda)
        } 
    }
}


log_lik_gbql_U <- function(A_U, p, M) {
    q_mean <- t(M) %*% p
    return(as.vector(A_U %*% log(q_mean)))
}
log_lik_gbql_L <- function(A_L, G_L, M) {
    q_mean <- G_L %*% M
    return(rowSums(A_L * log(q_mean)))
}

### Individual method log likelihood
gbql_log_lik <- function(post_samples, A_U, A_L, G_L) {
    C <- ncol(A_U)
    log_lik_list <- lapply(1:length(post_samples), function(i) {
        chain_samples <- post_samples[[i]]
        param_names <- colnames(chain_samples)
        p_samples <- chain_samples[,grepl("p", param_names)]
        M_samples <- chain_samples[,grepl("M", param_names)]
        N <- nrow(A_U)
        n <- nrow(A_L)
        log_lik_chain <- matrix(NA, nrow = nrow(chain_samples), ncol = N + n)
        for(s in 1:nrow(log_lik_chain)) {
            p <- p_samples[s,]
            M_vec <- M_samples[s,]
            M_mat <- M_vec
            dim(M_mat) <- c(C, C)
            ### Adjust M_matrix incase there are 0
            M_mat[M_mat == 0] <- 1e-16
            M_mat <- M_mat / rowSums(M_mat)
            log_lik_chain[s,1:N] <- log_lik_gbql_U(A_U, p, M_mat)
            log_lik_chain[s,(N+1):(N+n)] <- log_lik_gbql_L(A_L, G_L, M_mat)
        }
        return(log_lik_chain)
    })
    log_lik_mat <- do.call(rbind, log_lik_list)
    return(log_lik_mat)
}

#### Ensemble log likelihood
gbql_ensemble_log_lik <- function(post_samples, A_U, A_L, G_L) {
    C <- dim(A_U)[2]
    N <- dim(A_U)[1]
    n <- dim(A_L)[1]
    K <- dim(A_U)[3]
    log_lik_list <- lapply(1:length(post_samples), function(i) {
        chain_samples <- post_samples[[i]]
        param_names <- colnames(chain_samples)
        p_samples <- chain_samples[,grepl("p", param_names)]
        M_samples <- chain_samples[,grepl("M", param_names)]
        #### Create M_array
        M_array <- array(NA, dim = c(nrow(M_samples),C^2, K))
        for(k in 1:K) {
            M_array[,,k] <- M_samples[,as.numeric(substr(colnames(M_samples), 7,7)) == k]
        }
        
        log_lik_chain <- matrix(0, nrow = nrow(chain_samples), ncol = N + n)
        for(s in 1:nrow(log_lik_chain)) {
            p <- p_samples[s,]
            for(k in 1:K) {
                M_vec <- M_array[s,,k]
                M_mat <- M_vec
                dim(M_mat) <- c(C, C)
                ### Adjust M_matrix incase there are 0
                M_mat[M_mat == 0] <- 1e-16
                M_mat <- M_mat / rowSums(M_mat)
                log_lik_chain[s,1:N] <- log_lik_chain[s,1:N] + log_lik_gbql_U(A_U[,,k], p, M_mat)
                log_lik_chain[s,(N+1):(N+n)] <- log_lik_chain[s,(N+1):(N+n)] + log_lik_gbql_L(A_L[,,k], G_L, M_mat)
            }
        }
        return(log_lik_chain)
    })
    log_lik_mat <- do.call(rbind, log_lik_list)
    return(log_lik_mat)
}

uncalib_log_lik <- function(post_samples, A_U, A_L, G_L, delta = 1, eps = .001) {
    C <- ncol(A_U)
    log_lik_list <- lapply(1:length(post_samples), function(i) {
        chain_samples <- post_samples[[i]]
        param_names <- colnames(chain_samples)
        S <- nrow(chain_samples)
        v <- colSums(A_U)
        #p_samples <- chain_samples[,grepl("p", param_names)]
        #M_samples <- chain_samples[,grepl("M", param_names)]
        p_samples <- rdirichlet(S, v + delta)
        C <- ncol(A_U)
        M_mat <- (1-eps) * diag(1, C) + eps / C
        N <- nrow(A_U)
        n <- nrow(A_L)
        log_lik_chain <- matrix(NA, nrow = nrow(chain_samples), ncol = N + n)
        for(s in 1:nrow(log_lik_chain)) {
            p <- p_samples[s,]
            log_lik_chain[s,1:N] <- log_lik_gbql_U(A_U, p, M_mat)
            log_lik_chain[s,(N+1):(N+n)] <- log_lik_gbql_L(A_L, G_L, M_mat)
        }
        return(log_lik_chain)
    })
    log_lik_mat <- do.call(rbind, log_lik_list)
    return(log_lik_mat)
}


uncalib_ensemble_log_lik <- function(post_samples, A_U, A_L, G_L, delta = 1, sens = .95) {
    C <- dim(A_U)[2]
    N <- dim(A_U)[1]
    n <- dim(A_L)[1]
    K <- dim(A_U)[3]
    log_lik_list <- lapply(1:length(post_samples), function(i) {
        chain_samples <- post_samples[[i]]
        S <- nrow(chain_samples)
        v <- rowSums(sapply(1:K, function(k) colSums(A_U[,,k])))
        #p_samples <- chain_samples[,grepl("p", param_names)]
        #M_samples <- chain_samples[,grepl("M", param_names)]
        p_samples <- rdirichlet(S, v + delta)
        M_mat <- matrix(NA, nrow = C, ncol = C)
        for(i in 1:C) {
            for(j in 1:C) {
                M_mat[i,j] <- ifelse(i==j, sens, (1-sens)/(C-1))
            }
        } 
        N <- nrow(A_U)
        n <- nrow(A_L)
        log_lik_chain <- matrix(0, nrow = nrow(chain_samples), ncol = N + n)
        for(s in 1:nrow(log_lik_chain)) {
            p <- p_samples[s,]
            for(k in 1:K) {
                log_lik_chain[s,1:N] <- log_lik_chain[s,1:N] + log_lik_gbql_U(A_U[,,k], p, M_mat)
                log_lik_chain[s,(N+1):(N+n)] <- log_lik_chain[s,(N+1):(N+n)] + log_lik_gbql_L(A_L[,,k], G_L, M_mat)
            }
        }
        return(log_lik_chain)
    })
    log_lik_mat <- do.call(rbind, log_lik_list)
    return(log_lik_mat)
}

quietly_get_waic <- function(...) "dummy"
quietly_get_waic_uncalib <- function(...) "dummy"

get_waic <- function(calibration, A_U, A_L, G_L, method = c("single_alg", "ensemble")[1]) {
    if(method == "single_alg") {
        log_lik_mat_calib <- gbql_log_lik(calibration,
                                          A_U = A_U,
                                          A_L = A_L,
                                          G_L = G_L)  
    }
    if(method == "ensemble"){
        log_lik_mat_calib <- gbql_ensemble_log_lik(calibration,
                                                   A_U = A_U,
                                                   A_L = A_L,
                                                   G_L = G_L) 
    }
    waic_calib <- waic(log_lik_mat_calib)$estimates[3,1] 
    return(waic_calib)
}


get_waic_uncalib <- function(calibration, A_U, A_L, G_L, method = c("single_alg", "ensemble")[1]) {
    if(method == "single_alg") {
        log_lik_mat_calib <- uncalib_log_lik(calibration,
                                             A_U = A_U,
                                             A_L = A_L,
                                             G_L = G_L)  
    }
    if(method == "ensemble"){
        log_lik_mat_calib <- uncalib_ensemble_log_lik(calibration,
                                                      A_U = A_U,
                                                      A_L = A_L,
                                                      G_L = G_L) 
    }
    waic_uncalib <- waic(log_lik_mat_calib)$estimates[3,1] 
    return(waic_uncalib)
}


.onLoad <- function(lib, pkg) {
    quietly_get_waic  <<- purrr::quietly(get_waic)
    quietly_get_waic_uncalib  <<- purrr::quietly(get_waic_uncalib)
}

