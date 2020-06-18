context("Make sure CalibratedVA runs as expected")

test_that("M-shrinkage works with single cause labels", {
    set.seed(123)
    N_U <- 1000
    N_L <- 100
    C <- 3
    p <- c(.5, .35, .15)
    M <- matrix(NA, nrow = C, ncol = C)
    sens <- .6
    for(i in 1:C) {
        for(j in 1:C) {
            M[i,j] <- ifelse(i==j, sens, (1-sens)/(C-1))
        }
    }
    q <- as.vector(t(M) %*% p)
    A_U <- t(rmultinom(N_U, 1, q))
    A_L <- G_L <- matrix(0, nrow = N_L, ncol = C)
    for(i in 1:N_L) {
        cause <- sample(1:C, 1, prob = p)
        A_L[i,] <- t(rmultinom(1, 1, prob = M[cause,]))
        G_L[i,cause] <- 1
    }
    calibratedva_out <- calibratedva(A_U, A_L, G_L, causes = as.character(1:C),
                                    method = "mshrink", ndraws = 500, burnin = 1, thin = 1)
    samples <- calibratedva_out$samples
    p_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,1:3]
    }))
    p_rowsums <- rowSums(p_samples)
    M_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,4:12]
    }))
    M_rowsums <- do.call(c, lapply(1:nrow(M_samples), function(i) {
        M <- M_samples[i,]
        dim(M) <- c(3,3)
        return(rowSums(M))
    }))
    expect_equal(class(samples), "mcmc.list")
    expect_equal(nrow(samples[[1]]), 500)
    expect_true(max(abs(p_rowsums - 1)) < 1e-8)
    expect_true(max(abs(M_rowsums - 1)) < 1e-8)
})

test_that("M-shrinkage works with multi cause labels", {
    set.seed(123)
    N_U <- 1000
    N_L <- 100
    C <- 3
    p <- c(.5, .35, .15)
    M <- matrix(NA, nrow = C, ncol = C)
    sens <- .6
    for(i in 1:C) {
        for(j in 1:C) {
            M[i,j] <- ifelse(i==j, sens, (1-sens)/(C-1))
        }
    }
    q <- as.vector(t(M) %*% p)
    A_U <- gtools::rdirichlet(N_U, 5*q)
    G_L <- gtools::rdirichlet(N_L, p)
    ### First half of G_L will be single cause labels
    G_L[1:50,] <- t(rmultinom(50, 1, p))
    A_L <- matrix(0, nrow = N_L, ncol = C)
    for(i in 1:N_L) {
        cause <- sample(1:C, 1, prob = G_L[i,])
        A_L[i,] <- gtools::rdirichlet(1,5*M[cause,])
    }
    calibratedva_out <- calibratedva(A_U, A_L, G_L, causes = as.character(1:C),
                                     method = "mshrink", ndraws = 500, burnin = 1, thin = 1)
    samples <- calibratedva_out$samples
    p_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,1:3]
    }))
    p_rowsums <- rowSums(p_samples)
    M_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,4:12]
    }))
    M_rowsums <- do.call(c, lapply(1:nrow(M_samples), function(i) {
        M <- M_samples[i,]
        dim(M) <- c(3,3)
        return(rowSums(M))
    }))
    expect_equal(class(samples), "mcmc.list")
    expect_equal(nrow(samples[[1]]), 500)
    expect_true(max(abs(p_rowsums - 1)) < 1e-8)
    expect_true(max(abs(M_rowsums - 1)) < 1e-8)
})

test_that("M-shrinkage works with ensemble multi cause labels", {
    set.seed(123)
    N_U <- 1000
    N_L <- 100
    C <- 3
    K <- 2
    p <- c(.5, .35, .15)
    sens <- c(.6, .7)
    M_list <- lapply(1:K, function(k) {
        M <- matrix(NA, nrow = C, ncol = C)
        for(i in 1:C) {
            for(j in 1:C) {
                M[i,j] <- ifelse(i==j, sens[k], (1-sens[k])/(C-1))
            }
        } 
        return(M)
    })
    A_U <- lapply(1:K, function(k) {
        q <- t(M_list[[k]]) %*% p
        gtools::rdirichlet(N_U, 5*q)
    })
    G_L <- gtools::rdirichlet(N_L, p)
    ### First half of G_L will be single cause labels
    G_L[1:50,] <- t(rmultinom(50, 1, p))
    A_L <- lapply(1:K, function(k) {
        return(matrix(0, nrow = N_L, ncol = C))
    })
    for(i in 1:N_L) {
        cause <- sample(1:C, 1, prob = G_L[i,])
        for(k in 1:K) {
            A_L[[k]][i,]<- gtools::rdirichlet(1,5*M_list[[k]][cause,])
        }
    }
    
    calibratedva_out <- calibratedva(A_U, A_L, G_L, causes = as.character(1:C),
                                     method = "mshrink", ndraws = 500, burnin = 1, thin = 1)
    samples <- calibratedva_out$samples
    p_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,1:3]
    }))
    p_rowsums <- rowSums(p_samples)
    M_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,4:21]
    }))
    start_index <- seq(1, C^2 * K, by = C^2)
    end_index <- seq(C^2, C^2 * K, by = C^2)
    M_rowsums <- do.call(c, lapply(1:nrow(M_samples), function(i) {
        do.call(c, lapply(1:K, function(k) {
            M <- M_samples[i,start_index[k]:end_index[k]]
            dim(M) <- c(3,3)
            return(rowSums(M))
        }))
    }))
    expect_equal(class(samples), "mcmc.list")
    expect_equal(nrow(samples[[1]]), 500)
    expect_true(max(abs(p_rowsums - 1)) < 1e-8)
    expect_true(max(abs(M_rowsums - 1)) < 1e-8)
})

test_that("M-shrinkage works with different rounding factor", {
    set.seed(123)
    N_U <- 1000
    N_L <- 100
    C <- 3
    K <- 2
    p <- c(.5, .35, .15)
    sens <- c(.6, .7)
    M_list <- lapply(1:K, function(k) {
        M <- matrix(NA, nrow = C, ncol = C)
        for(i in 1:C) {
            for(j in 1:C) {
                M[i,j] <- ifelse(i==j, sens[k], (1-sens[k])/(C-1))
            }
        } 
        return(M)
    })
    A_U <- lapply(1:K, function(k) {
        q <- t(M_list[[k]]) %*% p
        gtools::rdirichlet(N_U, 5*q)
    })
    G_L <- gtools::rdirichlet(N_L, p)
    ### First half of G_L will be single cause labels
    G_L[1:50,] <- t(rmultinom(50, 1, p))
    A_L <- lapply(1:K, function(k) {
        return(matrix(0, nrow = N_L, ncol = C))
    })
    for(i in 1:N_L) {
        cause <- sample(1:C, 1, prob = G_L[i,])
        for(k in 1:K) {
            A_L[[k]][i,]<- gtools::rdirichlet(1,5*M_list[[k]][cause,])
        }
    }
    
    calibratedva_out <- calibratedva(A_U, A_L, G_L, causes = as.character(1:C),
                                     method = "mshrink", ndraws = 500, burnin = 1, thin = 1,
                                     pseudo_samplesize = 50)
    samples <- calibratedva_out$samples
    p_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,1:3]
    }))
    p_rowsums <- rowSums(p_samples)
    M_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,4:21]
    }))
    start_index <- seq(1, C^2 * K, by = C^2)
    end_index <- seq(C^2, C^2 * K, by = C^2)
    M_rowsums <- do.call(c, lapply(1:nrow(M_samples), function(i) {
        do.call(c, lapply(1:K, function(k) {
            M <- M_samples[i,start_index[k]:end_index[k]]
            dim(M) <- c(3,3)
            return(rowSums(M))
        }))
    }))
    expect_equal(class(samples), "mcmc.list")
    expect_equal(nrow(samples[[1]]), 500)
    expect_true(max(abs(p_rowsums - 1)) < 1e-8)
    expect_true(max(abs(M_rowsums - 1)) < 1e-8)
})

test_that("M-shrinkage gives good estimate of p with ensemble multi cause labels", {
    set.seed(123)
    N_U <- 1000
    N_L <- 100
    C <- 3
    K <- 2
    p <- c(.5, .35, .15)
    sens <- c(.6, .7)
    M_list <- lapply(1:K, function(k) {
        M <- matrix(NA, nrow = C, ncol = C)
        for(i in 1:C) {
            for(j in 1:C) {
                M[i,j] <- ifelse(i==j, sens[k], (1-sens[k])/(C-1))
            }
        } 
        return(M)
    })
    A_U <- lapply(1:K, function(k) {
        q <- t(M_list[[k]]) %*% p
        gtools::rdirichlet(N_U, 5*q)
    })
    G_L <- gtools::rdirichlet(N_L, p)
    ### First half of G_L will be single cause labels
    G_L[1:50,] <- t(rmultinom(50, 1, p))
    A_L <- lapply(1:K, function(k) {
        return(matrix(0, nrow = N_L, ncol = C))
    })
    for(i in 1:N_L) {
        cause <- sample(1:C, 1, prob = G_L[i,])
        for(k in 1:K) {
            A_L[[k]][i,]<- gtools::rdirichlet(1,5*M_list[[k]][cause,])
        }
    }
    library(future)
    plan(multisession)
    calibratedva_out <- calibratedva(A_U, A_L, G_L, causes = as.character(1:C),
                                     method = "mshrink", ndraws = 5000, burnin = 1000, thin = 1)
    samples <- calibratedva_out$samples
    p_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,1:3]
    }))
    p_est <- colMeans(p_samples)
    expect_true(max(abs(p_est- p)) < .05)
})

test_that("p-shrinkage works with single cause labels", {
    set.seed(123)
    N_U <- 1000
    N_L <- 100
    C <- 3
    p <- c(.5, .35, .15)
    M <- matrix(NA, nrow = C, ncol = C)
    sens <- .6
    for(i in 1:C) {
        for(j in 1:C) {
            M[i,j] <- ifelse(i==j, sens, (1-sens)/(C-1))
        }
    }
    q <- as.vector(t(M) %*% p)
    A_U <- t(rmultinom(N_U, 1, q))
    A_L <- G_L <- matrix(0, nrow = N_L, ncol = C)
    for(i in 1:N_L) {
        cause <- sample(1:C, 1, prob = p)
        A_L[i,] <- t(rmultinom(1, 1, prob = M[cause,]))
        G_L[i,cause] <- 1
    }
    calibratedva_out <- calibratedva(A_U, A_L, G_L, causes = as.character(1:C),
                                     method = "pshrink", ndraws = 500, burnin = 1, thin = 1)
    samples <- calibratedva_out$samples
    p_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,1:3]
    }))
    p_rowsums <- rowSums(p_samples)
    M_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,4:12]
    }))
    M_rowsums <- do.call(c, lapply(1:nrow(M_samples), function(i) {
        M <- M_samples[i,]
        dim(M) <- c(3,3)
        return(rowSums(M))
    }))
    expect_equal(class(samples), "mcmc.list")
    expect_equal(nrow(samples[[1]]), 500)
    expect_true(max(abs(p_rowsums - 1)) < 1e-8)
    expect_true(max(abs(M_rowsums - 1)) < 1e-8)
})

test_that("p-shrinkage works with multi cause labels", {
    set.seed(123)
    N_U <- 1000
    N_L <- 100
    C <- 3
    p <- c(.5, .35, .15)
    M <- matrix(NA, nrow = C, ncol = C)
    sens <- .6
    for(i in 1:C) {
        for(j in 1:C) {
            M[i,j] <- ifelse(i==j, sens, (1-sens)/(C-1))
        }
    }
    q <- as.vector(t(M) %*% p)
    A_U <- gtools::rdirichlet(N_U, 5*q)
    G_L <- gtools::rdirichlet(N_L, p)
    ### First half of G_L will be single cause labels
    G_L[1:50,] <- t(rmultinom(50, 1, p))
    A_L <- matrix(0, nrow = N_L, ncol = C)
    for(i in 1:N_L) {
        cause <- sample(1:C, 1, prob = G_L[i,])
        A_L[i,] <- gtools::rdirichlet(1,5*M[cause,])
    }
    calibratedva_out <- calibratedva(A_U, A_L, G_L, causes = as.character(1:C),
                                     method = "pshrink", ndraws = 500, burnin = 1, thin = 1)
    samples <- calibratedva_out$samples
    p_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,1:3]
    }))
    p_rowsums <- rowSums(p_samples)
    M_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,4:12]
    }))
    M_rowsums <- do.call(c, lapply(1:nrow(M_samples), function(i) {
        M <- M_samples[i,]
        dim(M) <- c(3,3)
        return(rowSums(M))
    }))
    expect_equal(class(samples), "mcmc.list")
    expect_equal(nrow(samples[[1]]), 500)
    expect_true(max(abs(p_rowsums - 1)) < 1e-8)
    expect_true(max(abs(M_rowsums - 1)) < 1e-8)
})

test_that("p-shrinkage works with ensemble multi cause labels", {
    set.seed(123)
    N_U <- 1000
    N_L <- 100
    C <- 3
    K <- 2
    p <- c(.5, .35, .15)
    sens <- c(.6, .7)
    M_list <- lapply(1:K, function(k) {
        M <- matrix(NA, nrow = C, ncol = C)
        for(i in 1:C) {
            for(j in 1:C) {
                M[i,j] <- ifelse(i==j, sens[k], (1-sens[k])/(C-1))
            }
        } 
        return(M)
    })
    A_U <- lapply(1:K, function(k) {
        q <- t(M_list[[k]]) %*% p
        gtools::rdirichlet(N_U, 5*q)
    })
    G_L <- gtools::rdirichlet(N_L, p)
    ### First half of G_L will be single cause labels
    G_L[1:50,] <- t(rmultinom(50, 1, p))
    A_L <- lapply(1:K, function(k) {
        return(matrix(0, nrow = N_L, ncol = C))
    })
    for(i in 1:N_L) {
        cause <- sample(1:C, 1, prob = G_L[i,])
        for(k in 1:K) {
            A_L[[k]][i,]<- gtools::rdirichlet(1,5*M_list[[k]][cause,])
        }
    }
    
    calibratedva_out <- calibratedva(A_U, A_L, G_L, causes = as.character(1:C),
                                     method = "pshrink", ndraws = 500, burnin = 1, thin = 1)
    samples <- calibratedva_out$samples
    p_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,1:3]
    }))
    p_rowsums <- rowSums(p_samples)
    M_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,4:21]
    }))
    start_index <- seq(1, C^2 * K, by = C^2)
    end_index <- seq(C^2, C^2 * K, by = C^2)
    M_rowsums <- do.call(c, lapply(1:nrow(M_samples), function(i) {
        do.call(c, lapply(1:K, function(k) {
            M <- M_samples[i,start_index[k]:end_index[k]]
            dim(M) <- c(3,3)
            return(rowSums(M))
        }))
    }))
    expect_equal(class(samples), "mcmc.list")
    expect_equal(nrow(samples[[1]]), 500)
    expect_true(max(abs(p_rowsums - 1)) < 1e-8)
    expect_true(max(abs(M_rowsums - 1)) < 1e-8)
})

test_that("p-shrinkage works with different rounding factor", {
    set.seed(123)
    N_U <- 1000
    N_L <- 100
    C <- 3
    K <- 2
    p <- c(.5, .35, .15)
    sens <- c(.6, .7)
    M_list <- lapply(1:K, function(k) {
        M <- matrix(NA, nrow = C, ncol = C)
        for(i in 1:C) {
            for(j in 1:C) {
                M[i,j] <- ifelse(i==j, sens[k], (1-sens[k])/(C-1))
            }
        } 
        return(M)
    })
    A_U <- lapply(1:K, function(k) {
        q <- t(M_list[[k]]) %*% p
        gtools::rdirichlet(N_U, 5*q)
    })
    G_L <- gtools::rdirichlet(N_L, p)
    ### First half of G_L will be single cause labels
    G_L[1:50,] <- t(rmultinom(50, 1, p))
    A_L <- lapply(1:K, function(k) {
        return(matrix(0, nrow = N_L, ncol = C))
    })
    for(i in 1:N_L) {
        cause <- sample(1:C, 1, prob = G_L[i,])
        for(k in 1:K) {
            A_L[[k]][i,]<- gtools::rdirichlet(1,5*M_list[[k]][cause,])
        }
    }
    
    calibratedva_out <- calibratedva(A_U, A_L, G_L, causes = as.character(1:C),
                                     method = "pshrink", ndraws = 500, burnin = 1, thin = 1,
                                     pseudo_samplesize = 50)
    samples <- calibratedva_out$samples
    p_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,1:3]
    }))
    p_rowsums <- rowSums(p_samples)
    M_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,4:21]
    }))
    start_index <- seq(1, C^2 * K, by = C^2)
    end_index <- seq(C^2, C^2 * K, by = C^2)
    M_rowsums <- do.call(c, lapply(1:nrow(M_samples), function(i) {
        do.call(c, lapply(1:K, function(k) {
            M <- M_samples[i,start_index[k]:end_index[k]]
            dim(M) <- c(3,3)
            return(rowSums(M))
        }))
    }))
    expect_equal(class(samples), "mcmc.list")
    expect_equal(nrow(samples[[1]]), 500)
    expect_true(max(abs(p_rowsums - 1)) < 1e-8)
    expect_true(max(abs(M_rowsums - 1)) < 1e-8)
})


test_that("p-shrinkage gives good estimate of p with ensemble multi cause labels", {
    set.seed(123)
    N_U <- 1000
    N_L <- 100
    C <- 3
    K <- 2
    p <- c(.5, .35, .15)
    sens <- c(.6, .7)
    M_list <- lapply(1:K, function(k) {
        M <- matrix(NA, nrow = C, ncol = C)
        for(i in 1:C) {
            for(j in 1:C) {
                M[i,j] <- ifelse(i==j, sens[k], (1-sens[k])/(C-1))
            }
        } 
        return(M)
    })
    A_U <- lapply(1:K, function(k) {
        q <- t(M_list[[k]]) %*% p
        gtools::rdirichlet(N_U, 5*q)
    })
    G_L <- gtools::rdirichlet(N_L, p)
    ### First half of G_L will be single cause labels
    G_L[1:50,] <- t(rmultinom(50, 1, p))
    A_L <- lapply(1:K, function(k) {
        return(matrix(0, nrow = N_L, ncol = C))
    })
    for(i in 1:N_L) {
        cause <- sample(1:C, 1, prob = G_L[i,])
        for(k in 1:K) {
            A_L[[k]][i,]<- gtools::rdirichlet(1,5*M_list[[k]][cause,])
        }
    }
    library(future)
    plan(multisession)
    calibratedva_out <- calibratedva(A_U, A_L, G_L, causes = as.character(1:C),
                                     method = "pshrink", ndraws = 5000, burnin = 1000, thin = 1)
    samples <- calibratedva_out$samples
    p_samples <- do.call(rbind, lapply(1:3, function(chain) {
        samples[[chain]][,1:3]
    }))
    p_est <- colMeans(p_samples)
    expect_true(max(abs(p_est- p)) < .05)
})
test_that("M-shrinkage cross-validation works with ensemble multi cause labels", {
    set.seed(123)
    N_U <- 1000
    N_L <- 100
    C <- 3
    K <- 2
    p <- c(.5, .35, .15)
    sens <- c(.6, .7)
    M_list <- lapply(1:K, function(k) {
        M <- matrix(NA, nrow = C, ncol = C)
        for(i in 1:C) {
            for(j in 1:C) {
                M[i,j] <- ifelse(i==j, sens[k], (1-sens[k])/(C-1))
            }
        } 
        return(M)
    })
    A_U <- lapply(1:K, function(k) {
        q <- t(M_list[[k]]) %*% p
        gtools::rdirichlet(N_U, 5*q)
    })
    G_L <- gtools::rdirichlet(N_L, p)
    ### First half of G_L will be single cause labels
    G_L[1:50,] <- t(rmultinom(50, 1, p))
    A_L <- lapply(1:K, function(k) {
        return(matrix(0, nrow = N_L, ncol = C))
    })
    for(i in 1:N_L) {
        cause <- sample(1:C, 1, prob = G_L[i,])
        for(k in 1:K) {
            A_L[[k]][i,]<- gtools::rdirichlet(1,5*M_list[[k]][cause,])
        }
    }
    log10vec <- c(-3, -2, seq(-1, 2, length = 3))
    alpha_vec <- 10^log10vec
    calibratedva_out <- tune_calibratedva(A_U, A_L, G_L, causes = as.character(1:C),
                                          method = "mshrink", alpha_vec = alpha_vec,
                                          ndraws = 50, burnin = 1, thin = 1)
    expect_equal(names(calibratedva_out),
                 c("final_model", "alpha_final", "lambda_final",
                   "waic_final", "waic_uncalib", "waic_df"))
    samples <- calibratedva_out$final_model$samples
    expect_equal(class(samples), "mcmc.list")
    waic_df <- calibratedva_out$waic_df
    expect_equal(class(waic_df), "data.frame")
    expect_equal(nrow(waic_df), length(alpha_vec))
})

test_that("p-shrinkage cross-validation works with ensemble multi cause labels", {
    set.seed(123)
    N_U <- 1000
    N_L <- 100
    C <- 3
    K <- 2
    p <- c(.5, .35, .15)
    sens <- c(.6, .7)
    M_list <- lapply(1:K, function(k) {
        M <- matrix(NA, nrow = C, ncol = C)
        for(i in 1:C) {
            for(j in 1:C) {
                M[i,j] <- ifelse(i==j, sens[k], (1-sens[k])/(C-1))
            }
        } 
        return(M)
    })
    A_U <- lapply(1:K, function(k) {
        q <- t(M_list[[k]]) %*% p
        gtools::rdirichlet(N_U, 5*q)
    })
    G_L <- gtools::rdirichlet(N_L, p)
    ### First half of G_L will be single cause labels
    G_L[1:50,] <- t(rmultinom(50, 1, p))
    A_L <- lapply(1:K, function(k) {
        return(matrix(0, nrow = N_L, ncol = C))
    })
    for(i in 1:N_L) {
        cause <- sample(1:C, 1, prob = G_L[i,])
        for(k in 1:K) {
            A_L[[k]][i,]<- gtools::rdirichlet(1,5*M_list[[k]][cause,])
        }
    }
    log10vec <- c(-3, -2, seq(-1, 2, length = 3))
    lambda_vec <- 10^log10vec
    calibratedva_out <- tune_calibratedva(A_U, A_L, G_L, causes = as.character(1:C),
                                        method = "pshrink", lambda_vec = lambda_vec,
                                        ndraws = 50, burnin = 1, thin = 1)
    expect_equal(names(calibratedva_out),
                 c("final_model", "alpha_final", "lambda_final",
                   "waic_final", "waic_uncalib", "waic_df"))
    samples <- calibratedva_out$final_model$samples
    expect_equal(class(samples), "mcmc.list")
    waic_df <- calibratedva_out$waic_df
    expect_equal(class(waic_df), "data.frame")
    expect_equal(nrow(waic_df), length(lambda_vec))
})