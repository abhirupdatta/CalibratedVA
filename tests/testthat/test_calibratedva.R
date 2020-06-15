context("Make sure CalibratedVA runs as expected")

test_that("M-shrinkage works with single cause labels", {
    set.seed(123)
    N_U <- 1000
    N_L <- 100
    C <- 3
    p <- c(.5, .35, .15)
    A_U <- t(rmultinom(N_U, 1, p))
    M <- matrix(NA, nrow = C, ncol = C)
    sens <- .6
    for(i in 1:C) {
        for(j in 1:C) {
            M[i,j] <- ifelse(i==j, sens, (1-sens)/(C-1))
        }
    }
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
    expect_equal(class(samples), "mcmc.list")
    expect_equal(nrow(samples[[1]]), 500)
    expect_true(max(abs(p_rowsums - 1)) < 1e-8)
})

test_that("Ensemble M-shrinkage works with single cause labels", {
    set.seed(123)
    N_U <- 1000
    N_L <- 100
    C <- 3
    p <- c(.5, .35, .15)
    A_U <- t(rmultinom(N_U, 1, p))
    M <- matrix(NA, nrow = C, ncol = C)
    sens <- .6
    for(i in 1:C) {
        for(j in 1:C) {
            M[i,j] <- ifelse(i==j, sens, (1-sens)/(C-1))
        }
    }
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
    expect_equal(class(samples), "mcmc.list")
    expect_equal(nrow(samples[[1]]), 500)
    expect_true(max(abs(p_rowsums - 1)) < 1e-8)
})

test <- future_lapply(1:3, FUN = rnorm(3))
