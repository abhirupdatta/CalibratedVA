initialize.M <- function(T.mat) {
    ### T.mat is a matrix where T_{ij} is # true COD is cause i and predicted COD is cause j
    #return(prop.table(T.mat+0.1, 1))
    return(matrix(1/nrow(T.mat),nrow(T.mat),nrow(T.mat)))
    ## adding the small number 0.1 to make sure all entries of M are +ve
    # if(sum(T==0)){
    #     M=0.94*diag(nrow(T.mat))+0.02
    #     }else{
    #         M=t(t(T.mat)/rowSums(T.mat))
    #     }
    # M
}

initialize.p <- function(v) {
    ### v is a vector with predicted counts for each COD
    #return(v / sum(v))
    rep(1/length(v),length(v))
}


sample.B <- function(M, p, v) {
    C <- length(p)
    B <- matrix(NA, ncol = C, nrow = C)
    for(i in 1:C){
        prob <- M[,i] * p
        B[,i] <- rmultinom(n = 1, size = v[i], prob = prob)
    }
    return(B)
}

sample.M <- function(B, gamma, epsilon, T.mat) {
    C <- nrow(T.mat)
    M <- matrix(NA, ncol = C, nrow = C)
    for(j in 1:nrow(M)){
        alpha <- B[j,] + gamma * epsilon + T.mat[j,]
        alpha[j] <- alpha[j] + gamma
        M[j,] <- rdirichlet(1, alpha)
    }
    return(M)
}

sample.M2 <- function(B, gamma.vec, epsilon, T.mat) {
    C <- nrow(T.mat)
    M <- matrix(NA, ncol = C, nrow = C)
    for(j in 1:nrow(M)){
        alpha <- B[j,] + gamma.vec[j] * epsilon + T.mat[j,]
        alpha[j] <- alpha[j] + gamma.vec[j]
        M[j,] <- rdirichlet(1, alpha)
    }
    return(M)
}

sample.p <- function(B, delta){
    alpha <- rowSums(B) + delta
    p.samp <- rdirichlet(1, alpha)
    return(as.vector(p.samp))
}

log.pi <- function(gamma.val, e, alpha, beta, M, tol=.00001) {
    C <- ncol(M)
    first.term <-
        C * log(gamma(C * gamma.val * e + gamma.val + tol)) -
        C * (C - 1) * log(gamma(gamma.val * e + tol)) -
        C * log(gamma(gamma.val * e + gamma.val + tol))
    second.term <- (alpha - 1) * log(gamma.val + tol) - beta * gamma.val
    third.term <-
        sum((gamma.val * e + gamma.val) * log(diag(M) + tol)) +
        sum(gamma.val * e * log(M + tol)) -
        sum(gamma.val * e * log(diag(M) + tol))
    return(first.term + second.term + third.term)
}

log.pi2 <- function(gamma.val, e, alpha, beta, M, M.row, tol=.00001) {
    C <- ncol(M)
    gammafn.q1 <- C * gamma.val * e + gamma.val + tol
    gammafn.q2 <- gamma.val * e + tol
    gammafn.q3 <- gamma.val * e + gamma.val + tol
    first.term <-
        log(gamma(gammafn.q1)) -
        (C - 1) * log(gamma(gammafn.q2)) -
        log(gamma(gammafn.q3) + tol)
    second.term <- (alpha - 1) * log(gamma.val + tol) - beta * gamma.val
    M.vec <- M[M.row,]
    M.gamma <- M.vec[M.row]
    M.notgamma <- M.vec[-M.row]
    third.term <-
        (gamma.val * e + gamma.val) * log(M.gamma + tol) +
        sum(gamma.val * e * log(M.notgamma + tol))
    return(first.term + second.term + third.term)
}


sample.gamma <- function(gamma.t, epsilon, alpha, beta, M, tau) {
    ### will use log scale for numeric stability
    # gamma.prop <- rgamma(1, shape = gamma.t * tau, rate = tau)
    # g.num <- dgamma(gamma.t, shape = gamma.prop * tau, rate = tau, log = TRUE)
    # g.denom <- dgamma(gamma.prop, shape = gamma.t * tau, rate = tau, log = TRUE)
    gamma.prop <- rlnorm(1, meanlog = log(gamma.t), sdlog = tau)
    while(gamma.prop >= 100) {
        gamma.prop <- rlnorm(1, meanlog = log(gamma.t), sdlog = tau)
    }
    g.num <- dlnorm(gamma.t, meanlog = log(gamma.prop), sdlog = tau, log = TRUE)
    g.denom <- dlnorm(gamma.prop, meanlog = log(gamma.t), sdlog = tau, log = TRUE)
    pi.num <- log.pi(gamma.prop, epsilon, alpha, beta, M)
    pi.denom <- log.pi(gamma.t, epsilon, alpha, beta, M)
    log.alpha <- g.num + pi.num - g.denom - pi.denom
    accept.prob <- min(1, exp(log.alpha))
    u <- runif(1)
    if(u <= accept.prob) {
        return(gamma.prop)
    } else {
        return(gamma.t)
    }
}

sample.gamma2 <- function(gamma.vec, epsilon, alpha, beta, M, tau.vec, max.gamma) {
    ### will use log scale for numeric stability
    # gamma.prop <- rgamma(1, shape = gamma.t * tau, rate = tau)
    # g.num <- dgamma(gamma.t, shape = gamma.prop * tau, rate = tau, log = TRUE)
    # g.denom <- dgamma(gamma.prop, shape = gamma.t * tau, rate = tau, log = TRUE)
    sampled.gammas <- sapply(seq_along(gamma.vec), function(i) {
        gamma.t <- gamma.vec[i]
        tau <- tau.vec[i]
        gamma.prop <- rlnorm(1, meanlog = log(gamma.t), sdlog = tau)
        while(gamma.prop >= max.gamma) {
            gamma.prop <- rlnorm(1, meanlog = log(gamma.t), sdlog = tau)
        }
        g.num <- dlnorm(gamma.t, meanlog = log(gamma.prop), sdlog = tau, log = TRUE)
        g.denom <- dlnorm(gamma.prop, meanlog = log(gamma.t), sdlog = tau, log = TRUE)
        pi.num <- log.pi2(gamma.prop, epsilon, alpha, beta, M, M.row = i)
        pi.denom <- log.pi2(gamma.t, epsilon, alpha, beta, M, M.row = i)
        log.alpha <- g.num + pi.num - g.denom - pi.denom
        accept.prob <- min(1, exp(log.alpha))
        u <- runif(1)
        if(u <= accept.prob) {
            return(gamma.prop)
        } else {
            return(gamma.t)
        }
    })
    return(sampled.gammas)
}

#' @title Performs Gibbs sampling for calibration
#' @description \code{revamp.sampler} takes in estimation of the underlying
#' cause of death distribution from training data, as well as a transition
#' matrix based on calibration data. Along with the prior values,
#' it will return a list of posterior samples for parameters of interest
#'
#' @param v a C x 1 matrix, where C is the number of causes. Each row
#' should contain the number of estimated causes of death from
#' the training data
#' @param T.mat a C x C matrix where the i,j entry is the number of records based
#' on the calibration set where the true COD is i and the predicted COD is j
#' @param epsilon A numeric value for the epsilon in the prior
#' @param alpha A numeric value for the alpha in the prior
#' @param beta A numeric value for the beta in the prior
#' @param tau.vec A numeric vector for the logsd for the sampling distributions
#' of the gammas
#' @param delta A numeric value for the delta in the prior
#' @param gamma.init A numeric value for the starting value of gammas
#' @param ndraws The number of posterior samples to take
#' @param max.gamma The maximum value gamma is allowed to take in
#' posterior samples. Default is 75.
#'
#' @return A list of length \code{ndraws} where each entry in the list contains
#' the posterior draw for each parameter
#'
#' @import MCMCpack
#'
#' @export
revamp.sampler <- function(v, T.mat, epsilon, alpha, beta, tau.vec, delta,
                          gamma.init, ndraws, max.gamma = 75) {
    post.samples <- vector("list", ndraws)
    post.samples[[1]]$M <- initialize.M(T.mat)
    post.samples[[1]]$p <- initialize.p(v)
    post.samples[[1]]$B <- sample.B(post.samples[[1]]$M, post.samples[[1]]$p, v)
    post.samples[[1]]$gamma <- sample.gamma2(rep(gamma.init, nrow(T.mat)), epsilon,
                                             alpha, beta, post.samples[[1]]$M, tau.vec,
                                             max.gamma)
    for(i in 2:ndraws){
        post.samples[[i]]$M <- sample.M2(post.samples[[i-1]]$B, post.samples[[i-1]]$gamma,
                                        epsilon, T.mat)
        post.samples[[i]]$p <- sample.p(post.samples[[i-1]]$B, delta)
        post.samples[[i]]$B <- sample.B(post.samples[[i]]$M, post.samples[[i]]$p, v)
        post.samples[[i]]$gamma <- sample.gamma2(post.samples[[i-1]]$gamma,
                                                 epsilon, alpha, beta,
                                                 post.samples[[i]]$M, tau.vec,
                                                 max.gamma)
        #post.samples[[i]]$gamma <- gamma.init
        #if((i%%1000)==0) print(paste("Run", i, post.samples[[i]]$gamma))
        if((i%%10000)==0) print(paste("Draw", i))
    }
    return(post.samples)
}
