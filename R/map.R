gamma.h <- function(gamma.param, alpha, beta, epsilon, C) {
    #if(gamma.param > 75) {
    #    gamma.param <- 75
    #}
    if(gamma.param <= 0) {
        gamma.param <- .001
    }
    term.1 <- lgamma(C * gamma.param * epsilon + gamma.param)
    term.2 <- (C - 1) * lgamma(gamma.param * epsilon)
    term.3 <- lgamma(gamma.param * epsilon + gamma.param)
    return(term.1  - term.2 - term.3 + (alpha - 1) * log(gamma.param) - beta * gamma.param)
}

gamma.fn <- function(gamma.param, M, i, alpha, beta, epsilon) {
    C <- ncol(M)
    output <- 0
    for(j in 1:C){
        add.term <- (gamma.param * epsilon + gamma.param * (i == j) - 1) * log(M[i,j]) + gamma.h(gamma.param, alpha, beta, epsilon, C)
        output <- output + add.term
    } 
    return(-output)
}

calibva.update.em <- function(params, v, T.mat, epsilon, alpha, beta, delta){
    C <- length(v)
    p <- params[1:C]
    M <- matrix(params[(C+1):(C + C^2)], ncol = C, nrow = C, byrow = FALSE)
    gamma.vec <- params[(C + C^2 + 1):length(params)]
    
    ### Update B (E-step)
    B <- matrix(NA, ncol = C, nrow = C)
    for(i in 1:nrow(B)){
        for(j in 1:ncol(B)){
            B.ij.num <- v[j] * M[i,j] * p[i]
            B.ij.denom <- sum(M[,j] * p)
            B[i,j] <- B.ij.num /  B.ij.denom
        }
    }
    ### Update params (M-step)
    M.new <- M
    for(i in 1:C){
        for(j in 1:C){
            M.ij.num <- B[i,j] + T.mat[i,j] + gamma.vec[i] * epsilon + gamma.vec[i] * (i == j) - 1
            M.ij.denom <- sum(B[i,] + T.mat[i,]) + gamma.vec[i] * C * epsilon + gamma.vec[i] - C
            M.new[i,j] <- M.ij.num / M.ij.denom
            ### Deal with negativity
            if(M.new[i,j] < 0) {
                M.new[i,j] <- .001
            }
        }
        ### Normalize rows of M
        M.new[i,] <- M.new[i,] / sum(M.new[i,])
    }
    p.new <- p
    for(i in 1:C) {
        p.new.i.num <- sum(B[i,]) + delta - 1
        p.new.i.denom <- sum(B) + delta * C - C
        p.new[i] <- p.new.i.num / p.new.i.denom
        if(p.new[i] < 0){
            p.new[i] <- .001
        }
    }
    p.new <- p.new / sum(p.new)
    gamma.vec.new <- sapply(seq_along(gamma.vec), function(i) {
        #opt.gamma <- optim(par = gamma.vec[i], fn = gamma.fn, M = M.new, i = i)
        opt.gamma <- optim(par = gamma.vec[i], fn = gamma.fn, M = M.new, i = i,
                           alpha = alpha, beta = beta, epsilon = epsilon,
                           method = "Brent", lower = .1, upper = 50)
        #return(opt.gamma)
        return(opt.gamma$par)
    })
    new.params <- c(p.new, as.vector(M.new), gamma.vec.new)
    return(new.params)
}


#' @title Obtains the MAP estimates for parameters using SQUAREM
#' @description \code{calibva.map} takes in estimation of the underlying
#' cause of death distribution from training data, as well as a transition
#' matrix based on calibration data. It uses SQUAREM to obtain MAP estimates of the 
#' parameters p, M, and gamma 
#'
#' @param test.cod will be a vector of length N, with each entry as the estimated
#' COD (as a character)for indiv. i 
#' @param calib.cod is in the same format as test.cod, except for the calibration set
#' @param calib.truth is a character vector with the true COD for each subject in the
#' calibration set
#' @param causes is a character vector with the names of the causes you are interested in.
#' The order of the output vector p will correspond to this vector
#' @param epsilon A numeric value for the epsilon in the prior
#' @param alpha A numeric value for the alpha in the prior
#' @param beta A numeric value for the beta in the prior
#' @param delta A numeric value for the delta in the prior
#' @param init.seed An optional numeric scalar, with the initial seed 
#'
#' @return a list, with each object giving the MAP estimates for that parameter
#'
#' @import SQUAREM
#'
#' @export
calibva.map <- function(test.cod, calib.cod = NULL, calib.truth = NULL, causes,
                        epsilon = .001, alpha = 5, beta = .5, delta = 1, 
                        init.seed = NULL, tol = 1e-8, maxiter = 50E3) {
    v <- sapply(causes, function(c) sum(test.cod == c))
    C <- length(causes)
    T.mat <- matrix(NA, nrow = C, ncol = C)
    
    for(i in 1:C){
        for(j in 1:C){
            if(!is.null(calib.cod) & !is.null(calib.truth)){
                T.mat[i,j] <- sum(calib.truth == causes[i] & calib.cod == causes[j])
            } else {
                T.mat[i,j] <- 0
            }
        }
    }
    if(is.null(init.seed)){
        set.seed(123)
        init.seed <- sample(-1e6:1e6, 1, replace = F)
    } 
    set.seed(init.seed)
    ### Initialize parameters
    p0 <- CalibratedVA:::initialize.p(v)
    M0 <- CalibratedVA:::initialize.M(T.mat)
    gamma.vec0 <- rep(1, C)
    params0 <- c(p0, as.vector(M0), gamma.vec0)
    em.results <- squarem(par = params0, fixptfn = calibva.update.em, v = v, T.mat = T.mat, 
                          epsilon = epsilon, alpha = alpha, beta = beta, delta = delta,
                          control = list(tol = tol, maxiter = maxiter))
    message(paste("SQUAREM finished with", em.results$fpevals, "Evaluations"))
    p.em <- em.results$par[1:C]
    M.em <- matrix(em.results$par[(C+1):(C + C^2)], ncol = C, nrow = C, byrow = FALSE)
    gamma.em <- em.results$par[(C + C^2 + 1):length(em.results$par)]
    return(list(p = p.em, M = M.em, gamma = gamma.em))
}
