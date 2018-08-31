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
        if(sum(prob) == 0) {
            prob <- rep(.0001, length(p))
        }
        B[,i] <- rmultinom(n = 1, size = v[i], prob = prob)
    }
    return(B)
}

create.U <- function(M.array, j.mat, causes) {
    C <- length(causes)
    U <- matrix(NA, nrow = nrow(j.mat), ncol = C)
    for(j in 1:nrow(U)) {
        for(i in 1:ncol(U)){
            ### j is indexing the vector of combinations (each row in j.mat)
            ### i is indexing the causes
            U[j,i] <- prod(sapply(1:ncol(j.mat), function(k) {
                j.k <- which(causes == j.mat[j, k])
                return(M.array[i, j.k, k])
            }))
        }
    }
    return(U)
}

sample.B.ensemble <- function(U, p, v) {
    C <- length(p)
    B <- matrix(NA, ncol = C, nrow = length(v))
    for(i in 1:length(v)){
        prob <- U[i,] * p
        if(sum(prob) == 0) {
            prob <- rep(.0001, length(prob))
        }
        B[i,] <- rmultinom(n = 1, size = v[i], prob = prob)
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

sample.M.array.ensemble <- function(B, gamma.mat, epsilon, T.array, causes, j.mat) {
    M.array <- sapply(1:ncol(j.mat), function(k) {
        T.mat <- T.array[,,k]
        gamma.vec <- gamma.mat[,k]
        C <- nrow(T.mat)
        M <- matrix(NA, ncol = C, nrow = C)
        for(i in 1:nrow(M)){
            B.cause <- B[,i]
            j.k <- j.mat[,k]
            b.vec <- sapply(causes, function(c) {
                which.j.k <- which(j.k == c)
                return(sum(B.cause[which.j.k]))
            })
            alpha <- b.vec + gamma.vec[i] * epsilon + T.mat[i,]
            alpha[i] <- alpha[i] + gamma.vec[i]
            M[i,] <- rdirichlet(1, alpha)
        }
        return(M)
    }, simplify = "array")
    
    return(M.array)
}

sample.p <- function(B, delta){
    alpha <- rowSums(B) + delta
    p.samp <- rdirichlet(1, alpha)
    return(as.vector(p.samp))
}

sample.p.ensemble <- function(B, delta){
    alpha <- colSums(B) + delta
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
revamp.sampler <- function(test.cod, calib.cod, calib.truth, causes,
                           epsilon, alpha, beta, tau.vec, delta,
                           gamma.init, ndraws, max.gamma = 75) {
    v <- sapply(causes, function(c) sum(test.cod == c))
    C <- length(causes)
    T.mat <- matrix(NA, nrow = C, ncol = C)
    for(i in 1:C){
        for(j in 1:C){
            T.mat[i,j] <- sum(calib.truth == causes[i] & calib.cod == causes[j])
        }
    }
    post.samples <- vector("list", ndraws)
    post.samples[[1]]$M <- initialize.M(T.mat)
    for(i in 1:nrow(post.samples[[1]]$M)){
        if(sum(T.mat[i,]) > 0){
            post.samples[[1]]$M[i,] <- T.mat[i,] / sum(T.mat[i,])
        }
    }
    #post.samples[[1]]$p <- initialize.p(v)
    if(sum(v) > 0){
        post.samples[[1]]$p <- v / sum(v) 
    } else {
        post.samples[[1]]$p <- initialize.p(v) 
    }
    
    names(post.samples[[1]]$p) <- causes
    post.samples[[1]]$B <- sample.B(post.samples[[1]]$M, post.samples[[1]]$p, v)
    post.samples[[1]]$gamma <- sample.gamma2(rep(gamma.init, nrow(T.mat)), epsilon,
                                             alpha, beta, post.samples[[1]]$M, tau.vec,
                                             max.gamma)
    for(i in 2:ndraws){
        post.samples[[i]]$M <- sample.M2(post.samples[[i-1]]$B, post.samples[[i-1]]$gamma,
                                        epsilon, T.mat)
        post.samples[[i]]$p <- sample.p(post.samples[[i-1]]$B, delta)
        names(post.samples[[i]]$p) <- causes
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

#' @title collect CSMF estimates from the ReVAMP sampler
#' @description \code{revampCSMF} collects draws of the CSMF estimates
#' after running the \code{revamp.sampler} function, using the
#' user supplied burn-in and thinning parameters
#' 
#' @param revamp.samples a list returned from \code{revamp.sampler}
#' @param burnin an integer specifying the number of draws you wish to discard
#' before collecting draws. Default is 1,000
#' @param thin an integer specifying the amount the draws should be thinned by
#' 
#' @return a data frame with the values of the CSMF and the names of each cause
#' for each of the draws that are obtained
#' 
#' @export
revampCSMF <- function(revamp.samples, burnin = 1E3, thin = 5) {
    causes <- names(revamp.samples[[1]]$p)
    post.p <- lapply(seq(burnin, length(revamp.samples), by = thin), function(i) {
        draw.p <- revamp.samples[[i]]$p
        return(data.frame(p = draw.p, cause = causes, draw = i))
    })
    return(do.call(rbind, post.p))
} 

##########################
### test.cod.mat will be a N x K matrix, with entry i,j denoting estimated COD for indiv. i by alg. j
### calib.cod.mat is same as above, except for calibration set
### calib.truth is a vector with true COD
### causes should be a vector of the causes one is interested in

#' @title Performs Gibbs sampling for ensemble calibration 
#' @description \code{revamp.ensemble.sampler} takes in the top estimated COD
#' for each subject in the training data, as well as the calibration data.
#' Along with the prior values, it will return a list of posterior samples
#' for parameters of interest
#'
#' @param test.cod.mat will be a N x K matrix, with entry i,j denoting estimated
#' COD (as a character)for indiv. i by alg. j
#' @param calib.cod.mat is in the same format as test.cod.mat, except for the calibration set
#' @param calib.truth is a character vector with the true COD for each subject in the
#' calibration set
#' @param causes is a character vector with the names of the causes you are interested in.
#' The order of the output vector p will correspond to this vector
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
revamp.ensemble.sampler <- function(test.cod.mat, calib.cod.mat, calib.truth, causes,
                                    epsilon, alpha, beta, tau.vec, delta,
                                    gamma.init, ndraws, max.gamma = 75) {
    ### first get matrix of all combinations of j, where the j vector forms the rows
    K <- ncol(test.cod.mat)
    C <- length(causes)
    j.mat <- expand.grid(lapply(1:K, function(k) causes), stringsAsFactors = FALSE)
    ### columns of j.mat are in order of columns of test.cod.mat
    ### Now get the V vector
    v <- sapply(1:nrow(j.mat), function(i) {
        j <- as.character(as.numeric(j.mat[i,]))
        equal.rows <- sapply(1:nrow(test.cod.mat), function(r) identical(test.cod.mat[r,], j))
        v.j <- sum(equal.rows)
        return(v.j)
    })
    ### We will have an array of T matrices now
    T.array <- array(NA, dim = c(C, C, K))
    for(k in 1:K) {
        for(i in 1:C){
            for(j in 1:C) {
                T.array[i,j,k] <- sum(calib.truth == causes[i] & calib.cod.mat[,k] == causes[j])
            }
        }
    }
    ### Initialize array of M matrices
    M.array <- sapply(1:K, function(k) {
        T.mat <- T.array[,,k]
        M.mat <- initialize.M(T.mat)
        for(i in 1:nrow(T.mat)){
            if(sum(T.mat[i,]) > 0) {
                M.mat[i,] <- T.mat[i,] / sum(T.mat[i,])
            }
        }
        return(M.mat)
    }, simplify = 'array')
    
    post.samples <- vector("list", ndraws)
    
    post.samples[[1]]$M.array <- M.array
   # post.samples[[1]]$p <- rep(1 / length(causes), length(causes))
    init.p.list <- lapply(1:ncol(test.cod.mat), function(i) {
        sapply(causes, function(c) mean(test.cod.mat[,i] == c))
    })
    init.p <- Reduce("+", init.p.list) / length(init.p.list)
    init.p <- init.p / sum(init.p)
    post.samples[[1]]$p <- init.p
   # post.samples[[1]]$p <- sapply(causes, function(c) mean(test.cod.mat[,1] == c))
    names(post.samples[[1]]$p) <- causes
    ### Make U matrix
    U <- create.U(post.samples[[1]]$M.array, j.mat, causes)
    post.samples[[1]]$U <- U
    post.samples[[1]]$B <- sample.B.ensemble(U, post.samples[[1]]$p, v)
    post.samples[[1]]$gamma.mat <- matrix(NA, nrow = C, ncol = K)
    for(k in 1:K) {
        T.mat <- T.array[,,k]
        M.mat <- post.samples[[1]]$M.array[,,k]
        post.samples[[1]]$gamma.mat[,k] <- sample.gamma2(rep(gamma.init, nrow(T.mat)), epsilon,
                                                         alpha, beta, M.mat, tau.vec,
                                                         max.gamma) 
    }
    
    for(i in 2:ndraws){
        post.samples[[i]]$M.array <- sample.M.array.ensemble(post.samples[[i-1]]$B,
                                                             post.samples[[i-1]]$gamma.mat,
                                                             epsilon, T.array, causes, j.mat)
        post.samples[[i]]$p <- sample.p.ensemble(post.samples[[i-1]]$B, delta)
        names(post.samples[[i]]$p) <- causes
        U <- create.U(post.samples[[i]]$M.array, j.mat, causes)
        post.samples[[i]]$U <- U
        post.samples[[i]]$B <- sample.B.ensemble(U, post.samples[[i]]$p, v)
        post.samples[[i]]$gamma.mat <- matrix(NA, nrow = C, ncol = K)
        for(k in 1:K) {
            M.mat <- post.samples[[i]]$M.array[,,k]
            post.samples[[i]]$gamma.mat[,k] <- sample.gamma2(post.samples[[i-1]]$gamma.mat[,k],
                                                             epsilon, alpha, beta,
                                                             M.mat, tau.vec,
                                                             max.gamma) 
        }
        #post.samples[[i]]$gamma <- gamma.init
        #if((i%%1000)==0) print(paste("Run", i, post.samples[[i]]$gamma))
        #print(paste("Draw", i))
        if((i%%10000)==0) print(paste("Draw", i))
    }
    return(post.samples)
}

##########################
#' @title Wrapper function for implementing and the ReVAMP calibration
#' @description \code{revamp_with_va} trains a VA method (currently one of
#' Tariff, InterVA, or InSilicoVA) on training data (defaults to PHMRC data),
#' and then uses the calibration and test data to implement the ReVAMP 
#' hierarchical Bayesian model. This will perform the steps shown in the 
#' package vignette in one function call
#' 
#' @param train.data A data frame/matrix which is compatible for training the given VA method. 
#' Should have COD labels and be of moderate size
#' @param calibration.data A data frame/matrix in the same format as \code{train.data}.
#' @param test.data A data frame/matrix in the same format as \code{train.data}
#' and \code{calibration.data}, except without COD labels.
#' @param train.method A character string specifying the VA training method. As of now,
#' we support "Tariff", "InSilicoVA", and "InterVA". 
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
#' @param nchains The number of chains for which the ReVAMP sampling method will be run
#' @param train.seeds An optional numeric value containing the seed which should be set
#' before implementing the given training method
#' @param sampler.seeds An optional vector of numeric values containing the seeds that should be 
#' set for each chain of the Gibbs sampler. If given, should be of same length as
#' the value of \code{nchains}
#' @param ... Additional parameters for the given training method
#'
#' @return 
#' \item{causes }{ The order of the causes in the output of the Gibbs sampler}
#' \item{v }{ The C x 1 matrix used by the ReVAMP sampler}
#' \item{T.mat }{ The C x C calibration matrix used by the ReVAMP sampler}
#' \item{posterior.results}{ A list of length \code{nchains}, which each element
#' of the list itself being a list containing the output for the ReVAMP sampler
#' for that chain} \item{train.seed }{ The seed set before implementing the given training method}
#' \item{sampler.seeds }{ The seeds set before running each chain of the ReVAMP sampler}
#' @seealso \code{\link{revamp.sampler}}
#' 
#' @import openVA
#' 
revamp_with_va <- function(train.data, calibration.data, test.data, train.method,
                           epsilon, alpha, beta, tau.vec, delta,
                           gamma.init, ndraws, max.gamma = 75, nchains = 1, train.seed = NULL,
                           sampler.seeds = NULL, ...)  {
    if(length(train.method) != 1) {
        stop("train.method should be a character vector of length 1")
    }
    if(!(train.method %in% c('Tariff', 'InSilicoVA', 'InterVA'))) {
        stop("train.method must be one of Tariff, InSilicoVA, or InterVA")
    }
    if(is.null(train.seed)){
        train.seed <- sample(1:1e6, 1, replace = F)
    }
    set.seed(train.seed)
    VA.fit <- openVA::codeVA(data = rbind(calibration.data, test.data),
                             data.type = "customize", model = train.method,
                             data.train = train.data, ...)
    if(train.method == "Tariff") {
        all.predictions <- VA.fit$causes.test[,2]
    } else if (train.method == "InSilicoVA"){
        prediction.mat <- VA.fit$indiv.prob
        all.predictions <- colnames(prediction.mat)[apply(prediction.mat, 1, which.max)]
    } else {
        all.predictions <- sapply(VA.fit$VA, function(x) {
            wholeprob <- x$wholeprob
            return(names(wholeprob)[which.max(wholeprob)])
        })
    }
    
    allcauses <- sort(unique(c(calibration.data$Cause, train.data$Cause)))
    ncauses <- length(allcauses)
    calibration.predictions <- all.predictions[1:nrow(calibration.data)]
    
    v <- sapply(allcauses, function(c) sum(all.predictions == c))
    names(v) <- allcauses
    
    T.mat <- matrix(NA, nrow = ncauses, ncol = ncauses)
    calib.truecod <- calibration.data$Cause
    for(i in 1:ncauses){
        for(j in 1:ncauses){
            T.mat[i,j] <- sum(calib.truecod == allcauses[i] & calibration.predictions == allcauses[j])
        }
    }
    colnames(T.mat) <- rownames(T.mat) <- allcauses
    if(is.null(sampler.seeds)) {
        sampler.seeds <- sample(1:1e6, nchains, replace = F)
    }
    posterior.list <- lapply(1:nchains, function(i) {
        set.seed(sampler.seeds[i])
        posterior <- revamp.sampler(v = v, T.mat = T.mat, epsilon = epsilon,
                                    alpha=alpha, beta=beta, tau.vec=tau.vec, delta=delta,
                                    gamma.init=gamma.init, ndraws = ndraws)
        return(posterior)
    })
    names(posterior.list) <- paste0("chain", 1:nchains)
    return(list(causes = allcauses, v = v, T.mat = T.mat, posterior.results = posterior.list,
                train.seed = train.seed, sampler.seeds = sampler.seeds))
}
