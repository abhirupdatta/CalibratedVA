initialize.M <- function(T.mat) {
  ### T.mat is a matrix where T_{ij} is # true COD is cause i and predicted COD is cause j
  #return(prop.table(T.mat+0.1, 1))
  #return(matrix(1/nrow(T.mat),nrow(T.mat),nrow(T.mat)))
  M <- matrix(NA, nrow = nrow(T.mat), ncol = ncol(T.mat))
  for(i in 1:nrow(M)){
    T.row <- T.mat[i,]
    ### Normalize row of T
    if(sum(T.row) > 0){
      alpha <- T.row / sum(T.row) + 1
    } else {
      alpha <- rep(1, length(T.row))
    }
    M[i,] <- rdirichlet(1, alpha)
  }
  return(M)
}

initialize.p <- function(v) {
  ### v is a vector with predicted counts for each COD
  #return(v / sum(v))
  #rep(1/length(v),length(v))
  ### First normalize v
  v.norm <- v / sum(v)
  ### Then add normalized version of v to 1
  alpha <- v.norm + 1
  as.vector(rdirichlet(1, alpha))
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

sample.p.ensemble.lite <- function(B.array, delta){
  alpha <- apply(B.array,1,sum) + delta
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
    lgamma(gammafn.q1) -
    (C - 1) * lgamma(gammafn.q2) -
    lgamma(gammafn.q3) 
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
#' @description \code{calibva.sampler} takes in estimation of the underlying
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
#' @param ndraws The number of posterior samples to take within each chain
#' @param nchains The number of chains
#' @param init.seeds An optional numeric vector, of length nchains, with the initial seeds for each chain
#' @param max.gamma The maximum value gamma is allowed to take in
#' posterior samples. Default is 75.
#'
#' @return a mcmc.list of length \code{nchains}, where each element is a 
#' \code{mcmc.object} containing the posterior draws for a given chain
#'
#' @import MCMCpack
#'
#' @export
calibva.sampler <- function(test.cod, calib.cod = NULL, calib.truth = NULL, causes,
                           epsilon, alpha, beta, tau.vec, delta,
                           gamma.init, ndraws, nchains = 1, init.seeds = NULL,
                           max.gamma = 75, sample.gamma = TRUE,
                           gamma.final = NULL) {
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
  if(is.null(init.seeds)){
    set.seed(123)
    init.seeds <- sample(-1e6:1e6, nchains, replace = F)
  } else {
    if(length(init.seeds) != nchains) {
      stop("The init.seed vector needs to be of length nchains")
    }
  }
  posterior.list <- lapply(seq_along(init.seeds), function(chain) {
    seed <- init.seeds[chain]
    set.seed(seed)
    post.samples <- vector("list", ndraws)
    post.samples[[1]]$M <- initialize.M(T.mat)
    post.samples[[1]]$p <- initialize.p(v)
    
    # names(post.samples[[1]]$p) <- causes
    post.samples[[1]]$B <- sample.B(post.samples[[1]]$M, post.samples[[1]]$p, v)
    if(sample.gamma == FALSE) {
      if(is.null(gamma.final)) {
        stop("Please provide a fixed value of gamma")
      } else {
        post.samples[[1]]$gamma <- gamma.final
      }
      
    } else {
      post.samples[[1]]$gamma <- sample.gamma2(rep(gamma.init, nrow(T.mat)), epsilon,
                                               alpha, beta, post.samples[[1]]$M,
                                               tau.vec, max.gamma)  
    }
    
    for(i in 2:ndraws){
      post.samples[[i]]$M <- sample.M2(post.samples[[i-1]]$B, post.samples[[i-1]]$gamma,
                                       epsilon, T.mat)
      post.samples[[i]]$p <- sample.p(post.samples[[i-1]]$B, delta)
      # names(post.samples[[i]]$p) <- causes
      post.samples[[i]]$B <- sample.B(post.samples[[i]]$M, post.samples[[i]]$p, v)
      if(sample.gamma == FALSE) {
        post.samples[[i]]$gamma <- gamma.final
      } else {
        post.samples[[i]]$gamma <- sample.gamma2(post.samples[[i-1]]$gamma,
                                                 epsilon, alpha, beta,
                                                 post.samples[[i]]$M, tau.vec,
                                                 max.gamma)
      }
      
      #post.samples[[i]]$gamma <- gamma.init
      #if((i%%1000)==0) print(paste("Run", i, post.samples[[i]]$gamma))
      if((i%%10000)==0) message(paste("Chain", chain, "Draw", i))
    }
    ### Put everything into a matrix, to be converted into an mcmc object
    ### Number of params is 2 * C ^ 2 (M matrix and B matrix) + 2 * p (gamma vector & p vector)
    n.params <- 2 * C^2 + 2 * C
    post.samples.mat <- matrix(nrow = length(post.samples), ncol = n.params)
    for(i in 1:nrow(post.samples.mat)){
      samp <- post.samples[[i]]
      p.vec <- unname(samp$p)
      M.vec <- as.vector(samp$M)
      gamma.vec <- samp$gamma
      B.vec <- as.vector(samp$B)
      post.samples.mat[i,] <- c(p.vec, M.vec, gamma.vec, B.vec)
    }
    ### Column names is first cause names (with prefix p)
    ### then M (as.vector goes by column)
    cnames <- c(paste0("p[", 1:C, "]"),
                paste0(paste0("M[", rep(1:C, C), ","), rep(1:C, each = C), "]"),
                paste0("gamma[", 1:C, "]"),
                paste0(paste0("B[", rep(1:C, C), ","), rep(1:C, each = C), "]"))
    
    colnames(post.samples.mat) <- cnames
    return(mcmc(post.samples.mat))
  })
  posterior.list <- mcmc.list(posterior.list)
  return(posterior.list)
}

##########################
### test.cod.mat will be a N x K matrix, with entry i,j denoting estimated COD for indiv. i by alg. j
### calib.cod.mat is same as above, except for calibration set
### calib.truth is a vector with true COD
### causes should be a vector of the causes one is interested in

#' @title Performs Gibbs sampling for ensemble calibration 
#' @description \code{calibva.ensemble.sampler} takes in the top estimated COD
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
calibva.ensemble.sampler <- function(test.cod.mat, calib.cod.mat, calib.truth, causes,
                                    epsilon, alpha, beta, tau.vec, delta,
                                    gamma.init, ndraws, nchains = 1, init.seeds = NULL,
                                    max.gamma = 75) {
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
        #M.mat[i,] <- T.mat[i,] / sum(T.mat[i,])
        M.mat[i,] <- rdirichlet(1, T.mat[i,] + 1)
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
    post.samples[[i]]$U <- create.U(post.samples[[i]]$M.array, j.mat, causes)
    post.samples[[i]]$B <- sample.B.ensemble(post.samples[[i]]$U, post.samples[[i]]$p, v)
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
### test.cod.mat will be a N x K matrix, with entry i,j denoting estimated COD for indiv. i by alg. j
### calib.cod.mat is same as above, except for calibration set
### calib.truth is a vector with true COD
### causes should be a vector of the causes one is interested in

#' @title Performs Gibbs sampling for ensemble.lite calibration 
#' @description \code{calibva.ensemble.lite.sampler} takes in the top estimated COD
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
#' @param nchains The number of chains. Default is 1
#' @param init.seeds An optional numeric vector, of length nchains, with the initial seeds for each chain
#' @param max.gamma The maximum value gamma is allowed to take in
#' posterior samples. Default is 75.
#'
#' @return a mcmc.list of length \code{nchains}, where each element is a 
#' \code{mcmc.object} containing the posterior draws for a given chain
#'
#' @import MCMCpack
#'
#' @export
calibva.ensemble.lite.sampler <- function(test.cod.mat, calib.cod.mat, calib.truth, causes,
                                    epsilon, alpha, beta, tau.vec, delta,
                                    gamma.init, ndraws, nchains = 1, init.seeds = NULL,
                                    max.gamma = 75) {
  ### first get matrix of all combinations of j, where the j vector forms the rows
  K <- ncol(test.cod.mat)
  C <- length(causes)
  j.mat <- expand.grid(lapply(1:K, function(k) causes), stringsAsFactors = FALSE)
  ### columns of j.mat are in order of columns of test.cod.mat
  ### Now get the V matrix
  #v <- apply(test.cod.mat,2,table)[causes,]
  v <- matrix(NA, nrow = C, ncol = K)
  for(i in 1:nrow(v)){
      for(j in 1:ncol(v)){
          v[i,j] <- sum(test.cod.mat[,j] == causes[i])
      }
  }
  row.names(v) <- NULL
  
  ### We will have an array of T matrices now
  T.array <- array(NA, dim = c(C, C, K))
  for(k in 1:K) {
    for(i in 1:C){
      for(j in 1:C) {
        T.array[i,j,k] <- sum(calib.truth == causes[i] & calib.cod.mat[,k] == causes[j])
      }
    }
  }
  if(is.null(init.seeds)){
    set.seed(123)
    init.seeds <- sample(-1e6:1e6, nchains, replace = F)
  } else {
    if(length(init.seeds) != nchains) {
      stop("The init.seed vector needs to be of length nchains")
    }
  }
  posterior.list <- lapply(seq_along(init.seeds), function(chain) {
    seed <- init.seeds[chain]
    set.seed(seed)
    ### Initialize array of M matrices
    M.array <- sapply(1:K, function(k) {
      T.mat <- T.array[,,k]
      M.mat <- initialize.M(T.mat)
    }, simplify = 'array')
    
    post.samples <- vector("list", ndraws)
    
    post.samples[[1]]$M.array <- M.array
    # post.samples[[1]]$p <- rep(1 / length(causes), length(causes))
    init.p.list <- lapply(1:ncol(v), function(i) {
      initialize.p(v[,i])
    })
    init.p <- Reduce("+", init.p.list) / length(init.p.list)
    init.p <- init.p / sum(init.p)
    post.samples[[1]]$p <- init.p
    # post.samples[[1]]$p <- sapply(causes, function(c) mean(test.cod.mat[,1] == c))
   # names(post.samples[[1]]$p) <- causes
    
    B.array=sapply(1:K, function(k) {
      sample.B(post.samples[[1]]$M.array[,,k], post.samples[[1]]$p, v[,k])    
    }, simplify="array")
    post.samples[[1]]$B.array <- B.array
    
    post.samples[[1]]$gamma.mat <- matrix(NA, nrow = C, ncol = K)
    for(k in 1:K) {
      T.mat <- T.array[,,k]
      M.mat <- post.samples[[1]]$M.array[,,k]
      post.samples[[1]]$gamma.mat[,k] <- sample.gamma2(rep(gamma.init, nrow(T.mat)), epsilon,
                                                       alpha, beta, M.mat, tau.vec,
                                                       max.gamma) 
    }
    
    for(i in 2:ndraws){
      post.samples[[i]]$M.array <- sapply(1:k, function(k)
        sample.M2(post.samples[[i-1]]$B.array[,,k], 
                  post.samples[[i-1]]$gamma.mat[,k], epsilon, T.array[,,k]), simplify="array")
      
      post.samples[[i]]$p <- sample.p.ensemble.lite(post.samples[[i-1]]$B, delta)
      
      #names(post.samples[[i]]$p) <- causes
      
      post.samples[[i]]$B.array <- sapply(1:K, function(k) {
        sample.B(post.samples[[i]]$M.array[,,k], post.samples[[i]]$p, v[,k])    
      }, simplify="array")
      
      post.samples[[i]]$gamma.mat <- matrix(NA, nrow = C, ncol = K)
      for(k in 1:K) {
        M.mat <- post.samples[[i]]$M.array[,,k]
        post.samples[[i]]$gamma.mat[,k] <- sample.gamma2(post.samples[[i-1]]$gamma.mat[,k],
                                                         epsilon, alpha, beta, M.mat,
                                                         tau.vec, max.gamma) 
      }
      #post.samples[[i]]$gamma <- gamma.init
      #if((i%%1000)==0) print(paste("Run", i, post.samples[[i]]$gamma))
      #print(paste("Draw", i))
      if((i%%10000)==0) message(paste("Chain", chain, "Draw", i))
    }
    #return(post.samples)
    ### Put everything into a matrix, to be converted into an mcmc object
    ### Number of params is 2 * k * C ^ 2 (M matrix and B matrix for each algorithm)
    ###                     + C (p vector) + K * C (gamma vector for each algorithm)
    n.params <- 2 * K * C^2 + C + K * C
    post.samples.mat <- matrix(nrow = length(post.samples), ncol = n.params)
    for(i in 1:nrow(post.samples.mat)){
      samp <- post.samples[[i]]
      p.vec <- samp$p
      M.vec <- unlist(lapply(1:K, function(k) as.vector(samp$M.array[,,k])))
      gamma.vec <- unlist(lapply(1:K, function(k) samp$gamma.mat[,k]))
      B.vec <- unlist(lapply(1:K, function(k) as.vector(samp$B.array[,,k])))
      post.samples.mat[i,] <- c(p.vec, M.vec, gamma.vec, B.vec)
    }
    ### Column names is first cause names (with prefix p)
    ### then M (as.vector goes by column)
    p.names <- paste0("p[", 1:C, "]")
    M.names <- unlist(lapply(1:K, function(k) {
      paste0(paste0("M[", rep(1:C, C), ","), rep(1:C, each = C), ",", k, "]")
    }))
    gamma.names <- unlist(lapply(1:K, function(k) {
      paste0("gamma[", 1:C, ",", k, "]")
    }))
    B.names <- unlist(lapply(1:K, function(k) {
      paste0(paste0("B[", rep(1:C, C), ","), rep(1:C, each = C), ",", k, "]")
    }))
    cnames <- c(p.names,
                M.names,
                gamma.names,
                B.names)
    
    colnames(post.samples.mat) <- cnames
    return(mcmc(post.samples.mat))
  })
  return(mcmc.list(posterior.list))
}


#' @title Obtains individual COD predictions based on posterior samples
#' @description \code{calibvaIndPredictions} uses posterior samples from the CalibVA
#' hierarchical bayesian model to obtain the posterior probability for each individual 
#' dying of given causes, conditional on predictions from a single algorithm
#' 
#' @param calibva.samples a list of posterior samples obtained from \code{calibva.sampler}
#' @param test.cod a character vector containing the predicted COD from the algorithm
#' for each subject in the test set. This should be the same as the argument supplied to \code{calibva.sampler}
#' @param causes is a character vector with the names of the causes you are interested in.
#' This should be the same as the argument supplied to \code{calibva.sampler}
#' @param burnin an integer specifying the number of draws you wish to discard
#' before collecting draws. Default is 1,000
#' @param thin an integer specifying the amount the draws should be thinned by
#' 
#' @return a list with the first element \code{topCOD} giving the most likely
#' COD for each individual, based on posterior probability, and the second element
#' \code{ind.probabilities} giving a matrix such that entry i,j is the posterior
#' probability of individual i dying from cause j
#' 
#' @export
calibvaIndPredictions <- function(calibva.samples, test.cod, causes, burnin = 1E3, thin = 5) {
  prediction.mat <- matrix(NA, nrow = length(causes), ncol = length(causes))
  ### create a prediction mapping matrix where entry i, j denotes
  ### P(truth  = i | guess = j)
  for(i in 1:nrow(prediction.mat)) {
    for(j in 1:ncol(prediction.mat)){
      post.samples <- sapply(seq(burnin, length(calibva.samples), by = thin), function(draw) {
        x <- calibva.samples[[draw]]
        M <- x$M
        p <- x$p
        num <- M[i, j] * p[i]
        denom <- sum(M[,j] * p)
        return(num / denom)
      })
      prediction.mat[i,j] <- mean(post.samples)
    }
  }
  ind.predictions <- matrix(NA, nrow = length(test.cod), ncol = length(causes))
  which.cause <- sapply(test.cod, function(c) which(causes == c))
  for(i in 1:nrow(ind.predictions)) {
    ind.predictions[i,] <- prediction.mat[,which.cause[i]]
  }
  topCOD <- sapply(1:nrow(ind.predictions), function(i) {
    preds <- ind.predictions[i,]
    return(causes[which.max(preds)])
  })
  return(list(topCOD = topCOD, ind.probabilities = ind.predictions))
}

#' @title Obtains individual COD predictions based on posterior samples from the ensemble CalibVA method
#' @description \code{calibVAIndPredictions} uses posterior samples from the CalibVA
#' hierarchical bayesian model to obtain the posterior probability for each individual 
#' dying of given causes, conditional on predictions from multiple algorithms
#' 
#' @param calibVA.samples a list of posterior samples obtained from \code{calibVA.ensemble.sampler}
#' @param test.cod.mat will be a N x K matrix, with entry i,j denoting estimated
#' COD (as a character)for indiv. i by alg. j This should be the same as the argument 
#' supplied to \code{calibVA.ensemble.sampler}
#' @param causes is a character vector with the names of the causes you are interested in.
#' This should be the same as the argument supplied to \code{calibVA.ensemble.sampler}
#' @param burnin an integer specifying the number of draws you wish to discard
#' before collecting draws. Default is 1,000
#' @param thin an integer specifying the amount the draws should be thinned by
#' 
#' @return a list with the first element \code{topCOD} giving the most likely
#' COD for each individual, based on posterior probability, and the second element
#' \code{ind.probabilities} giving a matrix such that entry i,j is the posterior
#' probability of individual i dying from cause j
#' 
#' @export
calibVAEnsembleIndPredictions <- function(calibVA.samples, test.cod.mat, causes, burnin = 1E3, thin = 5) {
  K <- ncol(test.cod.mat)
  C <- length(causes)
  j.mat <- expand.grid(lapply(1:K, function(k) causes), stringsAsFactors = FALSE)
  ### create a prediction mapping matrix where entry i, j denotes
  ### P(truth  = i | guess = j)
  prediction.mat <- matrix(NA, nrow = length(causes), ncol = nrow(j.mat))
  ### predction.mat[i,j] = pr(g_r = i | a_r^(1) == j_1,...,a_r^(K) == j_K)
  for(i in 1:nrow(prediction.mat)) {
    for(j in 1:ncol(prediction.mat)){
      post.samples <- sapply(seq(burnin, length(calibVA.samples), by = thin), function(draw) {
        x <- calibVA.samples[[draw]]
        U <- x$U
        ### Index rows of U by j, columns by cause (i)
        p <- x$p
        num <- U[j, i] * p[i]
        denom <- sum(U[j,] * p)
        return(num / denom)
      })
      prediction.mat[i,j] <- mean(post.samples)
    }
  }
  ind.predictions <- matrix(NA, nrow = nrow(test.cod.mat), ncol = length(causes))
  for(r in 1:nrow(ind.predictions)) {
    a.r <- test.cod.mat[r,]
    which.j <- which(apply(j.mat, 1, function(x) identical(as.character(x), a.r)))
    ind.predictions[r,] <- prediction.mat[,which.j]
  }
  topCOD <- sapply(1:nrow(ind.predictions), function(i) {
    preds <- ind.predictions[i,]
    return(causes[which.max(preds)])
  })
  return(list(topCOD = topCOD, ind.probabilities = ind.predictions))
}

#' @title Obtains individual COD predictions based on posterior samples from the ensemble.lite CalibVA method
#' @description \code{calibVAIndPredictions} uses posterior samples from the CalibVA
#' hierarchical bayesian model to obtain the posterior probability for each individual 
#' dying of given causes, conditional on predictions from multiple algorithms
#' 
#' @param calibVA.samples a list of posterior samples obtained from \code{calibVA.ensemble.lite.sampler}
#' @param test.cod.mat will be a N x K matrix, with entry i,j denoting estimated
#' COD (as a character)for indiv. i by alg. j This should be the same as the argument 
#' supplied to \code{calibVA.ensemble.lite.sampler}
#' @param causes is a character vector with the names of the causes you are interested in.
#' This should be the same as the argument supplied to \code{calibVA.ensemble.lite.sampler}
#' @param burnin an integer specifying the number of draws you wish to discard
#' before collecting draws. Default is 1,000
#' @param thin an integer specifying the amount the draws should be thinned by
#' 
#' @return a list with the first element \code{topCOD} giving the most likely
#' COD for each individual, based on posterior probability, and the second element
#' \code{ind.probabilities} giving a matrix such that entry i,j is the posterior
#' probability of individual i dying from cause j
#' 
#' @export
calibVAEnsembleLiteIndPredictions <- function(calibVA.samples, test.cod.mat, causes, burnin = 1E3, thin = 5) {
  K <- ncol(test.cod.mat)
  C <- length(causes)
  j.mat <- expand.grid(lapply(1:K, function(k) causes), stringsAsFactors = FALSE)
  ### create a prediction mapping matrix where entry i, j denotes
  ### P(truth  = i | guess = j)
  prediction.mat <- matrix(NA, nrow = length(causes), ncol = nrow(j.mat))
  ### predction.mat[i,j] = pr(g_r = i | a_r^(1) == j_1,...,a_r^(K) == j_K)

  ### calculating U's
  print("calculating U's")
  l=seq(burnin, length(calibVA.samples), by = thin)
  U.array=array(0,c(nrow(j.mat),C,length(l)))
  for(i in 1:length(l)) {
    draw <- l[i]
    U.array[,,i] <- create.U(calibVA.samples[[draw]]$M.array, j.mat, causes)
    }
  print("Finished calculating U's")
  
  for(i in 1:nrow(prediction.mat)) {
    for(j in 1:ncol(prediction.mat)){
      post.samples <- sapply(1:length(l), function(r) {
        draw <- l[r]
        x <- calibVA.samples[[draw]]
        U <- U.array[,,r]
        ### Index rows of U by j, columns by cause (i)
        p <- x$p
        num <- U[j, i] * p[i]
        denom <- sum(U[j,] * p)
        return(num / denom)
      })
      prediction.mat[i,j] <- mean(post.samples)
    }
  }
  ind.predictions <- matrix(NA, nrow = nrow(test.cod.mat), ncol = length(causes))
  for(r in 1:nrow(ind.predictions)) {
    a.r <- test.cod.mat[r,]
    which.j <- which(apply(j.mat, 1, function(x) identical(as.character(x), a.r)))
    ind.predictions[r,] <- prediction.mat[,which.j]
  }
  topCOD <- sapply(1:nrow(ind.predictions), function(i) {
    preds <- ind.predictions[i,]
    return(causes[which.max(preds)])
  })
  return(list(topCOD = topCOD, ind.probabilities = ind.predictions))
}

#' @title Performs a MLE estimation for the likelihood portion of the CalibVA model
#' @description Uses the \link[Rsolnp]{solnp} function to obtain MLE estimates
#' of p and M, based on the likelihood portion of the CalibVA model
#' 
#' @param test.cod will be a vector of length N, with each entry as the estimated
#' COD (as a character)for indiv. i 
#' @param calib.cod is in the same format as test.cod, except for the calibration set
#' @param calib.truth is a character vector with the true COD for each subject in the
#' calibration set
#' @param causes is a character vector with the names of the causes you are interested in.
#' The order of the output vector p will correspond to this vector
#' 
#' @return a list with the first element \code{p} being MLE estimates
#' of the CSMF vector, and the second element \code{M} being the MLE
#' estimates of the M matrix described in the method
#' 
#' @import Rsolnp
#' 
mle.calibration <- function(test.cod, calib.cod, calib.truth, causes) {
  ### all arguments should be character vectors
  ### 
  C <- length(causes)
  v <- sapply(causes, function(c) sum(test.cod == c))
  T.mat <- matrix(NA, nrow = C, ncol = C)
  for(i in 1:nrow(T.mat)) {
    for(j in 1:ncol(T.mat)){
      T.mat[i,j] <- sum(calib.truth == causes[i] & calib.cod == causes[j])
    }
  }
  calib.negloglik <- function(theta) {
    ### theta will have first 4 elements be p1,...,p4, and then last 16 elements
    ### be M11,...,M14, M21,...,M44
    p <- theta[1:C]
    M <- matrix(theta[-(1:C)], nrow = C, ncol = C, byrow = TRUE)
    q <- t(M) %*% p
    v.negloglik <- -dmultinom(v, prob = q, log = TRUE)
    T.loglik <- sapply(1:nrow(T.mat), function(i) {
      return(dmultinom(T.mat[i,], prob = M[i,], log = TRUE))
    })
    T.negloglik <- -sum(T.loglik)
    return(T.negloglik + v.negloglik)
  }
  ### Initial vector is all equal
  theta.0 <- rep(1 / length(v), C + C^2)
  ### Make sure p and rows of M sum to 1
  eq.fun <- function(theta) {
    eq.vec <- sapply(seq(1, 1 + C^2, by = C), function(i) sum(theta[i:(i + C - 1)]))
    return(eq.vec)
  }
  
  ineq.fun <- function(theta) return(theta)
  length.repeat <- length(seq(1, 1 + C^2, by = C))
  mle.calib <- solnp(theta.0, calib.negloglik, eqfun = eq.fun, ineqfun = ineq.fun,
                     eqB = rep(1, length.repeat), ineqLB = rep(0.001, C + C^2), ineqUB = rep(1, C + C^2))
  p.final <- mle.calib$pars[1:C]
  M.final <- matrix(mle.calib$pars[-(1:C)], nrow = C, ncol = C, byrow = TRUE)
  return(list(p = p.final, M = M.final))
}

#' ##########################
#' #' @title Wrapper function for implementing and the CalibVA calibration
#' #' @description \code{calibVA_with_va} trains a VA method (currently one of
#' #' Tariff, InterVA, or InSilicoVA) on training data (defaults to PHMRC data),
#' #' and then uses the calibration and test data to implement the CalibVA 
#' #' hierarchical Bayesian model. This will perform the steps shown in the 
#' #' package vignette in one function call
#' #' 
#' #' @param train.data A data frame/matrix which is compatible for training the given VA method. 
#' #' Should have COD labels and be of moderate size
#' #' @param calibration.data A data frame/matrix in the same format as \code{train.data}.
#' #' @param test.data A data frame/matrix in the same format as \code{train.data}
#' #' and \code{calibration.data}, except without COD labels.
#' #' @param train.method A character string specifying the VA training method. As of now,
#' #' we support "Tariff", "InSilicoVA", and "InterVA". 
#' #' @param epsilon A numeric value for the epsilon in the prior
#' #' @param alpha A numeric value for the alpha in the prior
#' #' @param beta A numeric value for the beta in the prior
#' #' @param tau.vec A numeric vector for the logsd for the sampling distributions
#' #' of the gammas
#' #' @param delta A numeric value for the delta in the prior
#' #' @param gamma.init A numeric value for the starting value of gammas
#' #' @param ndraws The number of posterior samples to take
#' #' @param max.gamma The maximum value gamma is allowed to take in
#' #' posterior samples. Default is 75.
#' #' @param nchains The number of chains for which the CalibVA sampling method will be run
#' #' @param train.seeds An optional numeric value containing the seed which should be set
#' #' before implementing the given training method
#' #' @param sampler.seeds An optional vector of numeric values containing the seeds that should be 
#' #' set for each chain of the Gibbs sampler. If given, should be of same length as
#' #' the value of \code{nchains}
#' #' @param ... Additional parameters for the given training method
#' #'
#' #' @return 
#' #' \item{causes }{ The order of the causes in the output of the Gibbs sampler}
#' #' \item{v }{ The C x 1 matrix used by the CalibVA sampler}
#' #' \item{T.mat }{ The C x C calibration matrix used by the CalibVA sampler}
#' #' \item{posterior.results}{ A list of length \code{nchains}, which each element
#' #' of the list itself being a list containing the output for the CalibVA sampler
#' #' for that chain} \item{train.seed }{ The seed set before implementing the given training method}
#' #' \item{sampler.seeds }{ The seeds set before running each chain of the CalibVA sampler}
#' #' @seealso \code{\link{calibva.sampler}}
#' #' 
#' #' @import openVA
#' #' 
#' revamp_with_va <- function(train.data, calibration.data, test.data, train.method,
#'                            epsilon, alpha, beta, tau.vec, delta,
#'                            gamma.init, ndraws, max.gamma = 75, nchains = 1, train.seed = NULL,
#'                            sampler.seeds = NULL, ...)  {
#'   if(length(train.method) != 1) {
#'     stop("train.method should be a character vector of length 1")
#'   }
#'   if(!(train.method %in% c('Tariff', 'InSilicoVA', 'InterVA'))) {
#'     stop("train.method must be one of Tariff, InSilicoVA, or InterVA")
#'   }
#'   if(is.null(train.seed)){
#'     train.seed <- sample(1:1e6, 1, replace = F)
#'   }
#'   set.seed(train.seed)
#'   VA.fit <- openVA::codeVA(data = rbind(calibration.data, test.data),
#'                            data.type = "customize", model = train.method,
#'                            data.train = train.data, ...)
#'   if(train.method == "Tariff") {
#'     all.predictions <- VA.fit$causes.test[,2]
#'   } else if (train.method == "InSilicoVA"){
#'     prediction.mat <- VA.fit$indiv.prob
#'     all.predictions <- colnames(prediction.mat)[apply(prediction.mat, 1, which.max)]
#'   } else {
#'     all.predictions <- sapply(VA.fit$VA, function(x) {
#'       wholeprob <- x$wholeprob
#'       return(names(wholeprob)[which.max(wholeprob)])
#'     })
#'   }
#'   
#'   allcauses <- sort(unique(c(calibration.data$Cause, train.data$Cause)))
#'   ncauses <- length(allcauses)
#'   calibration.predictions <- all.predictions[1:nrow(calibration.data)]
#'   
#'   v <- sapply(allcauses, function(c) sum(all.predictions == c))
#'   names(v) <- allcauses
#'   
#'   T.mat <- matrix(NA, nrow = ncauses, ncol = ncauses)
#'   calib.truecod <- calibration.data$Cause
#'   for(i in 1:ncauses){
#'     for(j in 1:ncauses){
#'       T.mat[i,j] <- sum(calib.truecod == allcauses[i] & calibration.predictions == allcauses[j])
#'     }
#'   }
#'   colnames(T.mat) <- rownames(T.mat) <- allcauses
#'   if(is.null(sampler.seeds)) {
#'     sampler.seeds <- sample(1:1e6, nchains, replace = F)
#'   }
#'   posterior.list <- lapply(1:nchains, function(i) {
#'     set.seed(sampler.seeds[i])
#'     posterior <- revamp.sampler(v = v, T.mat = T.mat, epsilon = epsilon,
#'                                 alpha=alpha, beta=beta, tau.vec=tau.vec, delta=delta,
#'                                 gamma.init=gamma.init, ndraws = ndraws)
#'     return(posterior)
#'   })
#'   names(posterior.list) <- paste0("chain", 1:nchains)
#'   return(list(causes = allcauses, v = v, T.mat = T.mat, posterior.results = posterior.list,
#'               train.seed = train.seed, sampler.seeds = sampler.seeds))
#' }
