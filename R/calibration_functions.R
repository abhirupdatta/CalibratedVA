matrix_to_array <- function(A) {
    A_array <- array(NA, dim = c(dim(A), 1 ))
    A_array[,,1] <- A
    return(A_array)
}

round_preserve_sum <- function(x, digits = 2) {
    up <- 10 ^ digits
    x <- x * up
    y <- floor(x)
    indices <- tail(order(x-y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    return(y / up)
}

create_v <- function(A_U, power) {
    ### Does rounding for pseudo data
    A_U_dim <- dim(A_U)
    C <- A_U_dim[2]
    mult_factor <- round(1 / power)
    A_U_mult <- A_U * mult_factor
    A_U_int <- t(apply(A_U_mult, 1, function(x) round_preserve_sum(x, digits = 0)))
    d_U <- matrix(NA, nrow = nrow(A_U), ncol = mult_factor)
    for(i in 1:nrow(d_U)){
        rep_factor <- as.integer(A_U_int[i,])
        rep_factor[rep_factor < 0] <- 0
        rep_factor[C] <- mult_factor - sum(rep_factor[1:(C-1)])
        d_U[i,] <- rep(1:C, rep_factor)
    }
    v <- sapply(1:C, function(c) sum(d_U==c))
    return(v)
}

create_T <- function(A_L, G_L, C, power) {
    ### Single algorithm
    ### Can be extended to ensemble
    mult_factor <- round(1 / power)
    A_L_mult <- A_L * mult_factor
    A_L_int <- t(apply(A_L_mult, 1, function(x) round_preserve_sum(x, digits = 0)))
    d_L <- matrix(NA, nrow = nrow(A_L), ncol = mult_factor)
    for(i in 1:nrow(d_L)){
        rep_factor <- as.integer(A_L_int[i,])
        rep_factor[rep_factor < 0] <- 0
        rep_factor[C] <- mult_factor - sum(rep_factor[1:(C-1)])
        d_L[i,] <- rep(1:C, rep_factor)
    }
    
    ### Make pseudo latent variables for labeled set
    G_L_mult <- G_L * mult_factor
    G_L_int <- t(apply(G_L_mult, 1, function(x) round_preserve_sum(x, digits = 0)))
    z_L <- matrix(NA, nrow = nrow(G_L), ncol = mult_factor)
    for(i in 1:nrow(z_L)){
        rep_factor <- as.integer(G_L_int[i,])
        rep_factor[rep_factor < 0] <- 0
        rep_factor[C] <- mult_factor - sum(rep_factor[1:(C-1)])
        z_L[i,] <- rep(1:C, rep_factor)
    }
    T.mat <- matrix(NA, nrow = C, ncol = C)
    for(i in 1:nrow(T.mat)) {
        for(j in 1:ncol(T.mat)) {
            T.mat[i,j] <- sum(z_L == i & d_L == j)
        }
    }
    return(T.mat)
}

create_labeled_pseudodata <- function(A_L, C, power) {
    mult_factor <- round(1 / power)
    A_L_mult <- A_L * mult_factor
    A_L_int <- t(apply(A_L_mult, 1, function(x) round_preserve_sum(x, digits = 0)))
    # d_L <- matrix(NA, nrow = nrow(A_L), ncol = mult_factor)
    # for(i in 1:nrow(d_L)){
    #     rep_factor <- as.integer(A_L_int[i,])
    #     rep_factor[rep_factor < 0] <- 0
    #     rep_factor[C] <- mult_factor - sum(rep_factor[1:(C-1)])
    #     d_L[i,] <- rep(1:C, rep_factor)
    #     A_L_int[i,] <- sapply(1:C, function(c) sum(d_L[i,]==c))
    # }
    return(A_L_int)
}

sample.B <- function(M, p, v) {
    C <- length(p)
    B <- matrix(NA, ncol = C, nrow = C)
    for(j in 1:C){
        prob <- M[,j] * p
        if(sum(prob) == 0) {
            prob <- rep(.0001, length(p))
        }
        B[,j] <- rmultinom(n = 1, size = v[j], prob = prob)
    }
    return(B)
}

sample.M <- function(B, T.mat, C, gamma.vec, epsilon = .01, power) {
    M <- matrix(NA, ncol = C, nrow = C)
    for(i in 1:nrow(M)) {
        alpha <- (B[i,]  + T.mat[i,]) * power + gamma.vec * epsilon 
        alpha[i] <- alpha[i] + gamma.vec[i]
        M[i,] <- rdirichlet(1, alpha) 
    }
    return(M)
}

sample.p <- function(B, power, delta, C){
    alpha <- rep(NA, C)
    alpha <- (rowSums(B) * power) + delta
    p.samp <- rdirichlet(1, alpha)
    return(as.vector(p.samp))
}

sample.p.ensemble.power <- function(B.array, power, delta){
    alpha <- apply(B.array,1,sum) * power + delta
    p.samp <- rdirichlet(1, alpha)
    return(as.vector(p.samp))
}

initialize.M <- function(C) {
    M <- matrix(NA, nrow = C, ncol = C)
    for(i in 1:nrow(M)){
        M[i,] <- rdirichlet(1, rep(1, C))
    }
    return(M)
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

sample.gamma <- function(gamma.vec, epsilon, alpha, beta, M, tau.vec, max.gamma) {
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

### Function to sample individual T for M for r in L
sample.T <- function(M, C, G_L_mc, A_L_int_mc, T.mat.sc) {
    if(dim(G_L_mc)[1] == 0) {
        ### If we're in single cause scenario, just return T.mat
        return(T.mat.sc)
    } else {
        B <- array(NA, dim = c(C, C, nrow(A_L_int_mc)))
        for(j in 1:C) {
            probs <- sweep(G_L_mc, MARGIN=2, M[,j], `*`)
            for(r in 1:nrow(A_L_int_mc)) {
                prob <- probs[r,]
                if(sum(prob) == 0) {
                    prob <- rep(.0001, C)
                }
                B[,j,r] <- rmultinom(n = 1, size = A_L_int_mc[r, j], prob = prob)
            }
        }
        T.mat <- apply(B, 1:2, sum) + T.mat.sc
        return(T.mat)
    }
}