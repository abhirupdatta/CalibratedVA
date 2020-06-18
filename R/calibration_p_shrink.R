calibratedva_pshrink <- function(A_U, A_L = NULL, G_L = NULL, causes,
                                 power = 1 / 100,
                                 lambda = 1, ndraws = 10000,
                                 burnin = 1000, thin = 1,
                                 epsilon = .001,
                                 print.chains = FALSE,
                                 nchains = 3,
                                 init.seed = 123) {
    ### If A_U is a matrix (dimension 2), change to array
    if(length(dim(A_U)) == 2) {
        A_U <- matrix_to_array(A_U)
        A_L <- matrix_to_array(A_L)
        K <- 1
    } else {
        K <- dim(A_U)[3]
    }
    q=apply(A_U,2,mean)                                 
    C <- length(causes)
    ### Make pseudo data
    v <- matrix(NA, nrow = C, ncol = K)
    for(k in 1:K) {
        v[,k] <- create_v(A_U[,,k], power)
    }
    
    ### Make pseudo latent variables for labeled set 
    if(is.null(A_L) | is.null(G_L)) {
        T.array <- array(0, dim = c(C, C, K))
    } else {
        ### create T matrix for those with single cause
        A_L_int <- A_L
        for(k in 1:K) {
            A_L_int[,,k] <- create_labeled_pseudodata(A_L[,,k], C, power) 
        }
        ### Make pseudo latent variables for labeled set
        is_sc <- which(rowSums(G_L == 1) == 1)
        if(length(is_sc) > 0) {
            A_L_int_mc <- A_L_int[-is_sc,,]
            if(length(dim(A_L_int_mc)) == 2) {
                A_L_int_mc <- matrix_to_array(A_L_int_mc)
            }
            G_L_mc <- G_L[-is_sc,]
            T.array.sc <- array(NA, dim = c(C, C, K))
            for(k in 1:K) {
                T.array.sc[,,k] <- create_T(A_L[is_sc,,k], G_L[is_sc,], C, power)
            }
        } else {
            A_L_int_mc <- A_L_int
            G_L_mc <- G_L
            T.array.sc <- array(0, dim=c(C,C,K))
        }
    }
    posterior.list <- future_lapply(1:nchains, function(chain){
        #seed <- init.seeds[chain]
        #set.seed(seed)
        post.samples <- vector("list", ndraws)
        ### Initialize array of M matrices
        M.array <- sapply(1:K, function(k) {
            M.mat <- initialize.M(C)
        }, simplify = 'array')
        post.samples[[1]]$M.array <- M.array
        # post.samples[[1]]$p <- rep(1 / length(causes), length(causes))
        init.p.list <- lapply(1:ncol(v), function(k) {
            as.vector(rdirichlet(1, rep(1, C)))
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
        ### Sample B for L
        if(is.null(A_L) | is.null(G_L)) {
            post.samples[[1]]$T.array <- T.array
        } else {
            post.samples[[1]]$T.array <- sapply(1:K, function(k) {
                sample.T(post.samples[[1]]$M.array[,,k], C, G_L_mc, A_L_int_mc[,,k], T.array.sc[,,k])
            }, simplify = "array")
        }
        
        for(i in 2:ndraws){
            post.samples[[i]]$M.array <- sapply(1:K, function(k) 
                sample.M(post.samples[[i-1]]$B.array[,,k], 
                         T.mat = post.samples[[i-1]]$T.array[,,k],
                         C = C,
                         power = power,
                         gamma.vec = rep(1,C),
                         epsilon = epsilon), simplify="array")
            
            post.samples[[i]]$p <- sample.p.ensemble.power(post.samples[[i-1]]$B,
                                                           power,
                                                           delta = lambda * q)
            
            post.samples[[i]]$B.array <- sapply(1:K, function(k) {
                sample.B(post.samples[[i]]$M.array[,,k], post.samples[[i]]$p, v[,k])    
            }, simplify="array")
            
            if(is.null(A_L) | is.null(G_L)) {
                post.samples[[i]]$T.array <- T.array
            } else {
                post.samples[[i]]$T.array <- sapply(1:K, function(k) {
                    sample.T(post.samples[[i]]$M.array[,,k], C, G_L_mc, A_L_int_mc[,,k], T.array.sc[,,k])
                }, simplify = "array")
            }
            if((i%%10000)==0 & print.chains) message(paste("Chain", chain, "Draw", i))
        }
        #return(post.samples)
        ### Put everything into a matrix, to be converted into an mcmc object
        ### Number of params is K * C ^ 2 (M matrix for each algorithm) + C (p vector)
        n.params <- K * C^2 + C  
        post.samples.mat <- matrix(nrow = length(post.samples), ncol = n.params)
        for(i in 1:nrow(post.samples.mat)){
            samp <- post.samples[[i]]
            p.vec <- samp$p
            M.vec <- unlist(lapply(1:K, function(k) as.vector(samp$M.array[,,k])))
            post.samples.mat[i,] <- c(p.vec, M.vec)
        }
        ### Column names is first cause names (with prefix p)
        ### then M (as.vector goes by column)
        p.names <- paste0("p[", 1:C, "]")
        if(K == 1) {
            M.names <- paste0(paste0("M[", rep(1:C, C), ","), rep(1:C, each = C), "]")
        } else {
            M.names <- unlist(lapply(1:K, function(k) {
                paste0(paste0("M[", rep(1:C, C), ","), rep(1:C, each = C), ",", k, "]")
            }))
        }
        cnames <- c(p.names, M.names)
        colnames(post.samples.mat) <- cnames
        return(mcmc(post.samples.mat))
    }, future.seed = init.seed)
    post.samples.list  <- mcmc.list(posterior.list)
    post.samples.list <- window(post.samples.list, start = burnin, thin = thin)
    return(post.samples.list)
}


# calibratedva_pshrink <- function(A_U, A_L = NULL, G_L = NULL, causes,
#                                  lambda = 1, ndraws = 10000,
#                                  burnin = 1000, thin = 1,
#                                  power = 1 / 100,
#                                  epsilon = .001,
#                                  print.chains = TRUE,
#                                  init.seeds = c(123,1234,12345)) {
#     C <- length(causes)
#     q=apply(A_U,2,mean)
#     ### Make pseudo data
#     v <- create_v(A_U, power)
#     if(is.null(A_L) | is.null(G_L)) {
#         T.mat <- matrix(0, nrow = C, ncol = C)
#     } else {
#         T.mat <- create_T(A_L, G_L, C, power)
#     }
#     posterior.list <- lapply(seq_along(init.seeds), function(chain){
#         seed <- init.seeds[chain]
#         set.seed(seed)
#         post.samples <- vector("list", ndraws)
#         post.samples[[1]]$M <- initialize.M(C)
#         post.samples[[1]]$p <- as.vector(rdirichlet(1, rep(1, C)))
#         names(post.samples[[1]]$p) <- causes
#         post.samples[[1]]$B <- sample.B(post.samples[[1]]$M, post.samples[[1]]$p, v)
#         for(i in 2:ndraws){
#             post.samples[[i]]$M <- sample.M(B = post.samples[[i-1]]$B,
#                                             T.mat = T.mat,
#                                             C = C,
#                                             power = power,
#                                             gamma.vec = rep(1,C),
#                                             epsilon = epsilon)
#             post.samples[[i]]$p <- sample.p(post.samples[[i-1]]$B, power = power,  delta = lambda*q, C = C)
#             post.samples[[i]]$B <- sample.B(post.samples[[i]]$M, post.samples[[i]]$p, v= v)
#             if((i%%10000)==0 & print.chains) message(paste("Chain", chain, "Draw", i))
#         }
#         ### Put everything into a matrix, to be converted into an mcmc object
#         ### Number of params is *C ^ 2 (M matrix) +  p (p vector)
#         n.params <- C^2 + C 
#         post.samples.mat <- matrix(nrow = length(post.samples), ncol = n.params)
#         for(i in 1:nrow(post.samples.mat)){
#             samp <- post.samples[[i]]
#             p.vec <- unname(samp$p)
#             M.vec <- as.vector(samp$M)
#             post.samples.mat[i,] <- c(p.vec, M.vec)
#         }
#         ### Column names is first cause names (with prefix p)
#         ### then M (as.vector goes by column)
#         cnames <- c(paste0("p[", 1:C, "]"),
#                     paste0(paste0("M[", rep(1:C, C), ","), rep(1:C, each = C), "]"))
#         
#         colnames(post.samples.mat) <- cnames
#         return(mcmc(post.samples.mat))
#     })
#     post.samples.list  <- mcmc.list(posterior.list)
#     post.samples.list <- window(post.samples.list, start = burnin, thin = thin)
#     return(post.samples.list)
# }
# 
# 
# calibratedva_ensemble_pshrink <- function(A_U, A_L = NULL, G_L = NULL, causes,
#                                           lambda = 1, ndraws = 10000,
#                                           burnin = 1000, thin = 1,
#                                           delta = 1,
#                                           epsilon = .001,
#                                           print.chains = TRUE,
#                                           init.seeds = c(123,1234,12345)) {
#     
#     q=apply(A_U,2,mean)
#     C <- length(causes)
#     ### Make pseudo data for unlabeled set
#     K <- dim(A_U)[3]
#     ### Make alpha vector
#     if(length(alpha) == 1) {
#         alpha <- rep(alpha, K)
#     }
#     v <- matrix(NA, nrow = C, ncol = K)
#     for(k in 1:K) {
#         v[,k] <- create_v(A_U[,,k], power)
#     }
#     
#     ### Make pseudo latent variables for labeled set
#     if(is.null(A_L) | is.null(G_L)) {
#         T.array <- array(0, dim = c(C, C, K))
#     } else {
#         T.array <- array(NA, dim = c(C, C, K))
#         for(k in 1:K) {
#             T.array[,,k] <- create_T(A_L[,,k], G_L, C, power)
#         }
#     }
#     posterior.list <- lapply(seq_along(init.seeds), function(chain) {
#         seed <- init.seeds[chain]
#         set.seed(seed)
#         ### Initialize array of M matrices
#         M.array <- sapply(1:K, function(k) {
#             T.mat <- T.array[,,k]
#             M.mat <- initialize.M(C)
#         }, simplify = 'array')
#         
#         post.samples <- vector("list", ndraws)
#         
#         post.samples[[1]]$M.array <- M.array
#         # post.samples[[1]]$p <- rep(1 / length(causes), length(causes))
#         init.p.list <- lapply(1:ncol(v), function(k) {
#             as.vector(rdirichlet(1, rep(1, C)))
#         })
#         init.p <- Reduce("+", init.p.list) / length(init.p.list)
#         init.p <- init.p / sum(init.p)
#         post.samples[[1]]$p <- init.p
#         # post.samples[[1]]$p <- sapply(causes, function(c) mean(test.cod.mat[,1] == c))
#         # names(post.samples[[1]]$p) <- causes
#         
#         B.array=sapply(1:K, function(k) {
#             sample.B(post.samples[[1]]$M.array[,,k], post.samples[[1]]$p, v[,k])    
#         }, simplify="array")
#         post.samples[[1]]$B.array <- B.array
#         for(k in 1:K) {
#             T.mat <- T.array[,,k]
#             M.mat <- post.samples[[1]]$M.array[,,k]
#         }
#         
#         for(i in 2:ndraws){
#             post.samples[[i]]$M.array <- sapply(1:K, function(k) 
#                 sample.M(post.samples[[i-1]]$B.array[,,k], 
#                          T.mat = T.array[,,k],
#                          C = C,
#                          power = power,
#                          gamma.vec = rep(1,C),
#                          epsilon = epsilon), simplify="array")
#             
#             post.samples[[i]]$p <- sample.p.ensemble.power(post.samples[[i-1]]$B,
#                                                            power,
#                                                            delta=lambda*q)
#             
#             #names(post.samples[[i]]$p) <- causes
#             
#             post.samples[[i]]$B.array <- sapply(1:K, function(k) {
#                 sample.B(post.samples[[i]]$M.array[,,k], post.samples[[i]]$p, v[,k])    
#             }, simplify="array")
#             
#             if((i%%10000)==0 & print.chains) message(paste("Chain", chain, "Draw", i))
#         }
#         #return(post.samples)
#         ### Put everything into a matrix, to be converted into an mcmc object
#         ### Number of params is 2 * k * C ^ 2 (M matrix and B matrix for each algorithm)
#         ###                     + C (p vector) + K * C (gamma vector for each algorithm)
#         n.params <- K * C^2 + C 
#         post.samples.mat <- matrix(nrow = length(post.samples), ncol = n.params)
#         for(i in 1:nrow(post.samples.mat)){
#             samp <- post.samples[[i]]
#             p.vec <- samp$p
#             M.vec <- unlist(lapply(1:K, function(k) as.vector(samp$M.array[,,k])))
#             post.samples.mat[i,] <- c(p.vec, M.vec)
#         }
#         ### Column names is first cause names (with prefix p)
#         ### then M (as.vector goes by column)
#         p.names <- paste0("p[", 1:C, "]")
#         M.names <- unlist(lapply(1:K, function(k) {
#             paste0(paste0("M[", rep(1:C, C), ","), rep(1:C, each = C), ",", k, "]")
#         }))
#         cnames <- c(p.names,
#                     M.names)
#         
#         colnames(post.samples.mat) <- cnames
#         return(mcmc(post.samples.mat))
#     })
#     post.samples.list  <- mcmc.list(posterior.list)
#     post.samples.list <- window(post.samples.list, start = burnin, thin = thin)
#     return(post.samples.list)
# }
# 
# 
# calibratedva_multi_pshrink <- function(A_U, A_L = NULL, G_L = NULL, causes,
#                                        lambda = 1, ndraws = 10000,
#                                        burnin = 1000, thin = 1,
#                                        delta = 1, gamma = 1, power = 1 / 100,
#                                        epsilon = .001,
#                                        init.seeds = c(123,1234,12345)) {
#     C <- length(causes)
#     q=colMeans(A_U)
#     ### Make pseudo data for unlabeled data
#     v <- create_V(A_U, power)
#     
#     ### psuedo data for labeled data (if we have it)
#     if(is.null(A_L) | is.null(G_L)) {
#         T.mat <- matrix(0, nrow = C, ncol = C)
#     } else {
#         ### create T matrix for those with single cause
#         A_L_int <- create_labeled_pseudodata(A_L, C, power) 
#         ### Make pseudo latent variables for labeled set
#         is_sc <- which(rowSums(G_L == 1) == 1)
#         if(length(is_sc) > 0) {
#             A_L_int_mc <- A_L_int[-is_sc,]
#             G_L_mc <- G_L[-is_sc,]
#             T.mat.sc <- create_T(A_L[is_sc,], G_L[is_sc,], C, power)
#         } else {
#             A_L_int_mc <- A_L_int
#             G_L_mc <- G_L
#             T.mat.sc <- matrix(0, nrow = C, ncol = C)
#         }
#     }
#     posterior.list <- lapply(seq_along(init.seeds), function(chain){
#         seed <- init.seeds[chain]
#         set.seed(seed)
#         post.samples <- vector("list", ndraws)
#         post.samples <- vector("list", ndraws)
#         post.samples[[1]]$M <- initialize.M(C)
#         post.samples[[1]]$p <- as.vector(rdirichlet(1, rep(1, C)))
#         #names(post.samples[[1]]$p) <- causes
#         ### sample B for U
#         post.samples[[1]]$B <- sample.B(post.samples[[1]]$M, post.samples[[1]]$p, v)
#         ### Sample B for L
#         if(is.null(A_L) | is.null(G_L)) {
#             post.samples[[1]]$T.mat <- T.mat
#         } else {
#             post.samples[[1]]$T.mat <- sample.T(post.samples[[1]]$M, C,G_L_mc, A_L_int_mc, T.mat.sc)
#         }
#         for(i in 2:ndraws){
#             post.samples[[i]]$M <- sample.M(B = post.samples[[i-1]]$B,
#                                             T.mat = post.samples[[i-1]]$T.mat,
#                                             C = C,
#                                             power = power,
#                                             gamma.vec = rep(1,C),
#                                             epsilon = epsilon)
#             post.samples[[i]]$p <- sample.p(post.samples[[i-1]]$B, power = power, delta = lambda*q, C = C)
#             #names(post.samples[[i]]$p) <- causes
#             post.samples[[i]]$B <- sample.B(post.samples[[i]]$M, post.samples[[i]]$p, v= v)
#             if(is.null(A_L) | is.null(G_L)) {
#                 post.samples[[i]]$T.mat <- T.mat
#             } else {
#                 post.samples[[i]]$T.mat <- sample.T(post.samples[[i]]$M, C, G_L_mc, A_L_int_mc, T.mat.sc)
#             }
#             if((i%%5000)==0) message(paste("Chain", chain, "Draw", i))
#         }
#         # return.index <- seq(burnin, ndraws, by = thin)
#         ### Put everything into a matrix, to be converted into an mcmc object
#         ### Number of params is 2 * C ^ 2 (M matrix and B matrix) + 2 * p (gamma vector & p vector)
#         n.params <- C^2 + C 
#         post.samples.mat <- matrix(nrow = length(post.samples), ncol = n.params)
#         for(i in 1:nrow(post.samples.mat)){
#             samp <- post.samples[[i]]
#             p.vec <- unname(samp$p)
#             M.vec <- as.vector(samp$M)
#             post.samples.mat[i,] <- c(p.vec, M.vec)
#         }
#         ### Column names is first cause names (with prefix p)
#         ### then M (as.vector goes by column)
#         cnames <- c(paste0("p[", 1:C, "]"),
#                     paste0(paste0("M[", rep(1:C, C), ","), rep(1:C, each = C), "]"))
#         
#         colnames(post.samples.mat) <- cnames
#         return(mcmc(post.samples.mat))
#     })
#     posterior.list <- mcmc.list(posterior.list)
#     posterior.list <- window(posterior.list, start = burnin, thin = thin)
#     return(posterior.list)
# }
# 
# calibratedva_multi_ensemble_pshrink <- function(A_U, A_L = NULL, G_L = NULL,
#                                                 lambda = 1, causes, ndraws = 10000,
#                                                 burnin = 1000, thin = 1,
#                                                 delta = 1, gamma = 1, power = 1 / 100,
#                                                 epsilon = .001,
#                                                 print.chains = TRUE,
#                                                 init.seeds = c(123,1234,12345)) {
#     C <- length(causes)
#     q=apply(A_U,2,mean)
#     ### Make pseudo data for unlabeled data
#     K <- dim(A_U)[3]
#     ### Alpha vector
#     if(length(alpha) == 1) {
#         alpha <- rep(alpha, K)
#     }
#     v <- matrix(NA, nrow = C, ncol = K)
#     for(k in 1:K) {
#         v[,k] <- create_v(A_U[,,k], power)
#     }
#     ### psuedo data for labeled data (if we have it)
#     ### psuedo data for labeled data (if we have it)
#     if(is.null(A_L) | is.null(G_L)) {
#         T.array <- array(0, dim = c(C, C, K))
#     } else {
#         ### create T matrix for those with single cause
#         A_L_int <- A_L
#         for(k in 1:K) {
#             A_L_int[,,k] <- create_labeled_pseudodata(A_L[,,k], C, power) 
#         }
#         ### Make pseudo latent variables for labeled set
#         is_sc <- which(rowSums(G_L == 1) == 1)
#         if(length(is_sc) > 0) {
#             A_L_int_mc <- A_L_int[-is_sc,]
#             G_L_mc <- G_L[-is_sc,]
#             T.array.sc <- array(NA, dim = c(C, C, K))
#             for(k in 1:K) {
#                 T.array.sc[,,k] <-create_T(A_L[is_sc,,k], G_L[is_sc,,k], C, power)
#             }
#         } else {
#             A_L_int_mc <- A_L_int
#             G_L_mc <- G_L
#             T.array.sc <- array(0, dim=c(C,C,K))
#         }
#     }
#     posterior.list <- lapply(seq_along(init.seeds), function(chain){
#         seed <- init.seeds[chain]
#         set.seed(seed)
#         post.samples <- vector("list", ndraws)
#         M.array <- sapply(1:K, function(k) {
#             M.mat <- initialize.M(C)
#         }, simplify = 'array')
#         
#         post.samples[[1]]$M.array <- M.array
#         init.p.list <- lapply(1:ncol(v), function(k) {
#             as.vector(rdirichlet(1, rep(1, C)))
#         })
#         init.p <- Reduce("+", init.p.list) / length(init.p.list)
#         init.p <- init.p / sum(init.p)
#         post.samples[[1]]$p <- init.p
#         #names(post.samples[[1]]$p) <- causes
#         ### sample B for U
#         B.array <- sapply(1:K, function(k) {
#             sample.B(post.samples[[1]]$M.array[,,k], post.samples[[1]]$p, v[,k])
#         }, simplify="array")
#         post.samples[[1]]$B.array <- B.array
#         ### Sample B for L
#         if(is.null(A_L) | is.null(G_L)) {
#             post.samples[[1]]$T.array <- T.array
#         } else {
#             post.samples[[1]]$T.array <- sapply(1:K, function(k) {
#                 sample.T(post.samples[[1]]$M.array[,,k], C, G_L_mc, A_L_int_mc[,,k], T.array.sc[,,k])
#             }, simplify = "array")
#         }
#         for(i in 2:ndraws){
#             post.samples[[i]]$M.array <- sapply(1:K, function(k)
#                 sample.M(post.samples[[i-1]]$B.array[,,k],
#                          T.mat = post.samples[[i-1]]$T.array[,,k],
#                          C = C,
#                          power = power,
#                          gamma.vec = rep(1,C),
#                          epsilon = epsilon), simplify="array")
#             post.samples[[i]]$p <- sample.p.ensemble.power(post.samples[[i-1]]$B,
#                                                            power,
#                                                            delta=lambda*q)
#             #names(post.samples[[i]]$p) <- causes
#             post.samples[[i]]$B.array <- sapply(1:K, function(k) {
#                 sample.B(post.samples[[i]]$M.array[,,k], post.samples[[i]]$p, v[,k])
#             }, simplify="array")
#             if(is.null(A_L) | is.null(G_L)) {
#                 post.samples[[i]]$T.array <- T.array
#             } else {
#                 post.samples[[i]]$T.array <- sapply(1:K, function(k) {
#                     sample.T(post.samples[[i]]$M.array[,,k], C, G_L_mc, A_L_int_mc[,,k], T.array.sc[,,k])
#                 }, simplify = "array")
#             }
#             if((i%%5000)==0 & print.chains) message(paste("Chain", chain, "Draw", i))
#         }
#         ### Put everything into a matrix, to be converted into an mcmc object
#         ### Number of params is 2 * C ^ 2 (M matrix and B matrix) + 2 * p (gamma vector & p vector)
#         n.params <-  K * C^2 + C 
#         post.samples.mat <- matrix(nrow = length(post.samples), ncol = n.params)
#         for(i in 1:nrow(post.samples.mat)){
#             samp <- post.samples[[i]]
#             p.vec <- samp$p
#             M.vec <- unlist(lapply(1:K, function(k) as.vector(samp$M.array[,,k])))
#             post.samples.mat[i,] <- c(p.vec, M.vec)
#         }
#         ### Column names is first cause names (with prefix p)
#         ### then M (as.vector goes by column)
#         p.names <- paste0("p[", 1:C, "]")
#         M.names <- unlist(lapply(1:K, function(k) {
#             paste0(paste0("M[", rep(1:C, C), ","), rep(1:C, each = C), ",", k, "]")
#         }))
#         cnames <- c(p.names,
#                     M.names)
#         
#         colnames(post.samples.mat) <- cnames
#         return(mcmc(post.samples.mat))
#     })
#     posterior.list <- mcmc.list(posterior.list)
#     posterior.list <- window(posterior.list, start = burnin, thin = thin)
#     return(posterior.list)
# }
