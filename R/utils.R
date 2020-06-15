#' @title organizes the posterior samples for the CSMF parameters from
#' \code{calibva.sampler} into a tibble
#' @param calibva.samples a list returned from \code{calibva.sampler}
#' @param causes the cause vector input to CalibVA
#' @return a tibble with each row representing a draw from the posterior sample
#' of the CSMF for a given cause 
#' 
#' @importFrom ggmcmc ggs
#' 
#' @export
calibvaCSMFPosteriorSamples <- function(calibva.samples, causes) {
    C <- length(causes)
    P <- data.frame(
        Parameter = paste0("p[", 1:C, "]"),
        Label = causes
    )
    p_tibble <- ggmcmc::ggs(calibva.samples, family = "p", par_labels = P) %>%
        rename(cause = Parameter)
    return(p_tibble)
}

#' @title collect CSMF posterior summaries from the CalibVA sampler
#' 
#' @param calibvaCSMFPosteriorSamples a tibble returned from \code{calibvaCSMFPosteriorSamples}
#' @param percentile.L the lower percentile for a credible interval. Default .025
#' @param percentile.U the upper percentile for a credible interval. Default .975
#' 
#' @return a tibble with the posterior means, and confidence intervals of the CSMF and the names of each cause, for each of the posterior samples that are obtained
#' 
#' @import dplyr
#' 
#' @export
calibvaCSMFPosteriorSummary <- function(calibvaCSMFPosteriorSamples, percentile.L = .025, percentile.U = .975) {
    p_summary <-
        calibvaCSMFPosteriorSamples %>%
        group_by(cause, ParameterOriginal) %>%
        summarize(mean = mean(value), var = var(value),
                  ci.L = quantile(value, percentile.L), ci.U = quantile(value, percentile.U)) 
    return(p_summary)
} 

#' @title Obtain raw CSMF estimates from verbal autopsy guesses
#' @param test.cod will be a vector of length N, with each entry as the estimated
#' COD (as a character)for indiv. i 
#' @param causes is a character vector with the names of the causes you are interested i
#' @return a data frame with the causes as one column and the CSMF for each cause as 
#' the second column
#' 
#' @export
getRawCSMF <- function(test.cod, causes) {
    csmf <- sapply(causes, function(c) mean(test.cod == c))
    return(data.frame(cause = causes, csmf = csmf))
}


#' @title obtain the raw misclassification matrix using gold standard and VA COD
#' @param calib.cod will be a vector of length N, with each entry as the estimated
#' COD (as a character) for indiv. i in the calibration set
#' @param calib.truth is a character vector with the true COD for each subject in the
#' calibration set
#' @return an integer matrix T where (i,j)th entry is the number of subjects
#' in the calibration set with gold standard COD i and VA COD j
#' @export
rawMisclassificationMatrix <- function(calib.cod, calib.truth, causes) {
    C <- length(causes)
    T.mat <- matrix(NA, nrow = C, ncol = C, dimnames = list(causes, causes))
    for(i in 1:C){
        for(j in 1:C){
          T.mat[i,j] <- sum(calib.truth == causes[i] & calib.cod == causes[j])
        }
    }
    return(T.mat)
}

#' @title normalize the raw misclassification matrix by row
#' @param T.mat an integer matrix produced by \code{rawMisclassificationMatrix}
#' 
#' @return A numeric matrix, where the rows are the conditional misclassification
#' probabilities
#'
#' @export
normalizedMisclassificationMatrix <- function(T.mat) {
    C <- ncol(T.mat)
    M.mat <- matrix(NA, nrow = nrow(T.mat), ncol = ncol(T.mat))
    colnames(M.mat) <- colnames(T.mat)
    rownames(M.mat) <- rownames(T.mat)
    for(i in 1:C) {
        rowSum <- sum(T.mat[i,])
        M.mat[i,] <- T.mat[i,] / rowSum
    }
    return(M.mat)
}

#' @title Get the posterior mean and variance for the M matrix
#' @param calibva.samples a mcmc.list object with samples from the CalibVA sampler
#' @param causes the cause vector input to CalibVA
#' @param output.format either "matrix" or "tibble" 
#' 
#' @return Either a tibble with the entry wise posterior means or a list of matrices
#' with the posterior means stored as m_posterior_mean and the posterior variance
#' stored as m_posterior_variance 
#' 
#' @import ggmcmc
#' @import dplyr
#' 
#' @export
mMatrixPosteriorSummary <- function(calibva.samples, causes, output.format = c("matrix", "tibble")[1]) {
    if(!(output.format %in% c("matrix", "tibble"))) {
        stop("output.format must be either a matrix or tibble")
    }
    C <- length(causes)
    P <- data.frame(
        Parameter = paste0(paste0("M[", rep(1:C, C), ","), rep(1:C, each = C), "]"),
        Label = paste0(rep(causes, C), ",", rep(causes, each = C))
    )
    m_tibble <- ggmcmc::ggs(calibva.samples, family = "M", par_labels = P)
    m_summary <-
        m_tibble %>%
        group_by(Parameter, ParameterOriginal) %>%
        summarize(mean = mean(value), var = var(value))
    if(output.format == "tibble") {
        return(m_summary)
    } else {
        mean.M <- matrix(NA, nrow = C, ncol = C, dimnames = list(causes, causes))
        var.M <- matrix(NA, nrow = C, ncol = C, dimnames = list(causes, causes))
        for(i in 1:nrow(mean.M)){
            for(j in 1:ncol(mean.M)) {
                m.entry <- filter(m_summary, ParameterOriginal == paste0("M[", i, ",", j, "]"))
                mean.M[i,j] <- m.entry$mean
                var.M[i,j] <- m.entry$var
            }
        }
        return(list(m_posterior_mean = mean.M, m_posterior_var = var.M))
    }
}

acceptance.rate <- function(x) {
    naccept <- 0
    for(i in 2:length(x)) {
        if(x[i] != x[i-1]) {
            naccept <- naccept + 1
        }
    }
    return(naccept / (length(x) + 1))
}

#' @title computes the acceptance rates for each gamma parameter from a
#' \code{calibva.sampler} object for each chain
#' @param calibva.samples a list returned from \code{calibva.sampler}
#' @return a tibble giving the acceptance rate for each gamma parameter in a chain
#' 
#' @import ggmcmc
#' @import dplyr
#' 
#' @export
gamma_acceptance_rates <- function(calibva.samples) {
    gamma.tibble <- ggs(calibva.samples, family = "gamma")
    acceptance_rates <-
        gamma.tibble %>%
        group_by(Chain, Parameter) %>%
        summarize(rate = acceptance.rate(value))
    return(acceptance_rates)
}
