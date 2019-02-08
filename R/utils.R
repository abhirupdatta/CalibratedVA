#' @title collect CSMF posterior summaries from the CalibVA sampler
#' 
#' @param calibva.samples a list returned from \code{calibva.sampler}
#' @param causes the cause vector input to CalibVA
#' @param percentile.L the lower percentile for a credible interval. Default .025
#' @param percentile.U the upper percentile for a credible interval. Default .975
#' 
#' @return a tibble with the posterior means, and confidence intervals of the CSMF and the names of each cause
#' for each of the draws that are obtained
#' 
#' @export
calibvaCSMFPosteriorSummary <- function(calibva.samples, causes, percentile.L = .025, percentile.U = .975) {
    C <- length(causes)
    P <- data.frame(
        Parameter = paste0("p[", 1:C, "]"),
        Label = causes
    )
    p_tibble <- ggmcmc::ggs(calibva.samples, family = "p", par_labels = P)
    p_summary <-
        p_tibble %>%
        group_by(Parameter, ParameterOriginal) %>%
        summarize(mean = mean(value), var = var(value),
                  ci.L = quantile(value, percentile.L), ci.U = quantile(value, percentile.U)) %>%
        rename(cause = Parameter)
    return(p_summary)
} 

#' @title Obtain raw CSMF estimates from verbal autopsy guesses
#' @param test.cod will be a vector of length N, with each entry as the estimated
#' COD (as a character)for indiv. i 
#' @param causes is a character vector with the names of the causes you are interested i
#' 
#' @export
rawCSMF <- function(test.cod, causes) {
    csmf <- sapply(causes, function(c) mean(test.cod == c))
    return(data.frame(cause = causes, csmf = csmf))
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