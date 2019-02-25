## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  cache = TRUE
)

## ----load_pkgs-----------------------------------------------------------
library(openVA)
library(CalibratedVA)

## ----convert_data--------------------------------------------------------
child.raw <- read.csv(getPHMRC_url("child"))
child.clean <- ConvertData.phmrc(child.raw, phmrc.type = "child")$output

## ----split_data----------------------------------------------------------
countries <- ifelse(child.raw$site %in% c("Dar", "Pemba"), "Tanzania", "Other")
tanzania.data <- child.clean[countries == "Tanzania",]
train.data <- child.clean[countries == "Other",]
set.seed(851745)
calibration.indices <- sample(nrow(tanzania.data), 100, replace = F)
calibration.data <- tanzania.data[calibration.indices,]
test.data <- tanzania.data[-calibration.indices,]

## ----tariff_train--------------------------------------------------------
set.seed(123)
tariff.train <- codeVA(data = rbind(calibration.data, test.data),
                     data.type = "customize", model = "Tariff",
                     data.train = train.data, causes.train = "Cause")

## ----insilico_train------------------------------------------------------
set.seed(123)
insilico.train <- codeVA(data = rbind(calibration.data, test.data),
                         data.type = "customize", model = "InSilicoVA",
                         data.train = train.data, causes.train = "Cause",
                         jump.scale = 0.05, Nsim=5000, auto.length = FALSE)

## ----top_3_cod-----------------------------------------------------------
top.cod <- names(sort(table(tanzania.data$Cause), decreasing = TRUE))
top3.cod <- top.cod[top.cod != "14"][1:3]
change.cause <- function(cause) {
    cause <- as.character(cause)
    cause[!(cause %in% top3.cod)] <- "99"
    return(cause)
}
tariff.train.cod <- change.cause(getTopCOD(tariff.train)[,2])
insilico.train.cod <- change.cause(getTopCOD(insilico.train)[,2])
test.changedcod <- change.cause(test.data$Cause)
calibration.changedcod <- change.cause(calibration.data$Cause)

## ----separate_cod--------------------------------------------------------
tariff.train.cod.test <- tariff.train.cod[-(1:100)]
tariff.train.cod.calib <- tariff.train.cod[1:100]
insilico.train.cod.test <- insilico.train.cod[-(1:100)]
insilico.train.cod.calib <- insilico.train.cod[1:100]

## ----hyperparams---------------------------------------------------------
causes <- as.character(sort(unique(test.changedcod)))
epsilon <- .001
alpha <- 5
beta <- .5
tau <- .5
tau.vec <- rep(tau, length(causes))
delta <- 1
gamma.init <- 1
ndraws <- 50E3
nchains <- 3

## ----tariff_calibva, cache = TRUE, message = FALSE-----------------------
set.seed(123)
calibva.seeds <- sample(1e6, nchains, replace = F)
tariff.calibva <- calibva.sampler(test.cod = tariff.train.cod.test,
                                  calib.cod = tariff.train.cod.calib,
                                  calib.truth = calibration.changedcod, causes = causes,
                                  epsilon = epsilon, alpha=alpha, beta=beta,
                                  tau.vec=tau.vec, delta=delta,
                                  gamma.init=gamma.init, ndraws = ndraws,
                                  nchains = nchains,
                                  init.seeds = calibva.seeds)


## ----gamma_acceptance_rates_tariff---------------------------------------
gamma_acceptance_rates(tariff.calibva)

## ----tariff_thin_burnin--------------------------------------------------
tariff.calibva <- window(tariff.calibva, start = 10e3, thin = 10)

## ----insilico_calibva, cache = TRUE, mesage = FALSE----------------------
insilico.calibva <- calibva.sampler(test.cod = insilico.train.cod.test,
                                    calib.cod = insilico.train.cod.calib,
                                    calib.truth = calibration.changedcod, causes = causes,
                                    epsilon = epsilon, alpha=alpha, beta=beta,
                                    tau.vec=tau.vec, delta=delta,
                                    gamma.init=gamma.init, ndraws = ndraws,
                                    nchains = nchains,
                                    init.seeds = calibva.seeds)
insilico.calibva <- window(insilico.calibva, start = 10E3, thin = 10)

## ----extract_tariff_csmf-------------------------------------------------
library(tidyverse)
library(ggmcmc)
tariff.calibva.csmf.samples <- calibvaCSMFPosteriorSamples(tariff.calibva, causes = causes)

## ----tariff_trace_plot---------------------------------------------------
ggplot(tariff.calibva.csmf.samples, aes(x = Iteration, y = value, color = factor(Chain))) +
  geom_line() +
  facet_wrap(~cause)

## ----tariff_gelman_rubin-------------------------------------------------
tariff.calibva.rhat.df <- rename(tariff.calibva.csmf.samples, Parameter = cause)
ggs_Rhat(tariff.calibva.rhat.df)

## ----tariff_posterior_dens-----------------------------------------------
tariff.calibva.csmf.mean <- calibvaCSMFPosteriorSummary(tariff.calibva.csmf.samples)
tariff.csmf.test <- getRawCSMF(tariff.train.cod.test, causes = causes)
ggplot(tariff.calibva.csmf.samples, aes(x = value)) +
  stat_density(geom="line", position = "identity") +
  facet_wrap(~cause, ncol = 1, scales = "free_y") +
  geom_vline(data = tariff.calibva.csmf.mean, aes(xintercept = mean, colour = "CalibratedVA")) +
  geom_vline(data = tariff.csmf.test, aes(xintercept = csmf, colour = "Tariff")) +
  scale_color_discrete(name = "Method") +
  xlim(0, 1)

## ----tariff_m_matrix-----------------------------------------------------
tariff.posterior.M <- mMatrixPosteriorSummary(tariff.calibva, causes = causes)
tariff.posterior.M$m_posterior_mean
tariff.posterior.M$m_posterior_var

## ----tariff_m_matrix_empirical-------------------------------------------
T.emp <-rawMisclassificationMatrix(calib.cod = insilico.train.cod.calib,
                                   calib.truth = calibration.changedcod,
                                   causes = causes)
M.emp <- normalizedMisclassificationMatrix(T.emp)
M.emp

## ----tariff_map, message = FALSE-----------------------------------------
tariff.map <- calibva.map(test.cod = tariff.train.cod.test,
                          calib.cod = tariff.train.cod.calib,
                          calib.truth = calibration.changedcod, causes = causes,
                          epsilon = epsilon, alpha=alpha, beta=beta,
                          delta=delta)
tariff.map$p
tariff.map$M

## ----calibva_ensemble, cache = TRUE, message = FALSE---------------------
test.cod.mat <- matrix(c(tariff.train.cod.test,  insilico.train.cod.test), ncol = 2)
calib.cod.mat <- matrix(c(tariff.train.cod.calib,  insilico.train.cod.calib), ncol = 2)
ensemble.calibva <- calibva.ensemble.lite.sampler(test.cod.mat = test.cod.mat,
                                                  calib.cod.mat = calib.cod.mat,
                                                  calib.truth = calibration.changedcod,
                                                  causes = causes,
                                                  epsilon = epsilon, alpha=alpha, beta=beta,
                                                  tau.vec=tau.vec, delta=delta,
                                                  gamma.init=gamma.init, ndraws = ndraws,
                                                  nchains = nchains,
                                                  init.seeds = calibva.seeds)
ensemble.calibva <- window(ensemble.calibva, start = 10E3, thin = 10)

## ----obtain_csmf---------------------------------------------------------
tariff.calibva.csmf.samples <- calibvaCSMFPosteriorSamples(tariff.calibva, causes = causes)
tariff.calibva.csmf.mean <- calibvaCSMFPosteriorSummary(tariff.calibva.csmf.samples)

insilico.calibva.csmf.samples <- calibvaCSMFPosteriorSamples(insilico.calibva, causes = causes)
insilico.calibva.csmf.mean <- calibvaCSMFPosteriorSummary(insilico.calibva.csmf.samples)

ensemble.calibva.csmf.samples <- calibvaCSMFPosteriorSamples(ensemble.calibva, causes = causes)
ensemble.calibva.csmf.mean <- calibvaCSMFPosteriorSummary(ensemble.calibva.csmf.samples)

tariff.train.csmf <- getRawCSMF(tariff.train.cod.test, causes)
insilico.train.csmf <- getRawCSMF(insilico.train.cod.test, causes)

## ----csmf_accuracy-------------------------------------------------------
methods <- c("tariff_calibva",
             "insilico_calibva",
             "ensemble_calibva",
             "tariff_train",
             "insilico_train")
ptrue <- sapply(causes, function(c) mean(change.cause(tanzania.data$Cause) == c))
csmf.acc.df <- data.frame(csmf.acc = c(getCSMF_accuracy(tariff.calibva.csmf.mean$mean, ptrue),
                                       getCSMF_accuracy(insilico.calibva.csmf.mean$mean, ptrue),
                                       getCSMF_accuracy(ensemble.calibva.csmf.mean$mean, ptrue),
                                       getCSMF_accuracy(tariff.train.csmf$csmf, ptrue),
                                       getCSMF_accuracy(insilico.train.csmf$csmf, ptrue)),
                          method = methods)
csmf.acc.df

