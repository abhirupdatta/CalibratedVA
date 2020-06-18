#' Tariff top COD predictions for adult deaths in Tanzania
#' 
#' A dataset containing top COD predictions for adult deaths in Tanzania,
#' from the PHMRC dataset, using Tariff.
#'
#' @format A data frame with 2023 rows and 2 columns:
#' \describe{
#'   \item{ID}{The ID of the deceased}
#'   \item{cause}{The predicted COD. One of 34 possible causes from gs_text34}
#' }
"tariff_tanzania"

#' InSilicoVA probability of COD predictions for adult deaths in Tanzania
#' 
#' A matrix containing predicted probabilities for each COD for adult deaths in Tanzania,
#' from the PHMRC dataset, using InSilicoVA
#'
#' @format A matrix with 223 rows and 34 columns. Entry i,j is the predicted probability that
#' individual i died from cause j. The rownames are the IDs of the deceased.
"insilico_tanzania"

#' PHMRC gold-standard COD for adult deaths in Tanzania
#' 
#' @format A vector giving the PHMRC gold-standard COD for adult deaths in Tanzania.
#' The order corresonds to the rows in \code{link{tariff_tanzania}}  and \code{link{insilico_tanzania}} 
"gs_cod_tanzania"

#' Map from PHMRC COD to broad COD
#' 
#' Provides a mapping from the 34 causes of death in the PHMRC adult data set (gs_text34)
#' to 5 broad causes of death
#' 
#' @format A dataframe with 34 rows and 2 columns:
#' \describe{
#'   \item{cause}{The PHMRC COD}
#'   \item{braod_cause}{The broad cause that corresponds to the PHMRC COD}
#' }
"phmrc_adult_cause_map"