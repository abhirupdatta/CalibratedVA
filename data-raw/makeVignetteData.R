devtools::load_all()
library(openVA)

adult.raw <- read.csv(getPHMRC_url("adult"))

countries <- ifelse(adult.raw$site %in% c("Dar", "Pemba"), "Tanzania", "Other")
tanzania.data <- adult.raw[countries == "Tanzania",]
train.data <- adult.raw[countries == "Other",]
set.seed(851745)
calibration.indices <- sample(nrow(tanzania.data), 200, replace = F)
calibration.data <- tanzania.data[calibration.indices,]
test.data <- tanzania.data[-calibration.indices,]

set.seed(123)
tariff.train <- codeVA(data = rbind(calibration.data, test.data),
                       data.type = "PHMRC", model = "Tariff",
                       data.train = train.data, causes.train = "gs_text34",
                       phmrc.type = "adult")


set.seed(123)
insilico.train <- codeVA(data = rbind(calibration.data, test.data), data.type = "PHMRC",
                         model = "InSilicoVA",
                         data.train = train.data, causes.train = "gs_text34",
                         phmrc.type = "adult",
                         jump.scale = 0.05, convert.type = "fixed",
                         Nsim=10000, auto.length = FALSE)
### InSilicoVA probabilities
insilico_tanzania <- getIndivProb(insilico.train)
### Tariff top score
tariff_tanzania <- getTopCOD(tariff.train)
usethis::use_data(insilico_tanzania, tariff_tanzania)

### Finally get Tanzania GS COD
gs_cod_tanzania <- tanzania.data$gs_text34
usethis::use_data(gs_cod_tanzania)
