library(openVA)
library(CalibratedVA)

child.raw <- read.csv(getPHMRC_url("child"))
child.clean <- ConvertData.phmrc(child.raw, phmrc.type = "child")$output

countries <- ifelse(child.raw$site %in% c("Dar", "Pemba"), "Tanzania", "Other")
tanzania.data <- child.clean[countries == "Tanzania",]
train.data <- child.clean[countries == "Other",]
set.seed(851745)
calibration.indices <- sample(nrow(tanzania.data), 100, replace = F)
calibration.data <- tanzania.data[calibration.indices,]
test.data <- tanzania.data[-calibration.indices,]

set.seed(123)
tariff.train <- codeVA(data = rbind(calibration.data, test.data),
                       data.type = "customize", model = "Tariff",
                       data.train = train.data, causes.train = "Cause")


set.seed(123)
insilico.train <- codeVA(data = rbind(calibration.data, test.data),
                         data.type = "customize", model = "InSilicoVA",
                         data.train = train.data, causes.train = "Cause",
                         jump.scale = 0.05, Nsim=5000, auto.length = FALSE)

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