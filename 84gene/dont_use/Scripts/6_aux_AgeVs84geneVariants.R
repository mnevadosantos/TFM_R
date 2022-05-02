###########################################################################################################
##### This script analyzes the variable "age_at_initial_pathologic_diagnosis" in relation to genetic 
##### variants, using regresion with backward step adjust. (Script for running with RScript).
###########################################################################################################
if (!require("pacman")) install.packages("pacman")
p_load(corrplot, ggplot2, dplyr, tidyr, bootStepAIC, MASS, ROCR, car, caret)

setwd("~/Documentos/Olavide/TFM/TFM_R/84gene/Scripts/")
datos_global <- read.table("../Data/datos_global_clean.csv", header = T,
                           stringsAsFactors = F, sep = ";")
dropped_cols <- c(1:13, 15:21, 621)
datos_global <- datos_global[, -dropped_cols]
for (n in 2:600) {
  datos_global[, n] <- as.integer(datos_global[, n])
}
set.seed(123456)
training.ids <- createDataPartition(datos_global$age_at_initial_pathologic_diagnosis, p = 0.75, list = F)
# Regresion analysis:
regmod_all <- lm(formula = age_at_initial_pathologic_diagnosis~., datos_global[training.ids,])

regmod_none <- lm(formula = age_at_initial_pathologic_diagnosis~1, datos_global[training.ids,])

regback <- step(regmod_all, scope = list(lower=regmod_none, upper=regmod_all), direction = "backward")

##### Root mean square error ###################################################
rmse <- sqrt(mean((datos_global[training.ids,]$age_at_initial_pathologic_diagnosis - 
                    regback$fitted.values)^2))

pred <- predict.lm(regback, datos_global[-training.ids,], type = "response")
sqrt(mean((pred - datos_global[-training.ids,]$age_at_initial_pathologic_diagnosis)^2))
pred

save.image("../Data/modelback_learning3.RData")
