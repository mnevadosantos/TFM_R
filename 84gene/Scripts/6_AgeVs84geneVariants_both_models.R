###########################################################################################################
##### This script analyzes the variable "age_at_initial_pathologic_diagnosis" in relation to genetic 
##### variants, using different machine learning methods (regresion with backward step adjust).
###########################################################################################################
if (!require("pacman")) install.packages("pacman")
p_load(corrplot, ggplot2, dplyr, tidyr, bootStepAIC, MASS, ROCR, car, caret)

datos_global <- read.table("../Data/datos_global_clean.csv", header = T,
                           stringsAsFactors = F, sep = ";")
dropped_cols <- c(1:13, 15:21, 621)
datos_global <- datos_global[, -dropped_cols]

# Regresion analysis:
regmod_all <- lm(formula = age_at_initial_pathologic_diagnosis~., datos_global)

regmod_none <- lm(formula = age_at_initial_pathologic_diagnosis~1, datos_global)

regback <- step(regmod_all, scope = list(lower=regmod_none, upper=regmod_all), direction = "backward")
regstep <- step(regmod_none, scope = list(lower=regmod_none, upper=regmod_all), direction = "both")

# Cross-valitation LOOCV method:

library(boot)
cv_error_regback <- cv.glm(data = datos_global, glmfit = regback)
cv_error_regstep <- cv.glm(data = datos_global, glmfit = regstep)

save.image("../Data/both_models.RData")
