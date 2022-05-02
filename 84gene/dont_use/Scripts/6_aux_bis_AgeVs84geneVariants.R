###########################################################################################################
##### This script analyzes the variable "age_at_initial_pathologic_diagnosis" in relation to genetic 
##### variants, using regresion with backward and forward step adjust. (Script for running with RScript).
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
training.ids <- createDataPartition(datos_global$age_at_initial_pathologic_diagnosis, p = 0.70, list = F)
# Regresion analysis:
regmod_all <- lm(formula = age_at_initial_pathologic_diagnosis~., datos_global[training.ids,])

regmod_none <- lm(formula = age_at_initial_pathologic_diagnosis~1, datos_global[training.ids,])

regstep <- step(regmod_none, scope = list(lower=regmod_none, upper=regmod_all), direction = "both")

pred <- predict.lm(regstep, datos_global[-training.ids,], type = "response")
sqrt(mean((pred - datos_global[-training.ids,]$age_at_initial_pathologic_diagnosis)^2))

cuadrito <- ggplot(data = datos_global[-training.ids,], aes(x = age_at_initial_pathologic_diagnosis, y = pred)) + geom_point() + geom_smooth(method = "lm", level = 0.99) +
  theme (text = element_text(size=8)) + # Tamaño de fuente del grafico por defecto
  ggtitle("Age at initial pathological diagnosis", subtitle = "Model with both directions step adjust")  + 
  theme (plot.title = element_text(family="NimbusSan",
                                   size=rel(2), #Tamaño relativo de la letra del título
                                   vjust=2, #Justificación vertical, para separarlo del gráfico
                                   face="bold", #Letra negrilla. Otras posibilidades "plain", "italic", "bold" y "bold.italic"
                                   color="red", #Color del texto
                                   hjust = 0.5, # Centrado horizontal
                                   lineheight=1.5),
         panel.background = element_rect(fill = "#D5DEE3")) + 
  labs(x = "Observed values", y = "Predicted values") + 
  theme(axis.title.x = element_text(face="bold", vjust=-0.5, colour="orange", size=rel(1.5))) +
  theme(axis.title.y = element_text(face="bold", vjust=1.5, colour="blue", size=rel(1.5)))


save.image("../Data/modelstep_learning3.RData")

