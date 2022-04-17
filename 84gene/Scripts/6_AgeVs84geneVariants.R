###########################################################################################################
##### This script analyzes the variable "age_at_initial_pathologic_diagnosis" in relation to genetic 
##### variants, using different machine learning methods.
###########################################################################################################
if (!require("pacman")) install.packages("pacman")
p_load(corrplot, ggplot2, dplyr, tidyr, bootStepAIC, MASS, ROCR, car)

datos_global <- read.table("../Data/datos_global_clean.csv", header = T,
                           stringsAsFactors = T, sep = ";")
datos_global$X <- NULL
dropped_cols <- c(1:13, 15:21, 621)
datos_global <- datos_global[, -dropped_cols]
# No missing values:
sum(datos_global[is.na(datos_global$age_at_initial_pathologic_diagnosis),]) 

# Regresion analysis:
regmod_all <- lm(formula = age_at_initial_pathologic_diagnosis~., datos_global)
summary(regmod_all)
regmod_none <- lm(formula = age_at_initial_pathologic_diagnosis~1, datos_global)
summary(regmod_none)
regback <- step(regmod_all, scope = list(lower=regmod_none, upper=regmod_all), direction = "backward")
regstep <- step(regmod_none, scope = list(lower=regmod_none, upper=regmod_all), direction = "both")
summary(regstep)
confint(regstep)
summary(regback)
confint(regback)


# ROC curve:
valorpred <- predict.lm(regback, type = "response")
valpred <- as.data.frame((valorpred))
names(valpred) <- "pred_age_at_init_patdiag"
valpred$early[datos_global$age_at_initial_pathologic_diagnosis <= 70] <- 0
valpred$early[datos_global$age_at_initial_pathologic_diagnosis > 70] <- 1
predicho <- prediction(valpred$pred_age_at_init_patdiag, valpred$early) 
perf <- performance(predicho, "tpr", "fpr")
perf1 <- performance(predicho, measure = "acc")
edad_max <- sapply(perf1@y.values, which.max)
punto_corte <- sapply(perf1@x.values, "[", edad_max)
plot (perf, colorize = T)
abline (0 ,1 , lty =2 , col = " blue " )
auc <- performance(predicho, "auc")
cat ( " AUC = " , auc@y.values [[1]] , " \ n " )
plot(perf1, col = "darkred")
abline(v = punto_corte, lty = 2)
plot(fitted(regback), rstandard(regback)) # Possible collinearity
qqnorm(rstandard(regback))
qqline(rstandard(regback)) # Residuals follow normal distribution.

################################################################################
##### Confusion matrix #####
################################################################################
confusion <- data.frame(datos_global$age_at_initial_pathologic_diagnosis, valpred$pred_age_at_init_patdiag)
confusion$pred[confusion$valpred.pred_age_at_init_patdiag <= 55] <- 'young'
confusion$pred[confusion$valpred.pred_age_at_init_patdiag > 55 &
                 confusion$valpred.pred_age_at_init_patdiag <= 70] <- 'middleage'
confusion$pred[confusion$valpred.pred_age_at_init_patdiag > 70] <- 'old'
confusion$pred <- ordered(confusion$pred, levels = c('young', 'middleage', 'old'))
confusion$real[confusion$datos_global.age_at_initial_pathologic_diagnosis <= 55] <- 'young'
confusion$real[confusion$datos_global.age_at_initial_pathologic_diagnosis > 55 &
                 confusion$datos_global.age_at_initial_pathologic_diagnosis <= 70] <- 'middleage'
confusion$real[confusion$datos_global.age_at_initial_pathologic_diagnosis > 70] <- 'old'
confusion$real <- ordered(confusion$real, levels = c('young', 'middleage', 'old'))
tabla_conf <- table(confusion$real, confusion$pred, dnn = c('Real', 'Prediccion'))
prop.table(tabla_conf)
round(prop.table(tabla_conf, 1)*100, 2)
summary(tabla_conf)
#### Mosaic diagram ####
barplot(tabla_conf, legend = T, xlab = "Prediction",
        col = c("grey", "orange", "darkgreen"))
mosaicplot(tabla_conf, color = "orange", main = 'Model performance')
################################################################################
### Model graphics
################################################################################
ggplot(data = datos_global, aes(x = age_at_initial_pathologic_diagnosis, y = regback$fitted.values)) + geom_point() + geom_smooth(method = "lm", level = 0.99) +
  theme (text = element_text(size=8)) + # Tamaño de fuente del grafico por defecto
  ggtitle("Age at initial pathological diagnosis")  + 
  theme (plot.title = element_text(family="Comic Sans MS",
                                   size=rel(3), #Tamaño relativo de la letra del título
                                   vjust=2, #Justificación vertical, para separarlo del gráfico
                                   face="bold", #Letra negrilla. Otras posibilidades "plain", "italic", "bold" y "bold.italic"
                                   color="red", #Color del texto
                                   hjust = 0.5, # Centrado horizontal
                                   lineheight=1.5),
         panel.background = element_rect(fill = "#BDB76B")) + 
  labs(x = "Observed values", y = "Predicted values") + 
  theme(axis.title.x = element_text(face="bold", vjust=-0.5, colour="orange", size=rel(1.5))) +
  theme(axis.title.y = element_text(face="bold", vjust=1.5, colour="blue", size=rel(1.5)))

# The model includes 318 variants:
length(regback$model)
write.csv2(names(regback$model), file = '../Data/variantes_modelo.csv')

