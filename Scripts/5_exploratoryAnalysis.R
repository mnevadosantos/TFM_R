

datos_global <- read.csv("../Data/datos_global_clean.csv", sep = ";", header = T, stringsAsFactors = T)
datos_global$X <- NULL
summary(datos_global[, 2:19])
str(datos_global[, 2:19])
