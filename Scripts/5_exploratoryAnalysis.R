library(modeest)
library(raster)
library(moments)

datos_global <- read.csv("../Data/datos_global_clean.csv", sep = ";", header = T, stringsAsFactors = T)
datos_global$X <- NULL
summary(datos_global[, 2:19])
str(datos_global[, 2:19])

attach(datos_global)
summary(gleason_score)
mfv(gleason_score)
var(gleason_score)
sd(gleason_score)
table(gleason_score)
skewness(gleason_score)
kurtosis(gleason_score)
hist(gleason_score, 
     breaks = 5, 
     col = "darkgreen",
     xlab = "Gleason", 
     ylab = "Frequency", 
     main = "Histogram of Gleason's score")

summary(who)
mfv(who)
var(who)
sd(who)
table(who)
skewness(who)
kurtosis(who)
hist(who, 
     breaks = 5, 
     col = "darkgreen",
     xlab = "WHO's grade group", 
     ylab = "Frequency", 
     main = "Histogram of WHO's grade groups")
