library(modeest)
library(raster)
library(moments)

datos_global <- read.csv("../Data/datos_global_clean.csv", sep = ";", header = T, stringsAsFactors = T)
datos_global$X <- NULL
clin_index <- c(2:18, 620)
summary(datos_global[, clin_index])
str(datos_global[, clin_index])

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

summary(WHO_group)
mfv(WHO_group)
var(WHO_group)
sd(WHO_group)
table(WHO_group)
skewness(WHO_group)
kurtosis(WHO_group)
hist(WHO_group, 
     col = "darkgreen",
     xlab = "WHO's grade group", 
     ylab = "Frequency", 
     main = "Histogram of WHO's grade groups")

summary(laterality)
summary(history_other_malignancy)
summary(race)
levels(race)
datos_global$race[datos_global$race == "[Not Evaluated]"] <- NA

summary(age_at_initial_pathologic_diagnosis)
nrow(datos_global[is.na(age_at_initial_pathologic_diagnosis), ])
mfv(age_at_initial_pathologic_diagnosis)
var(age_at_initial_pathologic_diagnosis)
sd(age_at_initial_pathologic_diagnosis)
skewness(age_at_initial_pathologic_diagnosis)
kurtosis(age_at_initial_pathologic_diagnosis)
hist(age_at_initial_pathologic_diagnosis, 
     breaks = 6,
     col = "darkgreen",
     xlab = "Age at initial pathologic diagosis", 
     ylab = "Frequency", 
     main = "Histogram of age at diagnosis")

summary(clinical_stage)

summary(pathologic_T)
barplot(table(pathologic_T), 
        col = "darkgreen", 
        xlab = "Pathologic T",
        ylab = "Frequency", 
        main = "Pathologic T")

summary(clinical_N)
summary(clinical_M)

# Initial study will include the age at initial histologic diagnosis, the Gleason score,
# the WHO's grade group and the pathologic T.
