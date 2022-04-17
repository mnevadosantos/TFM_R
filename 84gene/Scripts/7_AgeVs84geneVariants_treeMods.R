###########################################################################################################
##### This script analyzes the variable "age_at_initial_pathologic_diagnosis" in relation to genetic 
##### variants, using different machine learning methods.
###########################################################################################################
if (!require("pacman")) install.packages("pacman")
p_load(caret, rpart, rpart.plot)

datos_global <- read.table("../Data/datos_global_clean.csv", header = T,
                           stringsAsFactors = T, sep = ";")
dropped_cols <- c(1:13, 15:21, 621)
datos_global <- datos_global[, -dropped_cols]

training.ids <- createDataPartition(datos_global$age_at_initial_pathologic_diagnosis, p = 0.7, list = F)
modtree <- rpart(age_at_initial_pathologic_diagnosis ~., 
                 data = datos_global[training.ids,],
                 method = "class", 
                 control = rpart.control(minsplit = 20, cp = 0.01))
prp(modtree, type = 2, extra = 103, nn = T, fallen.leaves = T, faclen = 4,
    varlen = 8, shadow.col = "gray")
