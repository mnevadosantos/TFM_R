###########################################################################################################
##### This script analyzes the variable "age_at_initial_pathologic_diagnosis" in relation to genetic 
##### variants, using different machine learning methods.
###########################################################################################################
if (!require("pacman")) install.packages("pacman")
p_load(caret, rpart, rpart.plot, randomForest)

datos_global <- read.table("../Data/datos_global_clean.csv", header = T,
                           stringsAsFactors = T, sep = ";")
dropped_cols <- c(1:13, 15:21, 621)
datos_global <- datos_global[, -dropped_cols]
set.seed(2022)

training.ids <- createDataPartition(datos_global$age_at_initial_pathologic_diagnosis, p = 0.75, list = F)
modtree <- rpart(age_at_initial_pathologic_diagnosis ~., 
                 data = datos_global[training.ids,],
                 method = "class", 
                 control = rpart.control(minsplit = 20, cp = 0.01))
prp(modtree, type = 2, extra = 103, nn = T, fallen.leaves = T, faclen = 4,
    varlen = 8, shadow.col = "gray")
pred <- predict(modtree, datos_global[-training.ids,], type = "class")
table(datos_global[-training.ids, "age_at_initial_pathologic_diagnosis"], pred,
      dnn = c("Actual", "Predicted"))
valor_pred <- as.data.frame(pred)
datos_test <- datos_global[-training.ids, ]
valor_pred$early[datos_test$age_at_initial_pathologic_diagnosis <= 70] <- 0
valor_pred$early[datos_test$age_at_initial_pathologic_diagnosis > 70] <- 1
predecir <- prediction(as.numeric(valor_pred$pred), valor_pred$early)
perf <- performance(predecir, "tpr", "fpr")
################################################################################
########## RANDOM FOREST ##########
################################################################################

modforest <- randomForest(x = datos_global[training.ids, 2:600], 
                          y = datos_global[training.ids, 1], 
                          ntree = 500, 
                          keep.forest = T)
pred2 <- predict(modforest, datos_global[-training.ids, ], type = "class")
valor_pred <- as.data.frame(pred2)
