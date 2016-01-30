
# In this script, various classifiers were built and their performance is compared. 
# the best classifier is picked and it is further tuned
# The selected classifier was also trained with reduced features to get performance
# all the models were saved at the end

load("../mirnaData.RData");

x = t(mirna1)
y = class


svmlinear_model <- train(x, as.factor(y),
                   method = "svmLinear",
                   trControl = trainControl(method = "cv"))

svmradial_model <- train(x, as.factor(y),
                   method = "svmRadial",
                   trControl = trainControl(method = "cv"))




ctree_model <- train(x, as.factor(y),
                   method = "ctree",
                   trControl = trainControl(method = "cv"))


knn_model <- train(x, as.factor(y),
                method = "knn",
                tuneLength = 1, trControl = trainControl(method = "cv"))



lda_model <- train(x, as.factor(y),
                      method = "lda",
                      tuneLength = 1, trControl = trainControl(method = "cv"))


nb_model <- train(x, as.factor(y),
                   method = "nb",
                   trControl = trainControl(method = "cv"))


allmodels <- list(svmLinear = svmlinear_model, 
                  svmRadial = svmradial_model,
                  ctree = ctree_model, 
                  knn = knn_model, 
                  lda = lda_model, 
                  logisticModelTree = LogisticModelTree_model, 
                  naivebayes = nb_model)

results <- resamples(allmodels)

# Plot the various classifier performance
bwplot(results)

################### SVM is picked. Tune SVM parameters to get the best performace. 


load("../mirnaData.RData");

x = t(mirna1)
y = class


svmradial_model <- train(x, as.factor(y),
                         method = "svmRadial",
                         tuneLength = 20, 
                         trControl = trainControl(method = "cv", repeats = 10))



######################## SVM with 132 features ######################################



load("../mirnaData_reduced132.RData");

x = t(mirna1)
y = class


svmradial_132_model <- train(x, as.factor(y),
                         method = "svmRadial",
                         tuneLength = 20, 
                         trControl = trainControl(method = "cv", repeats = 10))



################################# SVM with 60 features ################################


load("../mirnaData_reduced60.RData");

x = t(mirna1)
y = class


svmradial_60_model <- train(x, as.factor(y),
                         method = "svmRadial",
                         tuneLength = 20, 
                         trControl = trainControl(method = "repeatedcv", repeats =10))


################### SVM with center/scale preprocess ######################################



load("../mirnaData.RData");

x = t(mirna1)
y = class


svmradial_prep_model <- train(x, as.factor(y),
                         method = "svmRadial",
                         preProcess = c("center", "scale"),
                         tuneLength = 20, 
                         trControl = trainControl(method = "cv", repeats = 10))


############################# save all the models ###########################################


save(svmradial_model, svmradial_prep_model, svmradial_60_model, svmradial_132_model, file = "../svm_models.RData")


