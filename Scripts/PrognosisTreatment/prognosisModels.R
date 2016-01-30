
############ This script has the prognosis prediction models with and without the miRNA data

# clinical prgnosis using only 60 miRNA values

load("../clinicalData.RData");
source("../helper_functions.R");


cancer_index = array(NA)
for(i in 1:length(clinicalTreatment))
  cancer_index[i] = which(clinicalClass[i] == cancer)

treatment_index = array(NA)
for(i in 1:length(clinicalTreatment))
  treatment_index[i] = which(clinicalTreatment[i] == topTreatments)

dists = adist(topTreatments, partial = FALSE)
fit3D <- cmdscale(dists,eig=TRUE, k=3)

index3 = NULL
for(i in 1:length(clinicalTreatment))
  index3 = rbind(index3, fit3D$points[which(clinicalTreatment[i] == topTreatments)[1], ])

clinicalInfo2 = clinicalInfo
for(i in 1:length(clinicalTreatment)) {
  clinicalInfo2[i, 2] = which(clinicalInfo[i, 2] == unique(clinicalInfo[, 2]))[1];
  clinicalInfo2[i, 1] = which(clinicalInfo[i, 1] == unique(clinicalInfo[, 1]))[1];
}

clinicalInfo3 = matrix(0, length(clinicalTreatment), 3)
clinicalInfo3[, 1] = as.numeric(clinicalInfo2[, 1])
clinicalInfo3[, 2] = as.numeric(clinicalInfo2[, 2])
clinicalInfo3[, 3] = as.numeric(clinicalInfo2[, 3])


#################### model with clinical, treatment and miRNA values

x = cbind(clinicalMirna, clinicalInfo3, cancer_index, index3);
colnames(x) = c(colnames(clinicalMirna), "gender", "Ethnicity", "age", "cancer", 
                                  "treatment_vector1", "treatment_vector2", "treatment_vector3")
x = data.frame(x)
y = data.matrix(clinicalOutcome);

data <- data.frame(cbind(x, y))
newData <- SMOTE(y ~ ., data, perc.over = 100,perc.under=200)
x <- newData[, 1:(ncol(data)-1)]
y <- newData[, ncol(data)]

clinFit60 <- train(x, as.factor(y),
                 method = "svmRadial", 
                 metric = "ROC", tuneLength = 5, trControl = trainControl(method = "cv", classProbs = TRUE), prob.model = TRUE)


#################### model with only treatment option and no MIRNA


x1 = cbind(clinicalInfo3, cancer_index, index3);
colnames(x1) = c("gender", "Ethnicity", "age", "cancer", 
                "treatment_vector1", "treatment_vector2", "treatment_vector3")

x1 = data.frame(x1)
y = data.matrix(clinicalOutcome);


data <- data.frame(cbind(x1, y))
newData <- SMOTE(y ~ ., data, perc.over = 100,perc.under=200)
x1 <- newData[, 1:(ncol(data)-1)]
y <- newData[, ncol(data)]



clinFit60_nomirna <- train(x1, as.factor(y),
                   method = "svmRadial", 
                   metric = "ROC", tuneLength = 5, trControl = trainControl(method = "cv", classProbs = TRUE), prob.model = TRUE)




########### save files

write.table(cbind(x, y), "../weka_clinical60.csv", sep=",", row.names=FALSE)


############# predict values and write it for matlab analysis of ROC curve

output <- predict(clinFit60, x)
output <- as.numeric(as.factor(output))

output_nomirna <- predict(clinFit60_nomirna, x1)
output_nomirna <- as.numeric(as.factor(output_nomirna))

write.table(cbind(y, output, output_nomirna), "../matlab.csv", sep=",")

####################### Do the sensitivity spacifiity

perf_data <- array(NA, c(2,4))
data <- confusionMatrix(clinFit60)[[1]]
perf_data[, 1:2] <- perf(data)

accuracy = (data[1,1] + data[2,2])/sum(data)

data <- confusionMatrix(clinFit60_nomirna)[[1]]
perf_data[, 3:4] <- perf(data)

accuracy_nomirna = (data[1,1] + data[2,2])/sum(data)


topCoxTreatments = fit3D$points;

save(clinFit60, topCoxTreatments, topTreatments, index3, clinicalInfo3, clinicalInfo2, clinicalInfo, cancer, fit3D, file="../clinicalAnalysis.RData")


