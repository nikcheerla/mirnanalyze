# calculate sensitivity/specificty

load("../svm_models.RData")


# radial model

output <- array(NA, c(22, 6))

mat <- confusionMatrix(svmradial_model) 
data <- (mat[[1]])
output[, 1:2] <- perf(data)


mat <- confusionMatrix(svmradial_132_model) 
data <- (mat[[1]])
output[, 3:4] <- perf(data)


mat <- confusionMatrix(svmradial_60_model) 
data <- (mat[[1]])
output[, 5:6] <- perf(data)

output <- round(output, 2)
rownames(output) <- rownames(data)

write.table(output, "sen_spec.csv", sep=",")
