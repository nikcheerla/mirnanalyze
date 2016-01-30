
## Do the classification at different stages of the Embryonic Development Tree 


# remove the normal samples
which(class == "NORMAL") -> index
x <- t(mirna1[, -index])
yt <- class[-index]

# Train and validate the models at stage 1, 2, and 3 of the Embryonic Development Tree
matrices = list();
d = 1;
while(length(unique(yt)) > 2){
  
  
  for(i in 1:length(yt)) {
    yt[i] = par[[yt[i]]];
  }
  print("Unique classes: ");
  print(unique(yt));
  
  svmFit = train(x, as.factor(yt),
                 method = "svmRadial", preprocess = c("knnImpute"),
                 tuneLength = 10, trControl = trainControl(method = "cv"))
  
  accuracy = max(svmFit$results["Accuracy"])
  print("Accuracy: ");
  print(accuracy);
  cm = confusionMatrix(svmFit)
  data = cm[[1]]
  
  source('calculate.R')
  sensitivity <- data.frame(sensitivity)
  specificity <- data.frame(specificity)
  rownames(sensitivity) <- rownames(data)
  rownames(specificity) <- rownames(data)
  
  
  matrices[[5*d]] = unique(yt);
  matrices[[5*d + 1]] = svmFit;
  matrices[[5*d + 2]] = accuracy;
  matrices[[5*d + 3]] = sensitivity;
  matrices[[5*d + 4]]  = specificity;
  
  d = d + 1;
}


save(matrices, file="../stem.RData")




