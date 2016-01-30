

# This is the blood valodation for the Alzheimer's dataset. First, the common miRNA between the 
# alzheimers dataset and the 60 miRNA were selected. Then, a new SVM model is built using TCGA data and some normal samples
# from this dataset. Then, the rest of the samples were prediected. 68% specificity

# read the GSE data and remove the NA values
data = read.csv("GSE46579_data.csv", header = TRUE, fill = TRUE, row.names =1)
data = log(data);
is.na(data) <- sapply(data, is.infinite)
data <- runImpute(data, 0.2)

# get the common miRNA between the TCGA and the GSE
bv = NULL;
wantedData = data[1, ];
load("../mirnaData_reduced60.RData");
mirnames = NULL;
rownames(data) = paste("", tolower(rownames(data)));

for(mirname in rownames(mirna1)) {
  mrn  = mirname;
  mirname = strsplit(mirname, "|", fixed = TRUE)[[1]][1]
  #mirname = strsplit(mirname, "-", fixed = TRUE)[[1]]
  #mirname = paste(mirname[2:(length(mirname) - 1)], sep = "", collapse = "-")
  mirname = tolower(mirname)
  print(mirname)
  print(grep(rownames(data), pattern= mirname, value=TRUE)[1]);
  print(":::")
  wantedData = rbind(wantedData, data[grep(rownames(data), pattern= mirname)[1], ]);
  rownames(wantedData)[nrow(wantedData)] = mrn;
  
  bv = c(bv, length(grep(rownames(data), pattern= mirname)) == 0);
  
  mirnames = c(mirnames, grep(rownames(data), pattern= mirname, value=TRUE)[1]);
}

wantedData = wantedData[-1, ];


## training by adding some normal samples from GSE to the TCGA samples
x <- mirna1[!bv,]
y <- as.factor(class)
test <- wantedData[!bv,]

# normalization of the GSE samples

x <- t(scale(t(x)))
test <- t((scale(t(test))))

az_normal <- test
save(az_normal, file= "AZdata.RData")


x <- cbind(x, test[, c(1:16, 34:44)])

test <- test[, -c(1:16, 34:44)]
y <- as.factor(c(class, rep("NORMAL", 27)))

x <- t(x)
test <- t(test)


  model <- train(x, y,
                 method = "svmRadial",
                 tuneLength = 20, trControl = trainControl(method = "cv"))
  
  predict(model,test)
  predict(model, test) -> output


##############################################################################


