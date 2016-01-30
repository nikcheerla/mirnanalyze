
# validation with GSE59565 samples

load("../mirnaData_reduced60.RData")


for ( i in 1:23) { 
  k = 1439480 + i
  file_name = paste("GSE59565/GSM", k , "_sample", i, ".txt", sep="");
  data_t = read.table(file_name, sep = " ", fill = TRUE);
  if (i == 1) { 
    ndata <- data_t
    ndata[,2] <- as.numeric(ndata[,2])
  } else { 
    data_t[,2] <- as.numeric(data_t[,2])
    ndata <- cbind(ndata, as.matrix(data_t[,2])) 
    } 
}

rownames(ndata) <- ndata[,1]
ndata <- data.matrix(ndata)
ndata <- ndata[, -1]
colnames(ndata) <- paste("col", 1:23, sep = "")
ndata <- log2(ndata)
is.na(ndata) <- sapply(ndata, is.infinite)

# impute function

ndata <- runImpute(ndata, 0.2)
data = ndata

bv = NULL;
wantedData = data[1, ];
mirnames = NULL;
rowname_data = rownames(data)

for(mirname in rownames(mirna1)) {
  mrn  = mirname;
  mirname = strsplit(mirname, "|", fixed = TRUE)[[1]][1]
  mirname = tolower(mirname)
  print(mirname)
  print(grep(rowname_data, pattern= mirname, value=TRUE)[1]);
  print(":::")
  wantedData = rbind(wantedData, data[grep(rowname_data, pattern= mirname)[1], ]);
  rownames(wantedData)[nrow(wantedData)] = mrn;
  
  bv = c(bv, length(grep(rowname_data, pattern= mirname)) == 0);
  
  mirnames = c(mirnames, grep(rowname_data, pattern= mirname, value=TRUE)[1]);
}

wantedData = wantedData[-1, ];


#################### train and test

x <- mirna1[!bv,]
y <- as.factor(class)
test <- wantedData[!bv,]

load("../AZdata.RData")

x <- t(scale(t(x)))
test <- t((scale(t(test))))

load("../AZdata.RData")

common_elements <- common_rows(test, az_normal)
common_elements_x <- common_rows(x, az_normal)

x <- common_elements_x[[1]]
az_normal <- common_elements[[2]]
test <- common_elements[[1]]

x<- cbind(x, az_normal)
y<- as.factor(c(malignancy_classification, rep("N", 70)))


x <- t(x)
test <- t(test)


ctrl <- trainControl(method = "cv",
                     ## new option here:
                     sampling = "smote")
model <- train(x, y,
               method = "svmLinear",
               tuneLength = 20, trControl = ctrl)

predict(model,test)
predict(model, test) -> output


