# Validation using cell-free miRNA samples from GSE46355 data set

# read the GSE data and remove the NA values

data1 = read.table("GSE46355_data_1.txt", fill = TRUE, row.names =1)
data2 = read.table("GSE46355_data_2.txt", fill = TRUE, row.names =1)
  
data1 <- data1[-1, -1]
data2 <- data2[, -1]

data <- cbind(data1, data2)


for (i in 1:ncol(data)) { 
  data[, i] <- as.numeric(as.character(data[,i]))
}

is.na(data) <- sapply(data, is.infinite)
data <- runImpute(data, 0.2)

# get the common miRNA between the TCGA and the GSE
bv = NULL;
wantedData = data[1, ];
load("../mirnaData_reduced60.RData");
mirnames = NULL;

for (i in 1:nrow(data)) {
  rowname = rownames(data)[i]
  rowname = strsplit(rowname, "-", fixed = TRUE)[[1]]
  rowname = paste(rowname[1:(length(rowname) - 1)], sep = "", collapse = "-")
  rownames(data)[i] <- rowname
  }

  
rownames(data) = tolower(rownames(data))

for(mirname in rownames(mirna1)) {
  mrn  = mirname;
  mirname = strsplit(mirname, "|", fixed = TRUE)[[1]][1]
  #mirname = strsplit(mirname, "-", fixed = TRUE)[[1]]
  #mirname = paste(mirname[1:(length(mirname) - 1)], sep = "", collapse = "-")
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


## Normalize the GSE data. Train and TCGA model with common miRNA and validate with GSE data

x <- mirna1[!bv,]
y <- as.factor(class)
test <- wantedData[!bv,]


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
                 tuneLength = 5, trControl = ctrl)
  
predict(model,test)
predict(model, test) -> output


#########################################################################

 
save(test, file  = "AZdata.RData")
  
  

