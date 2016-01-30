# Validation using cell-free miRNA samples from GSE67075 data set


# Read data and preprocess to add the class 
source("../helper_functions.R")
data = read.table("GSE67075_data.txt", header = TRUE, fill = TRUE)

for (i in 2:ncol(data)) { 
  data[, i] <- as.numeric(as.character(data[,i]))
}

# 8 control
# 8 polyp
# 16 advanced adenoma
# 8 early
# 8 late stage

colnames(data) <- paste(
  c(rep("NORMAL",9),rep("NORMAL",8),rep("AD", 16), rep("ECRC",8), rep("LCRC", 8)), 
  c(1:8,1:8, 1:16, 1:8, 1:8),
  sep="")

class_data <- c(rep("NORMAL",8),rep("NORMAL",8),rep("AD", 16), rep("ECRC",8), rep("LCRC", 8))


# Load the miRNA data and find the common miRNA between TCGA dataset and the GSE67075 dataset

bv = NULL;
wantedData = data[1, ];
load("../mirnaData_reduced60.RData");
mirnames = NULL;
rowname_data <- array(NA)

for (i in 1:nrow(data)) { 
  arr = strsplit(as.character(data[i,1]), "-", fixed=TRUE)[[1]]
  arr = paste(arr[1:(length(arr) - 1)], sep="", collapse ="-")
  rowname_data[i] = tolower(arr);
  }
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
wantedData = wantedData[, -1];

##################### Train the classifier with TCGA and validate with GSE67075 data #############

# add a few cell-miRNA samples from Alzheimer's dataset to the TCGA classifier

load("../AZdata.RData")

x <- mirna1[!bv,]
y <- as.factor(malignancy_classification)
test <- wantedData[!bv,]

# Normalization of the dataset

x <- t(scale(t(x)))
test <- t((scale(t(test))))

common_elements <- common_rows(x, az_normal)
common_elements_1 <- common_rows(test, az_normal)

x <- common_elements[[1]]
test <- common_elements_1[[1]]
az_normal <- common_elements[[2]]

x <- cbind(x, az_normal)
y <- as.factor(c(malignancy_classification, rep("N", 70)))


x <- t(x)
test <- t(test)

ctrl <- trainControl(method = "cv",
                     ## new option here:
                     sampling = "smote")
model <- train(x, y,
               method = "svmRadial",
               tuneLength = 30, trControl = ctrl)

predict(model,test)
predict(model, test) -> output




