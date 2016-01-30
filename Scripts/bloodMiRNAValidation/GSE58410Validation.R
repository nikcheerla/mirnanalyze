# Validation using cell-free miRNA samples from GSE58410 data set

# read the GSE data and remove the NA values

read.table("GSE58410_data.csv", sep=",") -> data
rownames(data) <- data[,1]
data <- data[, -1]

bv = NULL;
wantedData = data[1, ];
load("../mirnaData_reduced60.RData");
mirnames = NULL;
rowname_data <- tolower(rownames(data))


# Find the common miRNA with TCGA

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

##################### Normalize, train and testing


x <- mirna1[!bv,];
y <- as.factor(class_detailed_4)
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
               tuneLength = 5, trControl = ctrl)

predict(model,test)
predict(model, test) -> output



