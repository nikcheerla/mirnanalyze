
# This script preprocess the TCGA data. 
# End result is an array of mirna values as well as their class at each stage of the embryonic development tree
# load the libraries first before running the script


#Cancer Datasets to analyze
cancer = c("ACC", "CESC", "CHOL", "DLBC", "ESCA", "KICH", "KIRP", "LGG", 
           "LIHC", "LUSC", "MESO", "OV", "PAAD", "PCPG", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCS", "UVM")


#keep track of mirna data, cancer classification data, number of samples of each type per dataset, and 
#whether data is benign or malignant
mirna_classification <- array(NA)
mirna_detailed_classification = array(NA)

normal_num <- cancer_num <- num_per_dataset <- array(NA)
mat = array(NA, c(length(cancer), 10))
malignancy_classification <- array(NA)
num = 1;
row_na <- array(NA)
col_na <- array(NA)


for(x in cancer) {
  
  #read from directory with cancer name
  name = paste("../cancer_downloaded_data/gdac.broadinstitute.org_",".miRseq_Mature_Preprocess.Level_3.2015040200.0.0/",".miRseq_mature_RPM_log2.txt", sep=x)
  name2 = paste("../cancer_downloaded_data/gdac.broadinstitute.org_",".miRseq_Mature_Preprocess.Level_3.2015040200.0.0/",".miRseq_mature_RPM.txt", sep=x)
  
  #read in mirna data
  if(file.exists(name)) {
    data_t =  read.table(name, header = TRUE);
    isLog = FALSE
  } else {
    data_t = read.table(name2, header = TRUE);
    data_t[, 2:ncol(data_t)] = log2(data_t[, 2:ncol(data_t)]);
    isLog = TRUE
  }
  
  data_t[, 1] = as.character(data_t[, 1]);
  cols = colnames(data_t)
  normalc = 0; cancerc = 0;
  
  #assign classification data based on name of sample
  #01, 02 is cancerous, 11 is normal
  for(p in 2:ncol(data_t)){
    split_list = unlist(strsplit(cols[p], '[.]'))
    
    if(substr(split_list[4], 1, 1) != "0") {
      mirna_classification <- c(mirna_classification, "NORMAL")
      temp <- paste("NORMAL_", x, sep="")
      mirna_detailed_classification = c(mirna_detailed_classification, temp)
      malignancy_classification <- c(malignancy_classification, "N")
      normalc = normalc + 1
    } else {
      mirna_classification <- c(mirna_classification, x)
      mirna_detailed_classification = c(mirna_detailed_classification, x)
      malignancy_classification <- c(malignancy_classification, "P")
      cancerc = cancerc + 1
    }
    normal_num[x] = normalc;
    cancer_num[x] = cancerc;
    num_per_dataset[x] = normalc + cancerc;
  }
  #attach data to larger dataset
  if(x != "ACC") {
    data_t = data_t[, 2:ncol(data_t)]
    mirna <- cbind(mirna, data_t)
  } else {
    mirna <- data_t
  }
  print(x);
  mat[num, 1:7] = summary(apply(data_t[, 2:ncol(data_t)], 1, FUN = "mean"));
  mat[num, 8] =  isLog;
  mat[num, 9] = normalc
  mat[num, 10] = cancerc;
  num = num + 1;
}
rownames(mat) = cancer;
colnames(mat) = c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max. ", "NA's", "log2 missing", "normal samples", "cancer samples")
write.table(mat, "../cancer_stats_summary.txt")

#######################################################################################################################################

##function to impute gene values below a certain threshold of tolerance and remove genes otherwise
runImpute <- function(mat, thresh = 0.25) {
  #preprocessing to remove genes with too many NA values
  bv = array(NA)
  for(i in 1:nrow(mat)){
    data = mat[i, ];
    # keep each row if the number of NA's is less than thresh% of the total number of samples
    bv[i] = (sum(is.na(data)) < thresh*length(data)) 
    if(i %% 100 == 0)
      print(i)
  }
  mat = mat[bv, ]
  #imputation step to fill in leftover NA's with predicted values
  impute.knn(data.matrix(mat[, -1]), rowmax = thresh, colmax=0.5)$data -> mat[, -1]
  mat<-na.omit(mat)
  return(mat)
}


# run the impute function
mirna = runImpute(mirna, 0.25);


#delete the first column that has the miran names (as mirna names are already row names now)
xt = as.character(mirna[, 1])
mirna1 = mirna[,-1]
rownames(mirna1) = xt;
class = mirna_classification[-1]
colnames(mirna1) = gsub("Mature", "", colnames(mirna1))

# create classes at the different devlopment stage levels of Embryonic Development Tree

yt <- class
  
  for(i in 1:length(yt)) {
    yt[i] = par[[yt[i]]];
  }

class_4 <- yt

for(i in 1:length(yt)) {
  yt[i] = par[[yt[i]]];
}

class_3 <- yt

for(i in 1:length(yt)) {
  yt[i] = par[[yt[i]]];
}

class_2 <- yt

  
for(i in 1:length(yt)) {
  yt[i] = par[[yt[i]]];
}

class_1 <- yt

  

######## save data for all 470 mirna #####################################################################################################

#save processed mirna data for each cancer and in total
malignancy_classification <- malignancy_classification[-1]
save(cancer, mirna1, class, class_1, class_2, class_3, class_4, malignancy_classification, file="../mirnaData.RData")
save(cancer, mirna1, class, class_1, class_2, class_3, class_4, malignancy_classification, file="../App/mirnaData.RData")

# save mirna data for each cancer seperately
for(x in cancer){
  tmp = mirna1[, which(class == x)]
  save(tmp, file=paste("../cancer_mirna_RData/", ".RData", sep = x));
}

#write.table(rownames(mirna1), "selectedMiRNA.txt");
data <- t(rbind(mirna1, class))
write.table(data, "../weka_data.csv", row.names = FALSE, sep=",")


######### save with reduced features (selected 132 from Weka feature selection algorithms) ######################################

read.table("features_reduced_1.csv", sep=",") -> features_data
as.character(features_data[,1]) -> features_data[,1]
grep(features_data[,1], pattern = "100 %", value=TRUE) -> features_subset



j=0
feature_index <- NULL
for (i in features_subset) { 
  j= j+1
  print(j)
  arr = strsplit(i, " ", fixed = TRUE)[[1]]
  grep(arr, pattern="hsa") -> index
  hsaname = arr[index]
  print(hsaname)
  grep(rownames(mirna1), pattern=hsaname) -> index
  print(index)
  print(":::")
  feature_index <- c(feature_index, index)
}

mirna1 <- mirna1[feature_index,]

data <- t(rbind(mirna1, class))
write.table(data, "../weka_data_reduced132.csv", row.names = FALSE, sep=",")
save(cancer, mirna1, class, class_1, class_2, class_3, class_4, malignancy_classification, file="../mirnaData_reduced132.RData")
save(cancer, mirna1, class, class_1, class_2, class_3, class_4, malignancy_classification, file="../App/mirnaData_reduced132.RData")


################ save with further reduces features (from weka) 60 features in total ######################################

read.table("features_reduced_2.csv", sep=",") -> features_data
as.character(features_data[,1]) -> features_data[,1]


j=0
feature_index <- NULL
for (i in features_data[,1]) { 
  j= j+1
  print(j)
  arr = strsplit(i, " ", fixed = TRUE)[[1]]
  grep(arr, pattern="hsa") -> index
  hsaname = arr[index]
  print(hsaname)
  grep(rownames(mirna1), pattern=hsaname) -> index
  print(index)
  print(":::")
  feature_index <- c(feature_index, index)
}

feature_index = feature_index[1:60]

mirna1 = mirna1[feature_index, ]

data <- t(rbind(mirna1, class))
write.table(data, "../weka_data_reduced60.csv", row.names = FALSE, sep=",")
save(cancer, mirna1, class, class_1, class_2, class_3, class_4,  malignancy_classification, file="../mirnaData_reduced60.RData")
save(cancer, mirna1, class, class_1, class_2, class_3, class_4, malignancy_classification, file="../App/mirnaData_reduced60.RData")

