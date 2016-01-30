library(rdrop2)

runImpute <- function(mat, thresh = 0.5) {
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

saveData <- function(data) {
  data <- t(data)
  # Create a unique file name
  fileName <- sprintf("%s_data.csv", as.integer(Sys.time()))
  # Write the file to the local system
  write.csv(
    x = data,
    file = file.path(fileName), 
    row.names = FALSE, quote = TRUE
  )
}

loadData <- function() {
  # Read all the files into a list
  files <- list.files(full.names = TRUE)
  data <- lapply(files, read.csv, stringsAsFactors = FALSE) 
  # Concatenate all data together into one data.frame
  data <- do.call(rbind, data)
  data
}


## Get the miRNA matching the ones in the classifier
cleanUp <- function(mirna1, data) { 
  
  bv <- NULL
  mirnames <- NULL
  wantedData <- NULL
  rownames(data) <- data[,1]
  rownames(data) = tolower(rownames(data));
  
  
  for(mirname in rownames(mirna1)) {
    mrn  = mirname;
    mirname = strsplit(mirname, "|", fixed = TRUE)[[1]][1]
    mirname = strsplit(mirname, "-", fixed = TRUE)[[1]]
    mirname = paste(mirname[1:3], sep = "", collapse = "-")
    mirname = tolower(mirname)
    #print(mirname)
    #print(grep(rownames(data), pattern= mirname, value=TRUE)[1]);
    #print(":::")
    wantedData = rbind(wantedData, data[grep(rownames(data), pattern= mirname)[1], ]);
    rownames(wantedData)[nrow(wantedData)] = mrn;
    #bv = c(bv, length(grep(rownames(data), pattern= mirname)) == 0);
  }
  
  print("return")
  return(wantedData)
  
}


###################################

#recommends treatment given mirna data and clinical info
recommend_treatment = function(mirna, cancer_type, info) {
  load("clinicalData.RData");
  load("clinicalAnalysis.RData");

  info[1:2] = toupper(info[1:2]);
  mirnaComplete = mirna;

  cancer_idx = which(cancer_type == cancer)[1]
  
  infoMod = c(which(info[1] == unique(clinicalInfo[, 1]))[1], which(info[2] == unique(clinicalInfo[, 2]))[1], info[3])
  infoMod = as.double(infoMod);
  best = topTreatments[1];
  bestval = 0;
  best_index = NULL;
  for (i in 1: nrow(topCoxTreatments)) {
   
    cox = topCoxTreatments[i,];
    patient = c(mirna, infoMod, cancer_idx, cox)
    
    prob = predict(object=clinFit60, newdata=data.frame(rbind(as.double(patient), as.double(patient))), type = "prob")[1, 2]
    if (prob >= bestval) {
      bestval = prob;
      best = topTreatments[i];
      best_index = i;
    }
  }
  

  return(list(treatment = best, prob = bestval));
}

common_rows <- function (x1, x2) { 
  
  rowx1 <- rownames(x1) 
  rowx2 <- rownames(x2)
  x1_index <- NULL
  x2_index <- NULL
  
  rowcommon <- intersect(rowx1, rowx2)
  
  for (i in 1: length(rowcommon)) { 
    grep(rowcommon[i], rowx1) -> index1
    grep(rowcommon[i], rowx2) -> index2
    
    x1_index <- c(x1_index, index1)
    x2_index <- c(x2_index, index2)
    
    
  }
  
  data = list()
  data = list(x1[x1_index,], x2[x2_index,])
  return(data)
}


recommend_three_treatments = function(mirna, cancer_type, info) {
  
  load("clinicalData.RData");
  load("clinicalAnalysis.RData");
  
  info[1:2] = toupper(info[1:2]);
  mirnaComplete = mirna;
  
  cancer_idx = which(cancer_type == cancer)[1]
  
  infoMod = c(which(info[1] == unique(clinicalInfo[, 1]))[1], which(info[2] == unique(clinicalInfo[, 2]))[1], info[3])
  infoMod = as.double(infoMod);
  
  probarr = NULL;
  for (i in 1: nrow(topCoxTreatments)) {
    cox = topCoxTreatments[i,];
    patient = c(mirna, infoMod, cancer_idx, cox)
    
    prob = predict(object=clinFit60, newdata=data.frame(rbind(as.double(patient), as.double(patient))), type = "prob")[1, 2]
    probarr=c(probarr, prob);
  }
  
  probs = sort(probarr, decreasing = TRUE)[1:3]
  ind = c(which(probarr == probs[1])[1], which(probarr == probs[2])[1], which(probarr == probs[3])[1]);
  
  results_array = NULL;
  results_array = c(probs, topCoxTreatments[ind,], ind);
  
  return(list(results_array = results_array, treatments = topTreatments[ind], probs = probs));
}




######################## semi supervised learning


  




