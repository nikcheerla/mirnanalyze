
library("caTools")
library("e1071")
library("Hmisc")
library("impute")
library("glmnet")
library("preprocessCore")
library("shiny")
library("caret")
library("mlbench")
library("DMwR")
library("rdrop2")

source("../development_tree.R")

cancer = c("ACC", "CESC", "CHOL", "DLBC", "ESCA", "KICH", "KIRP", "LGG", 
           "LIHC", "LUSC", "MESO", "OV", "PAAD", "PCPG", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCS", "UVM")

cancer_name = list()

cancer_name[["ACC"]] = "Adrenocortical carcinoma"        
cancer_name[["CESC"]] = "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma"  
cancer_name[["CHOL"]] = "Cholangiocarcinoma"      
cancer_name[["ESCA"]] = "Esophageal Carcinoma"
cancer_name[["KICH"]] = "Kidney Chromophobe"        
cancer_name[["KIRP"]] = "Kidney Renal Papillary cell carcinoma"        
cancer_name[["LGG"]] = "Brain Lower Grade Glioma"
cancer_name[["LIHC"]] = "Liver Hepatocellular Carcinoma"
cancer_name[["LUSC"]] = "Lung Squamous Cell Carcinoma"
cancer_name[["MESO"]] = "Mesothelioma"         
cancer_name[["OV"]] = "Ovarian Serous Cystadenocarcinoma" 
cancer_name[["PAAD"]] = "Pancreatic Adenocarcinoma" 
cancer_name[["PCPG"]] = "Pheochromocytoma and Paraganglioma"       
cancer_name[["SKCM"]] = "Skin Cutaneous Melanoma"
cancer_name[["STAD"]] = "Stomach Adenocarcinoma"      
cancer_name[["TGCT"]] = "Testicular Germ Cell Tumors"
cancer_name[["THCA"]] = "Thyroid Carcinoma"
cancer_name[["THYM"]] = "Thymoma"
cancer_name[["UCS"]] = "Uterine Carcinosarcoma"
cancer_name[["UVM"]] = "Uveal Melanoma"         
cancer_name[["DLBC"]] = "Diffuse Large B-cell Lymphoma"
cancer_name[["NORMAL"]] = "NORMAL"



### some helpful functions

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
 
 # sensitivity and specificity 

 perf <- function (data) { 
  
 sensitivity <- array(NA)
 specificity <- array(NA)

 for ( k  in 1:nrow(data)) { 
   
   sum_tn = 0;
   sum_fp = 0;
   sum_fn = 0;
   for ( i in 1: nrow(data) ) { 
     for ( j in 1: nrow(data)) { 
       
       if ( ( k != j) && (k != i) ) 
         sum_tn = sum_tn + data[i,j]
       
       if (( j == k) && (i != k))
         sum_fn = sum_fn + data[i,j]
       
       if (( i ==k) && ( j != k))
         sum_fp = sum_fp + data[i,j]
       
     }
   }
   
   specificity[k] = sum_tn/(sum_tn + sum_fp)
   sensitivity[k] = data[k,k]/(sum_fn + data[k,k])
   
  

 }
 output <- t( data.frame(rbind(sensitivity, specificity)))
 rownames(output) <- rownames(data)
 return(output)
 }
 
 
 
 
 