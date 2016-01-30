######## This script processes multiple files to extract clinical data of the patients
############ then creates a clinicalMirna matrix - a matrix of patients with complete clinical info, treatment
########## and mirna values

minSamplesPerTreatment = 1;
threshold = 40;

## do the analysis only with 60 miRNA
load("../mirnaData_reduced60.RData")
cancer = c("ACC" , "CESC" , "CHOL", "ESCA", "KICH", "KIRP", "LGG", "LIHC", "LUSC", "MESO", "OV", "PAAD", "PCPG", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCS", "UVM");

clinicalClass = array(NA);
clinicalMirna = NULL;
clinicalTreatment = array(NA);
clinicalOutcome = array(NA);
clinicalName = array(NA);
clinicalInfo = c("F", "FREE", 3000);
therapies = NULL;

for(x in cancer){
  print(x);
  
  treatment_filename = paste("../cancer_clinical_data/", "/nationwidechildrens.org_clinical_drug_", ".txt", sep=tolower(x));
  treatment = read.table(treatment_filename, header = TRUE, sep = "\t", fill = TRUE);
  
  result_filename = paste("../cancer_clinical_data/", "/nationwidechildrens.org_clinical_follow_up_v4.0_", ".txt", sep=tolower(x));
  result = read.table(result_filename, header = TRUE, sep = "\t", fill = TRUE);
  
  radiation_filename = paste("../cancer_clinical_data/", "/nationwidechildrens.org_clinical_radiation_", ".txt", sep=x);
  radiation = read.table(radiation_filename, header = TRUE, sep = "\t", fill = TRUE);
  radiationType = radiation[3:nrow(radiation), ]
  
  options(warn=-1)
  patient_data_filename = paste("../cancer_clinical_data/", "/nationwidechildrens.org_clinical_patient_", ".txt", sep=x);
  patient_data = read.table(patient_data_filename, header = TRUE, sep = "\t", fill = TRUE);
  options(warn=0)
  
  patient_names = treatment[3:nrow(treatment), 2]
  patient_names = patient_names[!duplicated(patient_names)]
  treatment[, 2] = gsub("-", ".", treatment[, 2])
  radiation[, 2] = gsub("-", ".", radiation[, 2])
  patient_data[, 2] = gsub("-", ".", patient_data[, 2])
  result[, 2] = gsub("-", ".", result[, 2])
  patient_names = gsub("-", ".", patient_names)
  d = 0;
  data_t = array(NA);
  sample_idx = array(NA);
  treatment_type = array(NA);
  outcome = array(NA);
  for(i in which(class == x)){
    for(j in 1:length(patient_names) ){
      name = patient_names[j]
      if(grepl(name, colnames(mirna1)[i])){
        clinicalName = c(clinicalName, name);
         data_t = cbind(data_t, mirna1[, i]);
         sample_idx = c(sample_idx, i)
         
         arr = array(NA);
         for(k in which(treatment[, 2] == name)){
            type <- "";
            type = paste(type, gsub("-", "", gsub("\\s", "", tolower(treatment[k, 9]))), sep = "");
            type = paste(type, gsub("-", "", gsub("\\s", "", tolower(treatment[k, 7]))), sep = ", ");
            arr = c(arr, type);
         }
         for(k in which(radiation[, 2] == name)){
           type <- "";
           type = paste(type, "radiationtherapy", sep = "");
           type = paste(type, gsub("-", "", gsub("\\s", "", tolower(radiation[k, 6]))), sep = ", ");
           arr = c(arr, type);
         }
    
         b = FALSE;
         for(k in which(result[, 2] == name)){
           d = d + 1;
           txt = as.character(result[k, "treatment_outcome_first_course"]);
           outcome = c(outcome, txt);
           #print(txt)
           if(identical(txt, character(0))) {
             b = FALSE;
           }
           else {
             b = TRUE;
           }
           break; ## SY why is there a break here? 
         }
         if(!b){
           outcome = c(outcome, "[No outcome recorded]");
           #print(length(outcome));
           #print("outcome increased");
           #print(outcome)
         }
        
        b = FALSE;
        for(k in which(patient_data[, 2] == name)){
          d = d + 1;
          txt = as.character(patient_data[k, "gender"]);
          race = as.character(patient_data[k, "race"]);
          eth = as.character(patient_data[k, "ethnicity"]);
          if (eth != "NOT HISPANIC OR LATINO" && race == "WHITE") {
              race = "LATINO";
          }
          if(grepl(race, pattern = "[", fixed = TRUE)) {
            race = "OTHER";
          }
          age = as.numeric(as.character(patient_data[k, "age_at_initial_pathologic_diagnosis"]));
          if(identical(age, numeric(0)) || is.na(age)) {
            age = 60;
          }
          age = round(age, -1)
          clinicalInfo = rbind(clinicalInfo, c(txt, race, age));
          
         
          if(identical(txt, character(0))) {
            b = FALSE;
          }
          else {
            b = TRUE;
          }
          break;
        }
        
        if(!b){
          clinicalInfo = rbind(clinicalInfo, c("OTHER", race, age));
          
        }
         arr = unique(arr);
        
         arr = sort(arr);
         str = paste(arr, sep = "", collapse = "; ");
         treatment_type = c(treatment_type, str);
      }
      
    }
  }
  outcome = outcome[-1];
  arr = arr[-1];
  rownames(data_t) = NULL;
  data_t = data_t[, -1]
  sample_idx = sample_idx[-1]
  clinicalClass= c(clinicalClass, class[sample_idx]);
  clinicalMirna = cbind(clinicalMirna, data_t);
  clinicalTreatment = c(clinicalTreatment, treatment_type)
  clinicalOutcome = c(clinicalOutcome, outcome)
}

clinicalOutcome = gsub("Complete Remission/Response", "Remission", clinicalOutcome)
clinicalOutcome = gsub("Partial Remission/Response", "Remission", clinicalOutcome) # could be remission or disease
clinicalOutcome = gsub("Progressive Disease", "Disease", clinicalOutcome)
clinicalOutcome = gsub("Stable Disease", "Disease", clinicalOutcome)
clinicalOutcome = gsub("No Measureable Tumor or Tumor Markers", "Remission", clinicalOutcome)
clinicalOutcome = gsub("Normalization of Tumor Markers, but Residual Tumor Mass", "Disease", clinicalOutcome)

clinicalClass = clinicalClass[-1]
rownames(clinicalMirna) = NULL;
clinicalTreatment = clinicalTreatment[-1]
clinicalTreatment = clinicalTreatment[-1]
clinicalOutcome = clinicalOutcome[-1]
clinicalName = clinicalName[-1]
clinicalInfo = clinicalInfo[-1, ]
clinicalTreatment = na.omit(clinicalTreatment)


#calculate top treatment options and match other treatments to those
topTreatments = NULL;
for(x in clinicalTreatment){
  if(length(which(clinicalTreatment == x)) >= minSamplesPerTreatment && !grepl("\\[", x)) {
      topTreatments = c(topTreatments, x);
  }
}
  
topTreatments = unique(topTreatments)

cost = list();
cost["insertions"] = 5;
cost["deletions"] = 3;
cost["substitutions"] = 6;

bm = NULL;
for(i in 1:length(clinicalTreatment)){
  if(! (clinicalTreatment[i] %in% topTreatments)){
    best = topTreatments[1];
     for(x in topTreatments){
       if(adist(clinicalTreatment[i], best) > adist(clinicalTreatment[i], x)) {
         x = best;
       }
     }
    bm[i] = (adist(clinicalTreatment[i], best) <= threshold);
    clinicalTreatment[i] = best;
  } else {
    bm[i] = TRUE;
  }
  if(grepl("\\[", clinicalOutcome[i])){
    bm[i] = FALSE;
  }
  
  if(grepl("\\[", clinicalInfo[i])){
    bm[i] = FALSE;
  }
}



for(x in clinicalTreatment) {
  for(t in strsplit(x, "; "))
    therapies = c(therapies, t);
}
therapies = unique(therapies);

clinicalClass = clinicalClass[bm];
clinicalMirna = t(clinicalMirna);
clinicalMirna = clinicalMirna[bm, ];
rownames(clinicalMirna) = NULL;
clinicalTreatment = clinicalTreatment[bm];
clinicalOutcome = clinicalOutcome[bm];
clinicalName = clinicalName[bm];
clinicalInfo = clinicalInfo[bm, ];

colnames(clinicalMirna) <- rownames(mirna1)



clinicalTreatmentVectorized = array(NA);
for(i in 1:length(clinicalTreatment)){
  clinicalTreats = array(NA);
  for(j in 1:length(therapies)) {
    clinicalTreats = c(clinicalTreats, grepl(therapies[j], clinicalTreatment[i]));
  }
  clinicalTreatmentVectorized = rbind(clinicalTreatmentVectorized, clinicalTreats);
}

clinicalTreatmentVectorized[-1, ] -> clinicalTreatmentVectorized;
clinicalTreatmentVectorized[, -1] -> clinicalTreatmentVectorized;
rownames(clinicalTreatmentVectorized) = NULL;


uniqueTreatments = unique(clinicalTreatment)
topTreatments = unique(clinicalTreatment);

numPatients4Treatment <- NULL
numCancers4Treatment <- NULL
inum = 0 
too_few_index1 = NULL

for (i in topTreatments) { 
  
  inum = inum+1
  which(clinicalTreatment == i) -> index
  numPatients4Treatment[inum] = length(index) 
  
  if(numPatients4Treatment[inum] < 5) too_few_index1 <- c(too_few_index1, index)
  
  numCancers4Treatment[inum] <- length(unique(class[index]))
  
}

uniqueCancers = unique(clinicalClass)
numPatients4Cancer <- NULL
numtreatments4Cancer <- NULL

inum = 0
for (i in uniqueCancers) { 
  inum = inum+1
  
  which(clinicalClass == i) -> index
  numPatients4Cancer[inum] = length(index)
  numtreatments4Cancer[inum ] <- length(unique(clinicalTreatment[index]))
}

#### find the cancers with less than 5 patients
too_few_index <- NULL
inum = 0
for (i in uniqueCancers) { 
  inum = inum+1
  if(numPatients4Cancer[inum] < 5) 
    too_few_index <- c(too_few_index,  which(clinicalClass == i))
}


freq = NULL
for(i in 1:length(topTreatments)) {
  freq[i] = length(which(topTreatments[i] == clinicalTreatment));
}



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


####### remove the cancer types with less than 5 patients

too_few_index <- c(too_few_index, too_few_index1)
clinicalMirna = clinicalMirna[-too_few_index, ]
clinicalTreatment = clinicalTreatment[-too_few_index]
clinicalTreatmentVectorized = clinicalTreatmentVectorized[-too_few_index, ]
clinicalOutcome = clinicalOutcome[-too_few_index]
clinicalName = clinicalName[-too_few_index]
clinicalClass = clinicalClass[-too_few_index]
clinicalInfo = clinicalInfo[-too_few_index, ]

uniqueTreatments = unique(clinicalTreatment)
topTreatments = unique(clinicalTreatment);

numPatients4Treatment <- NULL
numCancers4Treatment <- NULL
inum = 0 

for (i in topTreatments) { 
  
  inum = inum+1
  which(clinicalTreatment == i) -> index
  numPatients4Treatment[inum] = length(index) 
  numCancers4Treatment[inum] <- length(unique(class[index]))
  
}

uniqueCancers = unique(clinicalClass)
numPatients4Cancer <- NULL
numtreatments4Cancer <- NULL

inum = 0
for (i in uniqueCancers) { 
  inum = inum+1
  
  which(clinicalClass == i) -> index
  numPatients4Cancer[inum] = length(index)
  numtreatments4Cancer[inum ] <- length(unique(clinicalTreatment[index]))
}


save(clinicalInfo, 
     clinicalClass, 
     clinicalName, 
     clinicalOutcome,
     clinicalTreatmentVectorized, 
     clinicalTreatment, 
     clinicalMirna, 
     cancer_name, 
     topTreatments,
     freq, file="../clinicalData.RData")

