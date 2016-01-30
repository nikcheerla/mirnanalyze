
# Treatment Recommendation Tool and it's analysis on the test data

# Recommend Treatment Function

recommend_treatment = function(mirna, cancer_type, info) {
  
  
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
  
  return(list(results_array = results_array, treatments = topTreatments[ind], probabilities = probs));
}


actual_array = array(NA, c(nrow(clinicalMirna), 4))
results_array_1 = array(NA, c(nrow(clinicalMirna), 4))
results_array_2 = array(NA, c(nrow(clinicalMirna), 4))
results_array_3 = array(NA, c(nrow(clinicalMirna), 4))
dist_array_1 = array(NA)
dist_array_2 = array(NA)
dist_array_3 = array(NA)
final_results = array(NA, c(nrow(clinicalMirna), 16))
final_results_short = array(NA, c(nrow(clinicalMirna), 14))

clinicalTreatmentIdx = array(NA)
for ( i in 1:length(clinicalTreatment)) { 
   clinicalTreatmentIdx[i] = which(clinicalTreatment[i]== topTreatments)
}
  

# find the distance from the top three recommended treatments to actual treatment. 
# compare that with the random treatments to actual treatments

for ( i in 1: nrow(clinicalMirna)) { 
  
  mirna = clinicalMirna[i,]
  cancer_type = clinicalClass[i]
  info = clinicalInfo[i,]
  print(i)
  results = recommend_treatment(mirna, cancer_type, info)

  results_array_1[i, ] = results$results_array[c(1, 4, 7, 10)]
  results_array_2[i, ] = results$results_array[c(2, 5, 8, 11)]
  results_array_3[i, ] = results$results_array[c(3, 6, 9, 12)]
  
  

  info[1:2] = toupper(info[1:2]);
  mirnaComplete = mirna;
  cancer_idx = which(cancer_type == cancer)[1]
  infoMod = c(which(info[1] == unique(clinicalInfo[, 1]))[1], which(info[2] == unique(clinicalInfo[, 2]))[1], info[3])
  infoMod = as.double(infoMod);
  cox = index3[i,]
  patient = c(mirna, infoMod, cancer_idx, cox)
  
  actual_array[i,1] = predict(object=clinFit60, newdata=data.frame(rbind(as.double(patient), as.double(patient))), type = "prob")[1, 2];
  actual_array[i, 2:4] = topCoxTreatments[which(clinicalTreatment[i] == topTreatments)[1], ];

  dist_1 = dist(rbind(results_array_1[i,2:4], actual_array[i, 2:4]))
  dist_2 = dist(rbind(results_array_2[i,2:4], actual_array[i,2:4]))
  dist_3 = dist(rbind(results_array_3[i,2:4], actual_array[i,2:4]))
  
  
  random_1 = sample(length(topTreatments), 1)
  random_2 = sample(length(topTreatments), 1)
  random_3 = sample(length(topTreatments), 1)
  
  random_1_dist = dist(rbind(random_1, actual_array[i,2:4]))
  random_2_dist = dist(rbind(random_2, actual_array[i,2:4]))
  random_3_dist = dist(rbind(random_3, actual_array[i,2:4]))
  
  
  final_results[i, ] = c(actual_array[i,1], 
                    results$results_array[1:3], 
                    dist_1, 
                    dist_2, 
                    dist_3,
                    min(dist_1, dist_2, dist_3),
                    random_1_dist, 
                    random_2_dist,
                    random_3_dist,
                    min(random_1_dist, random_2_dist, random_3_dist),
                    results$results_array[13:15], clinicalTreatment[i]
                    )
  
  final_results_short[i, ] = c(actual_array[i,1], 
                         results$results_array[1:3], 
                         results$results_array[13:15], 
                         clinicalTreatmentIdx[i], 
                         topTreatments[results$results_array[13]], 
                         topTreatments[results$results_array[14]], 
                         topTreatments[results$results_array[15]], 
                         clinicalTreatment[i], 
                         clinicalOutcome[i], 
                         clinicalClass[i]
                         )
  
}


# Find the re-occurrence cases in test samples that have treatments with >0.70 remission

idx = array(NA)
for (i in 1:nrow(final_results_short)) { 
  if ((final_results_short[i, 1] < 0.3) &&
      ((final_results_short[i, 2] > 0.7) || (final_results_short[i, 3] > 0.7) || (final_results_short[i, 4] > 0.7))) {
 
    idx <- c(idx, i);
  
  }
  
}
idx <- idx[-1]


remission_results = final_results_short[idx, ]
  
final_results = final_results[which(clinicalOutcome == "Remission"),]
