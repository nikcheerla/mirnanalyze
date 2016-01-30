
treats = array(NA);
for(i in 1:length(clinicalClass)) {
  print(i);
  treats = c(treats, recommend_treatment(clinicalMirna[1, ], clinicalClass[1], clinicalInfo[1, ])$prob);
}

