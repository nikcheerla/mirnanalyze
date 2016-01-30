
# sensitivity and specificity 

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
        sum_fp = sum_fp + data[i,j]
      
      if (( i ==k) && ( j != k))
        sum_fn = sum_fn + data[i,j]
      
    }
  }
  
  specificity[k] = sum_tn/(sum_tn + sum_fp)
  sensitivity[k] = data[k,k]/(sum_fn + data[k,k])
  
}