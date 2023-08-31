auc.calculator <- function(df){
  df <- df[order(df$FP),]
  TPR <- df$TP
  FPR <- df$FP
  ## TPR & FPR have to be ordered array
  if(is.unsorted(FPR)) stop('FPR has to be sorted')
  
  # Step 2.1: prepare to calculate AUC
  TPR <- c(0,TPR,1)
  FPR <- c(0,FPR,1)
  len <- length(TPR)
  
  # Step 2.2: calculate AUC ROC
  area <- 0
  for(i in 1:(len-1)){
    delta.x <- FPR[i+1]-FPR[i]
    y.lower <- TPR[i]
    y.upper <- TPR[i+1]
    area <- area + (y.upper + y.lower) * delta.x / 2
  }
  return(list(AUC=as.numeric(area), TPR=TPR[2:(len-1)], FPR=FPR[2:(len-1)]))
}
