################## cross validation for choosing the tuning parameters in block GLasso #########

CV4tuning <- function(x, p, lambda, n.fold){
  #Para
  n = nrow(x)
  M = ncol(x) / p
  cv_index = kfoldCV(n, n.fold)
  
  nloglik = rep(0, n.fold)
  for(i in 1:n.fold){
    index = cv_index[[i]]
    x_train = x[-index,]
    x_test = x[index,]
    S_train = cor(x_train) + 10^-3*diag(M*p)
    S_test = cor(x_test)
    result = bglasso(S = S_train, lambda = lambda, n = n, p = p)
    print(result$iter.Theta)
    Theta_pred = result$Theta
    
    nloglik[i] = -log(det(Theta_pred)) + sum(diag(S_test%*%Theta_pred))
  }
  
  avg_nloglik = mean(nloglik)
  
  return(avg_nloglik)
}
