kfoldCV <- function(n, k){
  # Paras:
  #   x is random objects:
  #     if type == 'distribution',  a n by m matrix for distributional objects; 
  #     if type == 'spd', a n by d by d array (a list of q by q matrices) for SPD matrices objects;
  #     if type == 'sphere', a n by 1 matrix for spherical objects.
  #   k is the number of folds.
  # Returns:
  #   A list with CV splited random objects. 
  
  if ((n %% k) != 0){
    warning("The number of samples cannot be divied by the number of folds")
  } 
  
  n_divided <- floor(n / k)
  
  index_list <- list()
  index_all <- c(1:n)
  index_left <- index_all
  
  for (i in 1:(k-1)){
    index_selected <- sample(index_left, n_divided, replace=FALSE)
    index_list[[i]] <- index_selected
    index_left <- setdiff(index_left, index_selected)
  }
  index_list[[k]] <- index_left
  
  return(index_list)
}
