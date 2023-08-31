###############
## matrix power
####################
matpower = function(a, alpha){
  a = round((a + t(a))/2,7); tmp = eigen(a)
  return(tmp$vectors%*%diag((tmp$values)^alpha)%*%t(tmp$vectors))
}

#############
## operator norm
################
onorm = function(a) 
  return(eigen(round((a+t(a))/2,8))$values[1])

################
## ## symmetry a asymmetric matrix
########################
sym = function(a) 
  return(round((a+t(a))/2,9))

#####################
## expand grids using x to d dimension
#########################
expand.grids <- function(x,d) {
  expand.grid(replicate(d, x, simplify=FALSE))
}

#########################
#  Moore-Penrose type power           
#  Taking power ignoring 0 eigenvalues;           
#  ignoring criterion=ignore 
#############################
mppower = function(matrix,power,ignore){
  matrix = sym(matrix)
  p = nrow(matrix)
  eig = eigen(matrix)
  eval = eig$values
  evec = eig$vectors
  m = length(eval[abs(eval)>ignore])
  D = diag(eval[1:m]^power)
  if(dim(D)[1]==0){
    return(matrix(0,nrow = p, ncol = p))
  }
  else{tmp = evec[,1:m]%*%diag(eval[1:m]^power)%*%t(evec[,1:m])
    return(tmp)
  }
}


#####################
## ridge power inverse
######################
riginv = function(matrix, eps){
  p = nrow(matrix)
  if(dim(matrix)[1]==0){
    return(matrix(0,nrow = p, ncol = p))
  }
  else{tmp = solve(matrix + eps*diag(p))
  return(tmp)
  }
}


## Screen the sample covariance matrix S

## Input:
## S: Sample Covariance Matrix: Mp times Mp
## lambda : tunning parameter

## Output:
## A: output adjacencey matrix (1: connected edges, 0: disconnected edges)
## K: number of connected components of the graphs
## index[[i]]: id for each connected component
screenS <-function(S, p, lambda){
  A <-matrix(0, nrow=p, ncol=p)
  M<-dim(S)[1]/p
  
  for(i in 1:(p-1)){
    index.i <-((i-1)*M+1):(i*M)
    for(j in (i+1):p){
      index.j <-((j-1)*M+1):(j*M)
      
      if(norm(S[index.i, index.j], type="F") > lambda){
        A[i,j] <-1
      }      
    }
  }
  A <-A + t(A) +diag(rep(1,p))
  g <-graph.adjacency(A)
  cout<-clusters(g)
  id <-cout$membership
  K <-cout$no
  index <-list()
  for(i in 1:K){
    index[[i]] <-which(id==i)
  }
  return(list(A=A, K=K, index=index))
}
