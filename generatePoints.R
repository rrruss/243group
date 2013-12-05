### using for loop to generate points
generatePoints <- function(N, glist){
  X <- sampleg(N)
  # sample N points from upper bound of g function
  U <- runif(N)
  # generate N points from uniform (0,1) distribution
  sampleX <- NULL
  # a vector to store X's that can be accepted
  for (i in 1:N){
    index_Upper <- which(glist$Upper[ ,'start'] <= X[i] & glist$Upper[ ,'end'] > X[i])
    index_Lower <- which(glist$Lower[ ,'start'] <= X[i] & glist$Lower[ ,'end'] > X[i])
    # in fact index_Upper and index_Lower have some relation
    upperX <- glist$Upper[index_Upper, 'm'] * (X[i] - glist$Upper[index_Upper, 'intersect']) + glist$Upper[index_Upper, 'b']
    # value of upper bound at X[i]
    lowerX <- glist$Lower[index_Lower, 'm'] * (X[i] - glist$Lower[index_Lower, 'start']) + glist$Upper[index_Lower, 'b']
    # value of lower bound at X[i]
    if (U[i] < exp(lowerX - upperX)){
      sampleX <- c(sampleX, X[i])
    }
    else{
      glist <- updateG(X[i], glist, f)
      if (U[i] <- glist$fx / exp(upperX)){
        sampleX <- c(sampleX, X[i])
        break
      }
      else{
        break
      } 
    }    
 }
 return(sampleX)   
}

### another approach, in a vectorize fashion
# the vectorize approach may be faster than the previous one, since it aviods using loop. 
# But if the rate of rejection is high, then evaluting upperX and lowerX at length of N may be a waste of time. 
# Because only the several beginning X's will be accpeted, most of X's will be rejected.
# At the beginning, N is large, the rejection rate may be high, then this approach may be less efficient.
# So we need to test using time of these two approachs and pick the faster one.
generatePoints <- function(N, glist){
  X <- sampleg(N)
  # sample N points from upperg function
  U <- runif(N)
  # generate N points from uniform (0,1) distribution
  sampleX <- NULL
  # a vector to store X's that can be accepted
  FindInxUpper <- function(x){ which(glist$Upper[ ,'start'] <= x & glist$Upper[ ,'end'] > x) }
  FindInxLower <- function(x){ which(glist$Lower[ ,'start'] <= x & glist$Lower[ ,'end'] > x) }
  index_Upper <- sapply(X, FindInxUpper)
  # index to show which part of upper bound function should be used to evaluate X
  index_Lower <- sapply(X, FindInxLower)
  # index to show which part of lower bound function should be used to evaluate X
  upperX <- glist$Upper[index_Upper, 'm'] * (X[i] - glist$Upper[index_Upper, 'intersect']) + glist$Upper[index_Upper, 'b']
  # value of upper bound at X, it is a vector
  lowerX <- glist$Lower[index_Lower, 'm'] * (X[i] - glist$Lower[index_Lower, 'start']) + glist$Upper[index_Lower, 'b']
  # value of lower bound at X, it is a vector
  index_test <- U < exp(lowerX - upperX)
  index <- which( index_test == 0)[1] # find the fisrt position in X to be rejected or to evaluate g(x)
  updateG(X[index], glist, f)
  if (U[index] <- f(X[index]) / exp(upperX[index])){
    sampleX <- X[1 : index] 
  }
  else{
    sampleX <- X[1:(index-1)]
  }               
  return(sampleX)   
}
