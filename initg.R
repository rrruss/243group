initG <- function(f, upperB = Inf, lowerB = -Inf) {
  if(upperB == Inf || lowerB == -Inf) {
    #If either bound is infinite, use 5 or -5 as starting point for upper or lower bounds
    u <- (max(lowerB, -6)+1)
    v <- (min(upperB, 6)-1)    
  } else {
    #If neither bound is infinite, divide interval into thirds
    u <- (2 * lowerB + upperB) / 3
    v <- (2 * upperB + lowerB) / 3
  }
  
  #Use the fPrime module to evaluate slopes at u and v
  uSlope <- fPrime(u)
  vSlope <- fPrime(v)
  
  #If f'(u) <= 0, slide u towards lowerB until f'(u) > 0
  while(uslope <= 0) {
    u <- (2 * u + max(3*u, lowerB)) / 3
    uSlope <- fPrime(u)
  }
  #If f'(v) >= 0, slide v towards upperB until f'(v) < 0
  while(vslope >= 0) {
    v <- (2 * v + min(3*v, upperB)) / 3
    vSlope <- fPrime(v)
  }
  uVal <- fVal(u)
  vVal <- fVal(v)
  
  #Find the intersection point of the tangent lines
  w <- (logf(v) - logf(u) - v*vSlope + u*uSlope)/(uSlope - vSlope)
  
  #Create g
  g <- matrix(c(lowerB, w, u, uSlope, uVal, w, upperB, v, vSlope, vVal), nrow=2, ncol=5)
  colnames(g) = c('start', 'end', 'intersect', 'm', 'b')
  return(g)
}