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
  #browser()
  library(numDeriv)
  #Use the fPrime module to evaluate slopes at u and v
  logf <- function(x) { log(f(x)) }
  uSlope <- grad(logf, x=u)
  vSlope <- grad(logf, x=v)
  
  #If f'(u) <= 0, slide u towards lowerB until f'(u) > 0
  while(uSlope <= 0) {
    u <- (2 * u + max(3*u, lowerB)) / 3
    uSlope <- grad(f, x=u)
  }
  #If f'(v) >= 0, slide v towards upperB until f'(v) < 0
  while(vSlope >= 0) {
    v <- (2 * v + min(3*v, upperB)) / 3
    vSlope <- grad(f, x=v)
  }
  uVal <- f(u)
  vVal <- f(v)
  
  #Find the intersection point of the tangent lines
  w <- (log(f(v)) - log(f(u)) - v*vSlope + u*uSlope)/(uSlope - vSlope)
  
  #Create g
  g <- rbind(c(lowerB, w, u, uSlope, uVal), c(w, upperB, v, vSlope, vVal))
  colnames(g) = c('start', 'end', 'intersect', 'm', 'b')
  return(g)
}