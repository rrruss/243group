initG <- function(f, upperB = Inf, lowerB = -Inf, initVals) {
  if (!is.null(initVals) & !(all(initVals > lowerB & initVals < upperB))) {
    stop("Initial values specified have to be strictly between the upper and lower bounds specified.")
  }
  if (!is.null(initVals)) {
    initVals <- sort(initVals)
    u <- initVals[1]
    v <- initVals[2]
  } else if(upperB == Inf || lowerB == -Inf) {
    #If either bound is infinite, use 5 or -5 as starting point for upper or lower bounds
    u <- (max(lowerB, -6)+1)
    v <- (min(upperB, 6)-1)    
  } else {
    #If neither bound is infinite, divide interval into thirds
    u <- (2 * lowerB + upperB) / 3
    v <- (2 * upperB + lowerB) / 3
  }
  # browser()
  #Use the fPrime module to evaluate slopes at u and v
  logf <- function(x) { log(f(x)) }
  uSlope <- grad(logf, x=u, method='simple')
  vSlope <- grad(logf, x=v, method='simple')
  
  #If f'(u) <= 0, slide u towards lowerB until f'(u) > 0
  while(uSlope <= 0) {
    u <- (2 * u + max(-abs(3*u), lowerB)) / 3
    uSlope <- grad(f, x=u, method='simple')
  }
  #If f'(v) >= 0, slide v towards upperB until f'(v) < 0
  while(vSlope >= 0) {
    v <- (2 * v + min(abs(3*v), upperB)) / 3
    vSlope <- grad(f, x=v, method='simple')
  }
  uVal <- log(f(u))
  vVal <- log(f(v))
  
  #Find the intersection point of the tangent lines
  w <- (log(f(v)) - log(f(u)) - v*vSlope + u*uSlope)/(uSlope - vSlope)
  
  #Create g
  g <- rbind(c(lowerB, w, u, uSlope, uVal), c(w, upperB, v, vSlope, vVal))
  colnames(g) <- c('start', 'end', 'intersect', 'm', 'b')
  gLower <- matrix(updateGLower(0, g, f), nrow = 1, ncol = 5)
  colnames(gLower) = c('start', 'end', 'intersect', 'm', 'b')
  return(list(Upper=g, Lower=gLower))
}