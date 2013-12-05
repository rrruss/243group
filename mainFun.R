# Stat 243 Fall 2013 Group Project
#Jeff Darling, Russell Chen, Xuecong Jia, Frank Cleary

#This is the file for the main function of the project

#Placeholder values for testing - will need user to input
n <- 20
f <- function(x) {
  dnorm(x)
}
upperB <- 1000
lowerB <- -1000


main <- function(f, n) {
  #Initialize g based on the given f and bounds
  g <- initG(f, upperB, lowerB)
  acceptedX <- NA
  while(length(acceptedX) < n) {
    if(is.na(acceptedX)) {
      #For the first run, replace NA with points generated
      acceptedX <- generatePoints(n - length(acceptedX), upperG, lowerG)
    } else {
      #For subsequent runs, lengthen existing vector
      acceptedX <- c(acceptedX, generatePoints(n-length(acceptedX)))
    }
  }
  return(acceptedX)
}

sample <- main(f, n)

######################initg.R######################
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
  uVal <- log(f(u))
  vVal <- log(f(v))
  
  #Find the intersection point of the tangent lines
  w <- (log(f(v)) - log(f(u)) - v*vSlope + u*uSlope)/(uSlope - vSlope)
  
  #Create g
  g <- rbind(c(lowerB, w, u, uSlope, uVal), c(w, upperB, v, vSlope, vVal))
  colnames(g) = c('start', 'end', 'intersect', 'm', 'b')
  return(list(Upper=g))
}
############################################

######################checkconcav.R######################
checkConcav <- function(fVal, g_upperVal){
  if(fVal > g_upperVal) {
    stop("f is not log-concave; this method will not work.")
  }
}
############################################

######################updateg.R######################

library(numDeriv)
updateGUpper <- function(x, g, f, fval) {
  ### Return g with an added row with intersect = x
  ### Calculates which other rows of g need updating and updates them
  ### 
  ### x: The point to update at
  ### f: The distribution being drawn from
  ### g: A matrix where each row's values are start, end, intersect, m and b
  ###   where start and end define the range g applies to, intersect is the
  ###   point x at which g is tangent to log(f(x)), and b is the value log(f(x))
  
  # find index of the function whose range includes x:
  toUpdate <- which(g[ ,'start'] <= x & g[ ,'end'] > x)
  logfval <- log(fval)
  logf <- function(x) { log(f(x)) }
  fprime <- grad(logf, x=x)
  
  # check if x is to the left or right of the intersection toUpdate
  # update either the g to the right or left.
  if (x < g[toUpdate, 'intersect']) {
    left <- toUpdate - 1
    gl <- g[left, ]
    newRangeLeft <- (logfval - gl['b'] + 
                       gl['m']*gl['intersect'] - fprime*x)/(gl['m'] - fprime)
    newRangeRight <- (g[toUpdate, 'b'] - logfval +
                        fprime*x - g[toUpdate, 'm']*g[toUpdate, 'intersect'])/(fprime - g[toUpdate, 'm'])
    g[left, 'end'] <- newRangeLeft
    g[toUpdate, 'start'] <- newRangeRight
    g <- rbind(g, c(newRangeLeft, newRangeRight, x, fprime, logfval))
  }
  if (x > g[toUpdate, 'intersect']) {
    right <- toUpdate + 1
    newRangeRight <- (logfval - g[right, 'b'] + 
                        g[right, 'm']*g[right, 'intersect'] - fprime*x)/(g[right, 'm'] - fprime)
    newRangeLeft <- (g[toUpdate, 'b'] - logfval + 
                       fprime*x - g[toUpdate, 'm']*g[toUpdate, 'intersect'])/(fprime - g[toUpdate, 'm'])
    g[right, 'start'] <- newRangeLeft
    g[toUpdate, 'end'] <- newRangeRight
    g <- rbind(g, c(newRangeLeft, newRangeRight, x, fprime, logfval))
  }
  g <- g[sort.list(g[ ,1], ), ]
  return(g)
}

updateGLower <- function(x, gu, f) {
  ### Return the lower bound glower, based on gu, the upper bound
  ### Returns a matrix where each row is contain start, end, intersect, m, and b
  ### start and end define the range that g lower is valid on
  ### intersect is not used
  ### m is the slope of the line
  ### b is the value of log(f(x)) at x = 'start', which can be used to
  ###   calculate the equation of the line.
  
  glower <- gu
  glower[ , 'start'] <- gu[ , 'intersect']
  glower[ , 'end'] <- c(gu[-1, 'intersect'], Inf)
  # calculate slope:
  glower[ , 'm'] <- c(diff(glower[ , 'b']), 0)/(glower[ , 'end'] - glower[ , 'start'])
  glower <- glower[-nrow(glower), ] # remove last row - it's meaningless
  glower # note the b column is the value of log(f) at 'start'
}

checkConcav <- function(fval, g_upperVal){
  if(fval > g_upperVal) {
    stop("f is not log-concave; this method will not work.")
  }
}

updateG <- function(x, glist, f){
  # Return a list with elements Upper (the upper bound of g in matrix form)
  #   and Lower (the lower bound of g in matrix form)
  
  # TODO: add ability to accept or reject x.
  # Will add after meeting Thrusday
  # Something like:
  # fx <- f(x)
  # if u < eval(g)/fx { accept x }
  # and maybe pass fx to updateGUpper
  index_Upper <- which(glist$Upper[ ,'start'] <= x & glist$Upper[ ,'end'] > x)
  upperX <- glist$Upper[index_Upper, 'm'] * (x - glist$Upper[index_Upper, 'intersect']) + glist$Upper[index_Upper, 'b']
  fval <- f(x)
  checkConcav(fval, upperX)
  gu <- updateGUpper(x, glist$Upper, f, fval)
  gLower <- updateGLower(x, gu, f)
  return(list(Upper=gu, Lower=gLower, fx=fval))
}
############################################

######################generatePoints.R######################
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
      if (U[i] < glist$fx / exp(upperX)){
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
############################################

######################sampleg.R######################
sampleSX <- function(g,n) {
  #this function draws n random samples from the exponentiated upper envelope
  #g is a matrix where
  #xk is a vector of the k abscissae x_i
  #hk is a vector of h(x) the k abscissae x_i
  #hpx is a vector of h'(x) of the k abscissae x_i
  #zk is the intersection points of the tangents
  xk <- g$x
  k <- length(xk)
  hk <- g$
    hpk <- (g$slope)*(g$x) + g$intercept
  zk <- c(g$start,g$end[k]) #note that there are (k+1) z values
  
  #calculate the areas of the k chunks (after exponentiating):
  #   areas <- rep(NA,k)
  #   #TO DO: vectorize11111111111111
  #   for (i in 1:k) {
  #     areas[i] <- (exp(zk[i+1])-exp(zk[i]))/hpx[i]
  #   }
  #here it is, a vectorized version of computing the areas:
  t <- exp(c(zk,1)) - exp(c(1,zk))
  areas <- t[-c(1,k+2)]/hpx
  
  #normalizing factor of the exponentiated upper hull:
  normFactor <- sum(areas)
  scum <- c(0,cumsum(areas)/normFactor)
  
  #   u <- runif(1)
  #   j <- 1
  #   while (u > scum[j+1]) {
  #     j <- j+1
  #   }
  #   sample <- zk[j] + log(1 + (hpx[j+1])*normFactor*(u-scum[j]) / (exp(zk[j])) )/hpx[j+1]
  #   return(sample)
  
  #try a vectorized version:
  u <- runif(n)
  whichChunk <- rep(1,n)
  for (k in 1:n) {
    while (u[k] > scum[whichChunk[k]+1])
      whichChunk[k] <- whichChunk[k] + 1
  }
  sample <- zk[whichChunk] + log(1 + (hpx[whichChunk+1])*normFactor*(u-scum[whichChunk]) / (exp(zk[whichChunk])) )/hpx[whichChunk+1]
  return(sample)
}
############################################

