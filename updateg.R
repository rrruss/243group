
library(numDeriv)
updateGUpper <- function(x, g, f) {
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
  fval <- log(f(x))
  logf <- function(x) { log(f(x)) }
  fprime <- grad(logf, x=x)
  
  # check if x is to the left or right of the intersection toUpdate
  # update either the g to the right or left.
  if (x < g[toUpdate, 'intersect']) {
    left <- toUpdate - 1
    newRangeLeft <- (fval - g[left, 'b'])/(g[left, 'm'] - fprime)
    newRangeRight <- (g[toUpdate, 'b'] - fval)/(fprime - g[toUpdate, 'm'])
    g[left, 'end'] <- newRangeLeft
    g[toUpdate, 'start'] <- newRangeRight
    g <- rbind(g, c(newRangeLeft, newRangeRight, x, fprime, fval))
  }
  if (x > g[toUpdate, 'intersect']) {
    right <- toUpdate + 1
    newRangeRight <- (fval - g[right, 'b'])/(g[right, 'm'] - fprime)
    newRangeLeft <- (g[toUpdate, 'b'] - fval)/(fprime - g[toUpdate, 'm'])
    g[right, 'start'] <- newRangeLeft
    g[toUpdate, 'end'] <- newRangeRight
    g <- rbind(g, c(newRangeLeft, newRangeRight, x, fprime, fval))
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

updateG <- function(x, glist, f){
  # Return a list with elements Upper (the upper bound of g in matrix form)
  #   and Lower (the lower bound of g in matrix form)
  
  # TODO: add ability to accept or reject x.
  # Will add after meeting Thrusday
  # Something like:
  # fx <- f(x)
  # if u < eval(g)/fx { accept x }
  # and maybe pass fx to updateGUpper
  gu <- updateGUpper(x, glist$Upper, f)
  gLower <- updateGLower(x, gu, f)

  return(list(Upper=gu, Lower=gLower))
}

# for testing:
f <- function(x) {
  dnorm(x)
}

source('initg.R')
glist <- initG(f)
glist
glist <- updateG(1.5, glist, f)
glist
glist <- updateG(1, glist, f)
glist