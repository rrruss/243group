
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
  
  # checkConcav fails to reject for f = dnorm(1/x), should maybe
  # check that start and end values make sense
  index_Upper <- which(glist$Upper[ ,'start'] <= x & glist$Upper[ ,'end'] > x)
  upperX <- glist$Upper[index_Upper, 'm'] * (x - glist$Upper[index_Upper, 'intersect']) + glist$Upper[index_Upper, 'b']
  fval <- f(x)
  checkConcav(fval, upperX)
  gu <- updateGUpper(x, glist$Upper, f, fval)
  gLower <- updateGLower(x, gu, f)
  return(list(Upper=gu, Lower=gLower, fx=fval))
}

# for testing:
f <- function(x) {
  dnorm(1/x)
}

source('initg.R')
glist <- initG(f)
glist
glist <- updateG(1.5, glist, f)
glist
glist <- updateG(-.1, glist, f)
glist

x2 <- seq(-6, 6, by=.01)
plot(x2, log(f(x2)), type='l', ylim=c(-20,0))
g <- glist$Upper
for (i in 1:nrow(g)){
  x1 <- x2[x2 < g[i, 'end'] & x2 > g[i, 'start']]
  lines(x1, g[i , 'm']*(x1 - g[i , 'intersect']) + g[i , 'b'], col='red')
  abline(v=g[i, 'end'], col='blue')
}
glow <- glist$Lower
for (i in 1:nrow(glow)) {
  x1 <- x2[x2 < glow[i, 'end'] & x2 > glow[i, 'start']]
  lines(x1, glow[i, 'm']*(x1 - glow[i, 'start']) + glow[i, 'b'], col='green')
}
