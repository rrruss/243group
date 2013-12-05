# For testing:
library(numDeriv)

g <- rbind(c(0,2,1,1,1), c(2,4,3,1,1),
          c(4,8,6,1,1), c(-2, 0, -1.5, 1, 1))
colnames(g) = c('start', 'end', 'intersect', 'm', 'b')
g
g = g[sort.list(g[ ,1], ), ]
x = 1.8
g[, 'start']

f <- function(x) {
  dnorm(x)
}


source('initg.R')
g <- initG(f)

updateg <- function(x, g, f){
  # find index of the function whose range includes x:
  toUpdate = which(g[ ,'start'] <= x & g[ ,'end'] > x)
  fval = log(f(x))
  fnumDer = grad(f, x=x)
  fprime = attr(fnumDer, 'gradient')
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
updateg(1.5, g, f)
grad(exp, x=1.5)
exp(1.5)