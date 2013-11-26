# For testing:

g <- rbind(c(0,2,1,1,1), c(2,4,3,1,1),
          c(4,8,6,1,1), c(-2, 0, -1.5, 1, 1))
colnames(g) = c('start', 'end', 'intersect', 'm', 'b')
g
g = g[sort.list(g[ ,1], ), ]
x = .8
g[, 'start']

f <- function(x) {
  exp(5 * x)
}

updateg <- function(x, g, f){
  # find index of the function whose range includes x:
  toUpdate = which(g[ ,'start'] <= x & g[ ,'end'] > x)
  fval = log(f(x))
  fnumDer = numericDeriv(quote(log(f(x))), 'x')
  fprime = attr(fnumDer, 'gradient')
  # check if x is to the left or right of the intersection toUpdate
  if (x < g[toUpdate, 'intersect']) {
    # i.e. x is to the left:
    leftg <- toUpdate - 1
    # calculate intersections:
    newRangeLeft <- (fval - g[leftg, 'b'])/(g[leftg, 'm'] - fprime)
    newRangeRight <- (g[toUpdate, 'b'] - fval)/(fprime - g[toUpdate, 'm'])
    g[leftg, 'end'] <- newRangeLeft
    g[toUpdate, 'start'] <- newRangeRight
    g <- rbind(g, c(newRangeLeft, newRangeRight, x, fprime, fval))
  }
  if (x > g[toUpdate, 'intersect']) {
    # i.e. x is to the right:
    rightg <- toUpdate + 1
    # calculate intersections:
    newRangeRight <- (fval - g[rightg, 'b'])/(g[rightg, 'm'] - fprime)
    newRangeLeft <- (g[toUpdate, 'b'] - fval)/(fprime - g[toUpdate, 'm'])
    g[rightg, 'end'] <- newRangeRight
    g[toUpdate, 'start'] <- newRangeLeft
    g <- rbind(g, c(newRangeLeft, newRangeRight, x, fprime, fval))
  }
  g <- g[sort.list(g[ ,1], ), ]
  return(g)
}
g