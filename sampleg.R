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
