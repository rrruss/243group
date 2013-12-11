# intersectnP <- function(xk,hk,hpx,xlb,xub) {
#   #this function finds the intersection points of the tangents
#   #presumes the following vectors of length k:
#   #xk is a vector of the k abscissae x_i
#   #hk is a vector of h(x) the k abscissae x_i
#   #hpx is a vector of h'(x) of the k abscissae x_i
#   #and two scalar values:
#   #z_0 = xlb lower bound for domain D
#   #z_k = xub upper bound for domain D
#   k <- length(xk)
#   zk <- rep(NA,k+1)
#   zk[1] <- xlb
#   zk[n] <- xub
#   #TO DO: vectorize the following
#   for (i in 2:(k-1)) {
#     zk[i] = xk[i] + (hk[i]-hk[i+1]+hpx[i+1]*(xk[i+1]-xk[i]))/(hpx[i+1]-hpx[i])
#   }
#   return(zk)
# }

sampleSX <- function(g,n) {
  #this function draws n random samples from the exponentiated upper envelope
  #g is a matrix where
  #xk is a vector of the k abscissae x_i
  #hk is a vector of h(x)=log(f(x)) the k abscissae x_i
  #hpx is a vector of h'(x) of the k abscissae x_i
  #zk is the intersection points of the tangents
  xk <- g[,'intersect']
  k <- length(xk)
  hk <- g[,'b']
  hpx <- g[,'m']
  zk <- c(g[,'start'],g[k,'end']) #note that there are (k+1) z values
  
  #calculate the areas of the k chunks (after exponentiating):
#   areas <- rep(NA,k)
#   #TO DO: vectorize
#   for (i in 1:k) {
#     areas[i] <- (exp(zk[i+1])-exp(zk[i]))/hpx[i]
#   }
  
  #to calculate the value of the upper envelope at the z points, interpolate using the tangent lines:
  xdummy <- c(0,xk)
  zkMinusxk <- zk - xdummy
  huz <- hpx*(zkMinusxk[-1]) +hk
  #the first value (at zk[1]) is calculated manually:
  huz0 <- huz[1] - hpx[1]*(zk[2]-zk[1])
  huz <- c(huz0,huz)
  
  #here it is, a vectorized version of computing the areas:
  t <- exp(c(huz,1)) - exp(c(1,huz))
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
  for (j in 1:n) {
    while (u[j] > scum[whichChunk[j]+1])
    whichChunk[j] <- whichChunk[j] + 1
  }
  #now for the inverse cdf, broken up into a few pieces for readability
  piece1 <- (hpx[whichChunk])*normFactor*(u-scum[whichChunk])
  piece2 <- log(exp(huz[whichChunk]) + piece1) - huz[whichChunk]
  sample <- zk[whichChunk] + piece2/hpx[whichChunk]
  
  #correcting the NaN's for samples with the uniform in the first chunk
  #this is a short for loop for hulls with many components
  for (m in which(whichChunk==1)) {
    piece11 <- (hpx[1])*normFactor*(scum[2]-u[m])
    piece12 <- log(exp(huz[2]) - piece11) - huz[2]
    sample[m] <- zk[2] + piece12/hpx[1]
  }
  return(sample)
}
