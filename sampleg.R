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
  
  #correcting for values of hpx approximately equal to 0
  hpzero <- which(sapply(hpx,function(x){all.equal(x,0)})==TRUE)
  areas[hpzero] <- exp(hk[hpzero])*(zk[hpzero+1]-zk[hpzero])
  
  #normalizing factor of the exponentiated upper hull:
  normFactor <- sum(areas)
  scum <- c(0,cumsum(areas)/normFactor)
  
  #generate uniforms to plug into inverse cdf
  u <- runif(n)  
  whichChunk <- sapply(u, function(x) {which(x < scum)[1]-1})
  
  #now for the inverse cdf, broken up into a few pieces for readability
  piece1 <- (hpx[whichChunk])*normFactor*(u-scum[whichChunk])
  piece2 <- log(exp(huz[whichChunk]) + piece1) - huz[whichChunk]
  sample <- zk[whichChunk] + piece2/hpx[whichChunk]
  
  #correcting the NaN's for samples with the uniform in the first chunk
  firstChunk <- which(whichChunk==1)
  piece11 <- (hpx[1])*normFactor*(scum[2]-u[firstChunk])
  piece12 <- log(exp(huz[2]) - piece11) - huz[2]
  sample[firstChunk] <- zk[2] + piece12/hpx[1]
  
  #correcting for values of hpx equal to 0
  for (b in hpzero) {
    chunkhpzero <- which(whichChunk==b)
    sample[chunkhpzero] <- zk[b] + (u[chunkhpzero]-scum[b])*normFactor*(zk[b+1]-zk[b])
  }
  
  return(sample)
}
