arSample <- function(f, n, showg = FALSE, upperB=Inf, lowerB=-Inf, initVals=NULL) {
  library(numDeriv)
  #Initialize g based on the given f and bounds
  glist <- initG(f, upperB, lowerB, initVals)
  acceptedX <- NULL
  while(length(acceptedX) < n) {
    updates <- generatePoints(min(n - length(acceptedX), 10*nrow(glist$Upper)), glist, f)
    #Append acceptedX with newly generated points
    acceptedX <- c(acceptedX, updates$sample)
    glist <- updates$g
  }
  if (showg) {
    plotg(f, glist)
  }
  return(acceptedX)
}