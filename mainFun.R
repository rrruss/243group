# Stat 243 Fall 2013 Group Project
#Jeff Darling, Russell Chen, Xuecong Jia, Frank Cleary

#This is the file for the main function of the project

n <- 20
f <- "dnorm"
upperB <- 1000
lowerB <- -1000
g <- initG(f, upperB, lowerB)

main <- function(f, n) {
  acceptedX <- NA
  while(length(acceptedX < n)) {
    if(is.na(acceptedX)) {
      acceptedX <- generatePoints(n - length(acceptedX))
    } else {
      acceptedX <- c(acceptedX, generatePoints(n-length(acceptedX)))
    }
  }
  return(acceptedX)
}

sample <- main(f, n)