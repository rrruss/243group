# Stat 243 Fall 2013 Group Project
#Jeff Darling, Russell Chen, Xuecong Jia, Frank Cleary

#This is the file for the main function of the project

#Placeholder values - will need user to input
n <- 20
f <- "dnorm"
upperB <- 1000
lowerB <- -1000


main <- function(f, n) {
  #Initialize g based on the given f and bounds
  g <- initG(f, upperB, lowerB)
  acceptedX <- NA
  while(length(acceptedX < n)) {
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