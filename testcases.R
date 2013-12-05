# Test normal distribution:
testdist <- function(f, testname, trueMean, trueVar, n=1e5){
  y <- mainFun(f, n)
  print mean(y)
  print var(y)
  if (mean(y) - trueMean > .1) {
    print(paste("Failed test", testname, 'with mean', mean(y)))
  }
  if (var(y) - trueVar > .1) {
    print(paste("Failed test", testname, 'with variance', var(y)))
  }
}
testdist(dnorm, 'normal', 0, 1)