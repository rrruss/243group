checkConcav <- function(fVal, g_upperVal){
  if(fVal > g_upperVal) {
    stop("f is not log-concave; this method will not work.")
  }
}