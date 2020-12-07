l2df <- function(L, len){
  
  x <- c();
  
  for (i in L){
    x <- c(x, rep(i, length.out=len))
  }
  mat <- matrix(x, nrow=len)
  df <- data.frame(mat)
  names(df)<- names(L)
  
  return(df)
}
