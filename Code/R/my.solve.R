my.solve <-
function(X){
  if (!is.matrix(X))  X <- matrix(X, nrow=sqrt(length(X)))
  Xinv<-ginv(X)
  return(Xinv)

  if(FALSE)
  {
  if (!is.matrix(X))  X <- matrix(X, nrow=sqrt(length(X)))
  ss <- svd(X)
  dd <- ss$d
  dd[abs(dd)<10^(-12)] <- 10^(-12)
  ###dd[abs(dd)<10^(-12)] <- 10^(-8)
  Xinv <- ss$u%*%diag(1/dd, nrow=nrow(X), ncol=nrow(X))%*%t(ss$v)
  return(Xinv)
  }
}
