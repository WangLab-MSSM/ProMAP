get.alpha <-
function(Y, X, Z, R, b, lambda1, lambda2, Cm, alpha0 = NULL)
{
  #### Y: nxKxL
  #### X: nxKxP
  #### Z: nxKxQ
  #### R: nxKxKxL
  #### b: nxQxL
  #### SigmaHat: KxK
  #### K is the maximum of the channel# through all the test, some NA for tests do not have K channels
  #### C: PxL 
  ####     C_m[i,j]=0 means the corresponding coefficient beta[i,j] is set to be zero in the model; C_m[i,j]=1 means the corresponding beta[i,j] is included in the MAP penalty; C_m[i,j]=2 means the corresponding beta[i,j] is not included in the MAP penalty; default(=NULL): C_m[i,j] are all set to be 1.
  #   Y = Ym.hat;
  #   X = Xm;
  #   Z = Zm;
  #   b = Ebi
  #   Cm = C.m;
  #   alpha0 = alpha.hat;
  
  
  nt = dim(Y)[1];
  L = dim(Y)[3];
  P = dim(X)[3];
  Q = dim(Z)[3];
  K = dim(Y)[2];
  ch.ind = apply(!is.na(Y),c(1,2),sum)!=0;
  
  Y.new = matrix(NA,nrow=nt*K,ncol=L);
  X.new = array(NA,dim=c(nt*K,P,L));
  
  for(i in 1:nt)
  {
    Yi = matrix(Y[i,ch.ind[i,],],sum(ch.ind[i,]),L);  
    Xi = matrix(X[i,ch.ind[i,],],sum(ch.ind[i,]),P);
    Zi = matrix(Z[i,ch.ind[i,],],sum(ch.ind[i,]),Q);
    bi = matrix(b[i,,],Q,L);
    Ri = array(R[i,ch.ind[i,],ch.ind[i,],],dim=c(sum(ch.ind[i,]),sum(ch.ind[i,]),L));
    Y.new[(i-1)*K+(1:K)[ch.ind[i,]],] = (apply(matrix(1:L),1,function(l){sqrt.inv(Ri[,,l])%*%(Yi -(Zi)%*%bi)[,l]}));              
    X.new[(i-1)*K+(1:K)[ch.ind[i,]],,] = vapply(matrix(1:L),function(l){sqrt.inv(Ri[,,l])%*%Xi},outer(1:sum(ch.ind[i,]),1:P));
  }
  
  ind.new = !is.na(Y.new[,1]);
  
  #### assuming all Ri keep the same across different l, so X.new[,,l]=X.new[,,1]
  
  X.m = as.matrix(X.new[ind.new,,1],ncol=P);
  Y.m = as.matrix(Y.new[ind.new,],ncol=L);
  
  lamL1 = lambda1;
  lamL2 = lambda2;
  C.m = Cm;
  rm.est = remMap(X.m, Y.m,lamL1, lamL2, phi0=alpha0, C.m);
  
  ooe = rm.est[[1]];### PxL matrix
  
  var.alpha= NA;### algorithm insert later
  
  return(list(est=ooe, var=var.alpha));
}
