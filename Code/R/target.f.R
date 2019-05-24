target.f <-
function(Ym, Xm, Zm, alpha, R, D, gamma0, gamma1, lambda1, lambda2, lambda3, lambda4, C.m, MLE=TRUE)
{  
  #### Ym: nxKxL
  #### Xm: nxKxP
  #### Zm: nxKxQ
  #### alpha: PxL
  #### b: nxQxL
  #### R: nxKxKxL
  #### D: QxQxL
  #### K is the maximum of the channel# through all the test, some NA for tests do not have K channels         
  #### C.m: PxL
  
  K = dim(Ym)[2];
  n = dim(Ym)[1];
  P = dim(Xm)[3];
  L = dim(Ym)[3];
  Q = dim(Zm)[3];
  
  ###    ch.ind = !is.na(Ym[,,1]);
  ###   D.inv = apply(D,3,solve);
  
  I1T=I2T=0;     
  for(i in 1:n)
  {
    #channel.all = (apply((Ym[i,,]),1,function(x){sum(is.na(x))})!=L);
    #channel.all = (apply((Xm[i,,]),1,function(x){sum(!is.na(x))})==P);
    
    ####index the not all-missing channels
    
    for(l in 1:L)
    {
      m_li=(is.na(Ym[i,,l]));
      if (sum(m_li)<K)
      {
#         Yli = matrix(Ym[i,!m_li,l],nrow=sum(!m_li));
#         Xi = matrix(Xm[i,!m_li,],nrow=sum(!m_li));
#         Zi = matrix(Zm[i,!m_li,],nrow=sum(!m_li));
#         Rli = matrix(R[i,!m_li,!m_li,l],nrow=sum(!m_li));
        Yli = matrix(Ym[i,,l],nrow=K);
        Xi = matrix(Xm[i,,],nrow=K);
        Zi = matrix(Zm[i,,],nrow=K);
        Rli = matrix(R[i,,,l],nrow=K);
        
        al = matrix(alpha[,l],nrow=P);
        Dli = matrix(D[,,l],nrow=Q);
        
        res.temp = (Yli - Xi%*%al)[!m_li];
        cov.temp = (Zi%*%Dli%*%t(Zi)+Rli)[!m_li,!m_li];
        cov.temp = matrix(cov.temp,sum(!m_li))
        I1.temp = -0.5*log(det(cov.temp)) - 0.5*t(res.temp)%*%my.solve(cov.temp)%*%(res.temp);
        I1T = I1T+I1.temp; 
        
      }else{
#         Xi = matrix(Xm[i,channel.all,],nrow=sum(channel.all));
#         Zi = matrix(Zm[i,channel.all,],nrow=sum(channel.all));
#         Rli = matrix(R[i,channel.all,channel.all,l],nrow=sum(channel.all));
        Xi = matrix(Xm[i,,],nrow=K);
        Zi = matrix(Zm[i,,],nrow=K);
        Rli = matrix(R[i,,,l],nrow=K);
        
        al = matrix(alpha[,l],nrow=P);
        Dli = matrix(D[,,l],nrow=Q);
        
        mu_i = Xi%*%al;
        
        cov.temp = Zi%*%Dli%*%t(Zi)+Rli;
        I2.temp = -gamma0-gamma1*mean(mu_i[!is.na(mu_i)])+gamma1^2*mean(cov.temp[!is.na(mu_i),!is.na(mu_i)])/2 #### Pei 7/24
        ########  I2.temp????
        I2T=I2T+I2.temp
        
        #Yi= Xi%*%alpha- gamma1/p*rowSums(Sigma.i)
        #I2 =  gamma0-gamma1*mean(Yi)               
      }

    }
  }
  pen = lambda1*sum(abs((C.m==1)*alpha))+
    lambda2*sum(apply((C.m==1)*alpha,1,function(x){sqrt(sum(x^2))}))+
    lambda3*sum(apply(D,3,function(x){sum(diag(x))}))+
    lambda4*sum(log(apply(D,3,det)));         
  
  return(list(result = I1T+I2T-pen,ll = I1T+I2T));
}
