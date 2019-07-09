nmarEM.MLE.fun <-
function(Ym, Xm,Zm, gamma0=0, gamma1, lambda1, lambda2, lambda3, lambda4, C.m, maxIter = 100, MLE=TRUE,tol = 0.0001)  
{
  ### Ym: nxKxL  observed data from n iTRAQ runs;
  ### Xm: nxKxP  for each iTRAQ run Xm[i,,]:
  ###  first row is 1 (indicates global mean);
  ###           2-(K-1)th row are subtype indicator of each sample in this iTRAQ run (reference sample is viewed as one unique subtype);
  ###           last row is for bimodel indicator (model bimodel as a fixed effect)   
  ### Zm: nxKxQfor iTRAQ application, Q=1
  ### gamma1: parameter for missing data mechanism
  ### alpha: PxL
  ### b: nxQxL
  ### R: nxKxKxL
  ### D: QxQxL
  ### K is the maximum of the channel# through all the test, some NA for tests do not have K channels       
  
  print(proc.time()); 
  #Ym.all = NULL;
  #alpha.all = NULL;
  #R.all = NULL;
  #S0.all = NULL; 
  #S1.all = NULL; 
  #D.all = NULL;
#       Ym= Y.obs;
#       Xm=X.temp;
#       Zm=Z.temp;
#       gamma0=gmhat0;
#       gamma1=gmhat1;
#       lambda1=l1;
#       lambda2=l2;
#       lambda3=l3;
#       lambda4=l4;
#       C.m=C.temp;
#       maxIter=It;
#       MLE=TRUE
#       tol=0.0001
  ################################################################################
  ep=0.1;
  ### threshold of initial variance estimator
  
  n=dim(Ym)[1];
  K=dim(Ym)[2];
  L=dim(Ym)[3];
  P=dim(Xm)[3];
  Q=dim(Zm)[3];
  
  Zi = matrix(1,K,Q);
  
  ############ initial values
  
  SS = array(0,dim=c(n,K,K,L));  
  D = array(0,dim=c(Q,Q,L));  
  R = array(0,dim=c(n,K,K,L));  
  b0 = array(0,dim=c(n,Q,L));
  alpha.hat = matrix(0,P,L)
  
  #   for(l in 1:L)
  #   {
  #     Y.temp = Ym[,,l];
  #     SIG.l = cov(Y.temp, use = "pairwise.complete.obs");
  #     ### SIG.l sample covariance matrix
  #     
  #     delta2.l = max(mean(SIG.l[upper.tri(SIG.l)]), ep);
  #     ### delta2.l average covariance     
  #     
  #     D.l = diag(delta2.l, Q);
  #     D[,,l] = D.l;
  #     
  #     R.l = diag(  sapply(diag(SIG.l)-delta2.l,function(x){max(x, ep)}) );
  #     
  #     SS.l = R.l + Zi%*%D.l%*%t(Zi);
  #     for(i in 1:n)
  #     {
  #       SS[i,,,l] = SS.l;
  #       R[i,,,l] = R.l;
  #     }
  #     
  #     Y.temp[is.na(apply(Y.temp,1,sum)),] = NA;    #### throw out the partial missingness
  #     idnObs=(rowSums(is.na(Y.temp))!=K);    #### all obs in K channel
  #     alpha.l = get.alpha0(Y.temp[idnObs,], Xm[idnObs,,], array(Zm[idnObs,,],dim=c(sum(idnObs),dim(Zm)[2],dim(Zm)[3])), 
  #                           SS[1,,,l], matrix(b0[idnObs,,l],nrow=sum(idnObs)) );   
  #     alpha.hat[,l] = alpha.l$est; 
  #     
  #   }
  
  
  Y.new0 = matrix(Ym,ncol=L);
  Y.new = Y.new0[apply(!is.na(Y.new0),1,sum)!=0,]
  X.new0 = matrix(Xm,ncol=P);
  X.new = X.new0[apply(!is.na(Y.new0),1,sum)!=0,]
  #   Y.impu = Y.new;
  #   Y.impu = apply(Y.new,c(2),function(x){xx = x;xx[is.na(x)]=min(x,na.rm=TRUE);return(xx)})
  
  data.raw = list(Y.new);
  names(data.raw) = 'x';
  data.raw$y = rep(1:n,K)[apply(!is.na(Y.new0),1,sum)!=0];
#   Y.new2 = pamr.knnimpute(data.raw ,k = 10, rowmax = 0.8, colmax = 0.8);
  Y.new2 = impute.knn(data ,k = 10, rowmax = 0.9, colmax = 0.9, maxp = 1500);

  Y.impu = Y.new2$data;
  
  for(l in 1:L)
  {
    Y.temp = rep(NA,n*K);  
    Y.temp[apply(!is.na(Y.new0),1,sum)!=0] = Y.impu[,l];
    Y.temp = matrix(Y.temp,n,K);
    SIG.l = cov(Y.temp, use = "pairwise.complete.obs");
    ### SIG.l sample covariance matrix
    
    delta2.l = max(mean(SIG.l[upper.tri(SIG.l)]), ep);
    ### delta2.l average covariance     
    
    D.l = diag(delta2.l, Q);
    D[,,l] = D.l;
    
  }
  
  proc.time();
  
  c0 = C.m;
  c0[C.m==1]=0;  
  
  rm.temp = remMap(X.new, Y.impu, lamL1 = lambda1 ,lamL2 =lambda2 , C.m = c0);
  
#   rm.temp0 = remMap(X.new,Y.impu,lamL1 = lambda1 ,lamL2 =lambda2 , C.m = c0);
#   rm.temp = remMap(X.new,Y.impu,lamL1 = lambda1 ,lamL2 =lambda2 , phi0=rm.temp0[[1]], C.m = C.temp);

  proc.time();
  # rm.temp = remMap(X.new,Y.impu,lamL1 = lambda1 ,lamL2 =lambda2 , C.m = C.temp);
  
  alpha.ini = rm.temp[[1]];
  alpha.hat = alpha.ini;
  res.ini = matrix(NA,n*K,L);
  res.ini[apply(!is.na(Y.new0),1,sum)!=0,] = Y.impu-X.new%*%alpha.ini;
  res.ini = array(res.ini,dim=c(n,K,L));
  # var.ini = apply(res.ini,c(2),function(x){var(c(x))});
  var.ini = sapply(apply(res.ini,c(2),function(x){var(c(x),na.rm=TRUE)})-mean(D),function(x){max(x,ep)});  
  
  RK = Diagonal(n = 4,x = var.ini);
  
  R = array(rep(rep(matrix(RK),each = n),L),dim=c(n,K,K,L));    
  
  Ym.hat=Ym;
  ############# constant   
  ###maxIter=100;
  
  #tol=0.0001;
  diff=999;
  fe.diff = 999; 
  iter=0;
  nllik=NULL;
  oll = NULL;
  ###dellis=delta2.hat
  ###siglis=c(sigma0.hat, sigma2.hat)
  ###alphalis=alpha.hat
  
  lliklis=NULL;
  
  
  Ebi = array(NA, dim=c(n,Q,L)); 
  Varbi = array(NA, dim=c(n,Q,Q,L));
  Eei = array(NA, dim=c(n,K,L));
  Varei = array(NA, dim=c(n,K,K,L));
  
  print(proc.time());
  
  ############# START ECM
  while (iter<maxIter & as.logical((diff>tol)+(fe.diff>0))) #for(iter in 0:maxIter),relative diff,feature identified;
  {
    iter = iter+1;        
    
    ######################
    #### E-step: 
    #time7 = proc.time();    
    for (i in 1:n)
    {    
      #kmax = max(apply(!is.na(Ym[i,,]),2,sum));
      ###lmax = which.max(apply(matrix(!is.na(Xm[i,,]),K,P),1,sum));
      
      for(l in 1:L)
      {
        m.temp=sum(!is.na(Ym[i,,l]));
        if (m.temp>0)
        {
          #### E(bi|yobs, theta.hat)
          #### E(ei|yobs, theta.hat)
          
          idxo = (!is.na(Ym[i,,l]));
          ki=sum(idxo);
          Yli = matrix(Ym[i,idxo,l],nrow=ki);
          Xli = matrix(Xm[i,idxo,],nrow=ki);
          Zli = matrix(Zm[i,idxo,],nrow=ki);
          Dli = matrix(D[,,l],nrow=Q);
          Sigma.li = Zli%*%Dli%*%t(Zli) + (R[i,idxo,idxo,l]);
          Wli  = my.solve(Sigma.li);
          
          Ebi[i,,l]=Dli%*%t(Zli)%*%Wli%*%(Yli-Xli%*%alpha.hat[,l]); 
          
          ###Eei[i,idxo,l]=as.vector(Yli- Xli%*%alpha.hat[,l] - Zli%*% Ebi[i,,l] );
          
          Varbi[i,,,l] = Dli - Dli%*%t(Zli)%*%Wli%*%Zli%*%Dli; 
          Varei[i,idxo,idxo,l] =  Zli%*%matrix(Varbi[i,,,l],nrow=Q)%*%t(Zli);
          
        } else {
          
          #### E(Yi|Mi=1, theta.hat)
          #### E(bi|Mi=1, theta.hat)
          #### E(ei|Mi=1, theta.hat)                    
          #### idxo = !is.na(Ym[i,,lmax]);
          idxo = apply(!is.na(Ym[i,,]),1,sum)!=0
          #idxo = (apply(matrix(is.na(Xm[i,,]),K,P),1,sum))==0;              
          Xli= matrix(Xm[i,idxo,],nrow=sum(idxo));
          Zli= matrix(Zm[i,idxo,],nrow=sum(idxo));
          Dli = matrix(D[,,l],nrow=Q);    
          Sigma.li = Zli%*%Dli%*%t(Zli) + R[i,idxo,idxo,l];
          Ebi[i,,l] = - gamma1*rowMeans(Dli%*%t(Zli));
          
          ###Eei[i,,l] = - gamma1* colMeans(R[i,idxo,idxo,l]); 
          ### Eei????
          
          Varbi[i,,,l]= D[,,l];
          Varei[i,idxo,idxo,l]= R[i,idxo,idxo,l];
          Ym.hat[i,idxo,l]= t(Xli%*%alpha.hat[,l] - gamma1*rowMeans(Sigma.li));
        }
      }
    }
    #time8 = proc.time()[[3]];        
    ###################################
    ###### M step     
    
    #### ??? we may need an iteration in the M-step between estimating alpha and R, but even without that iteration, the likelihood can still monotone increase. 
    
    alpha.result=get.alpha(Ym.hat, Xm, Zm, R, Ebi, lambda1, lambda2, C.m, alpha0 = alpha.hat); ### Ym.hat, R, Ebi are from last iteration;
    alpha.hat.new=matrix(alpha.result$est, ncol=L);
#     proc.time();  
    fe.diff = sum((alpha.hat.new!=0)!=(alpha.hat!=0));
    #fn.diff = sqrt(sum((alpha.hat.new[-(1:2),]-alpha.hat[-(1:2),])^2));
    
    alpha.hat = alpha.hat.new;
    
    #time9 = proc.time()[[3]];            
    for (l in 1:L)
    {
      D[,,l]  =  (2*lambda3*diag(Q)+t(Ebi[,,l])%*%(Ebi[,,l]) + apply(array(Varbi[,,,l],dim=c(n,Q,Q)),c(2,3),sum))/(n+2*lambda4);
    }  
    #    idxo = !is.na(Ym.hat[i,,l]);
    #    Yl= matrix(Ym.hat[,idxo,l],nrow=n);
    #    Xi= matrix(Xm[i,idxo,],nrow=sum(idxo));
    #    Zi= matrix(Zm[i,idxo,],nrow=sum(idxo));  
    #    Eei[i,idxo,l]=as.vector( Yli- Xi%*%alpha.hat[,l] - Zi%*% matrix(Ebi[i,,l]) );
    
    for(k in 1:K)
    {       
      idok = !is.na(Ym.hat[,k,]);
      Yk = matrix(Ym.hat[,k,],nrow=n);
      Xk = matrix(Xm[,k,],nrow=n);
      Zk = matrix(Zm[,k,],nrow=n);  
      random.k = vapply(matrix(1:n),function(i){matrix(Zk[i,],nrow=1)%*%matrix(Ebi[i,,],nrow=Q)},outer(1,1:L));
      Eei[,k,][idok]=as.vector( Yk- Xk%*%alpha.hat - t(matrix(random.k,nrow=L)) )[idok];
    }
    
    ###   R.l = cov(Eei[,,l], use = "pairwise.complete.obs");
    
    ###    
    df0= sum(!is.na(Eei[,1,]));
    df1= sum(!is.na(Eei[,-1,]));     
    sigma0.hat = (sum(Varei[,1,1,],na.rm=T) + sum(Eei[,1,]^2,na.rm=T))/df0;  
    sigma1.hat = (sum(diag( apply(Varei[,-1,-1,], c(2,3), sum, na.rm=T)))+ sum(Eei[,-1,]^2,na.rm=T))/df1;
    R0 = diag(c(sigma0.hat, rep(sigma1.hat,K-1)));   
    ###    
    
    ### faster than loop  
    R = array(rep(rep(R0,each = n),L),dim=c(n,K,K,L));    
    
    ###dellis=c(dellis, D)
    ###siglis=cbind(siglis, c(sigma0.hat, sigma2.hat))
    ###alphalis=cbind(alphalis, alpha.hat)
    
    nllik.temp=target.f(Ym=Ym, Xm=Xm, Zm=Zm, alpha=alpha.hat, R=R, D=D,  gamma0=gamma0, gamma1=gamma1, lambda1, lambda2, lambda3, lambda4, C.m, MLE=MLE);
    nllik=c(nllik, nllik.temp$result);
    oll = c(oll,nllik.temp$ll);
    ###lliklis=cbind(lliklis, nllikt);
    if (length(nllik)>=2) diff= abs(nllik[iter-1]-nllik[iter])/max(abs(nllik[iter-1]),1);
    
    #  ###############################################################
    #  Ym.all = c(Ym.all,list(Ym.hat));
    #  alpha.all = c(alpha.all,list(alpha.hat));
    #  R.all = c(R.all,list(R)); 
    #  D.all = c(D.all,list(D));
    
    #     write.table(c(fe.diff,fn.diff),file = paste('U_CV_diff_true',ii,'.txt',sep=''),append=TRUE);
    #     write.table(cbind(nllik,oll),file = paste('U_CV_llh_true',ii,'.txt',sep=''));
    #     write.table(alpha.hat,file = paste('U_CV_alpha_true',ii,'.txt',sep=''));
    #     write.table(c(sigma1.hat, sigma0.hat, D),file = paste('U_CV_var_true',ii,'.txt',sep=''));
    print(c(iter,sum(alpha.hat!=0),fe.diff,diff,max(D),proc.time()[3]));
    #print(c(iter,fe.diff,diff,proc.time()));
    
    if(max(D)>1000)  
      break;
  }
  
  #write.table(cbind(nllik,oll),file = paste('P:/work/aaa/T_CV_llh_true',iter,'.txt',sep=''));
  
  #proc.time();
  alpha.var = NA;
  
  return(list(nllik=nllik, alpha.hat=alpha.hat, alpha.var=alpha.var, sigma1.hat=sigma1.hat, sigma0.hat=sigma0.hat, D=D, oll=oll));
  
}
