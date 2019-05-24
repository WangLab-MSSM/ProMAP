ProMap <-
function(X.data, Y.data, num.channel, index.experiment, index.reference, lambda.1=lambda.1, lambda.2=lambda.2, nu.1=nu.1, nu.2=nu.2, V=V, CV.fold=CV.fold, max.iteration = max.iteration , C.matrix=NULL, tolerance=tolerance )
{
# lambda.1: lambda_1
# lambda.2: lambda_2
# nu.1: nu_1
# nu.2: nu_2
# CV.fold : CV. fold
# It: the maximum number of iteration

K = num.channel;
n = dim(Y.data)[1]/K;


## check step 1: sample numbers in X and Y and index should be equal
if(  ((dim(X.data)[1]==dim(Y.data)[1])+(dim(Y.data)[1]==length(index.experiment))+(dim(Y.data)[1]==length(index.reference)))!=3 )
   {
    stop("Dimension of data variables disagree.");
   }

## check step 2: index.experiment should be the same set as ( rep(1:num.channel,dim(Y.data)[1]) )
if( !identical(rep(1:n,each = K) , sort(index.experiment)) )
   {
    stop("Experiment index is incorrect.");
   }

## check step 3: index.reference[!is.na(index.reference)] should be the same set as 1:dim(Y.data)[1]??
if( !identical(1:n, sort(index.reference[!is.na(index.reference)] )) )
   {
    stop("Reference index is incorrect.");
   }

## check step 4: index.reference[!is.na(index.reference)] should be the same set as??index.experiment[!is.na(index.reference)]
if( !identical(index.experiment[!is.na(index.reference)] , index.reference[!is.na(index.reference)]) ) 
   {
    stop("Reference index and experiment index disagree.");
   }


if( CV.fold>V) 
{
  stop("Invalid CV fold.");
}


index.temp <- index.experiment-(!is.na(index.reference))/100;
### in this step we make reference a little bit different than the regular channels so that when we sort them we can always let the reference to be the first one in the experiment.

order.sample = c(t(matrix(order(index.temp),num.channel)));

X.array = array(as.matrix(X.data[order.sample,]),dim = c(n,K,dim(X.data)[2]));
Y.array = array(as.matrix(Y.data[order.sample,]),dim = c(n,K,dim(Y.data)[2]));



# Change data format
#X.array = array(as.matrix(X.data),dim = c(n,K,dim(X.data)[2]));
#Y.array = array(as.matrix(Y.data),dim = c(n,K,dim(Y.data)[2]));


########################
# Sep 0: Imputation
########################
# Define CV fold

P = dim(X.array)[3];
L = dim(Y.array)[3];
K = dim(Y.array)[2];
n = dim(Y.array)[1];

seed.temp = 1111;
set.seed(seed.temp);

# n.ind: index of CV fold
n.ind = (sample(1:n))[-((CV.fold-1)*n/V+1:(n/V))];
# nf: sample size of CV fold
nf = length(n.ind);

# X.temp: CV fold data of X
X.temp = X.array[n.ind,,];
# Y.temp: CV fold data of Y
Y.obs = Y.array[n.ind,,];

# C.matrix: It can have only values, 0, 1, and Inf
if( is.null(C.matrix) )
{
C.matrix = matrix(1,P,L);
C.matrix[1,] = C.matrix[2,] = 0;
C.temp <- 2 - C.matrix
C.temp[C.temp==-Inf] <- 0
}
if ( !is.null(C.matrix) ){ 
   if(  sum( !( C.matrix==0 | C.matrix==1 | C.matrix==Inf ) )>0 )
   {
    stop("C.matrix can have only 0, 1, or Inf.")
   }
}


# z.temp: coefficient of random effect
Z.temp = array(1,dim=c(length(n.ind),K,1));


#######
# Check missing values
total_sum<-matrix( 0 , nrow=1, ncol=dim(Y.obs[,1,])[2] )
for(i in 1:4)
{
total_sum<-total_sum+apply(Y.obs[,i,],2,sum)
}

if( sum(is.na(total_sum)) ==0  )
{
gmhat0 = Inf
gmhat1 = 0
}

if( sum(is.na(total_sum)) ==1  )
{
gmhat0 = -log(mean(is.na(Y.obs)))
gmhat1 = 0
}

###########

if( sum(is.na(total_sum)) >=2  )
{
# Average of channel-wise
Y.imp = apply(Y.obs,c(1,3),mean);
Y.imp = Y.imp[apply(!is.na(Y.imp),1,sum)!=0,];

# Imputation
minL = apply((apply(Y.obs,c(1,3),function(x){mean(x,na.rm=TRUE)})),2,function(x){mean(x,na.rm=TRUE)});
for(l in 1:L)
{
  Y.imp[,l][is.na(Y.imp[,l])] = minL[l];    
}

y.mi2 = apply(Y.imp,2,median);
mr.l = apply(as.matrix(is.na(Y.obs[,1,]),n,L),2,mean);


ind.mr = c(mr.l)!=0
gmhat0 = -lm(c(log(mr.l))[ind.mr]~c(y.mi2)[ind.mr])$coef[1];
gmhat1 = -lm(c(log(mr.l))[ind.mr]~c(y.mi2)[ind.mr])$coef[2];

if( length(unique(c(y.mi2)[ind.mr]))==1 )
{
gmhat0 = -log(mean(is.na(Y.obs)))
gmhat1 = 0
}

}



test.result.hat = nmarEM.MLE.fun(Ym= Y.obs, Xm=X.temp, Zm=Z.temp, gamma0=gmhat0, gamma1=gmhat1, lambda1=lambda.1, lambda2=lambda.2, lambda3=nu.1, lambda4=nu.2, C.m=C.temp, maxIter=max.iteration, MLE=TRUE,tol=tolerance);  

class(test.result.hat)<-"ProMap"
test.result.hat

}
