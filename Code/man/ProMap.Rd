\name{ProMap}
\alias{ProMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{   ProMap: A Regularized Linear Mixed Effect Model (Mixed remMap) 
%%  ~~function to do ... ~~
}
\description{
This package estimates a regularized linear mixed effect model for hub networks.
}
\usage{
ProMap <-
function(X.data, Y.data, num.channel, index.experiment, index.reference , lambda.1=lambda.1, lambda.2=lambda.2, nu.1=nu.1, nu.2=nu.2, V=V, CV.fold=CV.fold, max.teration = max.teration , C.matrix=NULL, tolerance=tolerance )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X.data}{  input data matrix with row n times K and column P. n indicates an experimental sample in each channel.  K indicates the number of channel. P indicates the number of upper biomarker.
%%     ~~Describe \code{Xarray} here~~
}
  \item{Y.data}{ input data matrix with row n times K and column L. n indicates an experimental sample in each channel.  K indicates the number of channel. L indicates the number of lower biomarker.
%%     ~~Describe \code{Yarray} here~~
}
  \item{num.channel}{ the number of channel
%%     ~~Describe \code{Yarray} here~~
}
  \item{index.experiment}{ the index of experiment
%%     ~~Describe \code{Yarray} here~~
}
  \item{index.reference}{ the index of reference
%%     ~~Describe \code{Yarray} here~~
}


  \item{lambda.1}{ the penalty parameter (\code{lambda_1}) for the L1 norm to control the total number of parameters in the coefficient matrix. It controls sparsity.%%     ~~Describe \code{lambda.1} here~~
}
  \item{lambda.2}{ the penalty parameter (\code{lambda_2}) for the L2 norm to group th coefficients in the same row together of the coefficient matrix. (between the same predictor and different responses) Thus, they are more likely to be included or excluded together in the model estimation, which indicates group effect.  %%     ~~Describe \code{lambda.2} here~~
}
  \item{nu.1}{  the first tunning parameter (\code{nu_1}) on the covariance matrix
}
  \item{nu.2}{  the second tunning parameter (\code{nu_2}) on the covariance matrix
}
  \item{V}{ the total partition number in CV
}


  \item{CV.fold}{ the number of cross validation used to select the value of penalty parameters
}
  \item{max.teration}{ the maximum number of iteration used to obtain the solution
}
  \item{tolerance}{ the interval between two consequent objective function values to stop the algorithm
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
   \item{nllik}{ the penalized log likelihood value }
  \item{alpha.hat}{ the estmated coeficient matrix PxL   }
  \item{sigma1.hat}{ the variance estimation of the error term of non-reference channel }
  \item{sigma0.hat}{ the variance estimation of the error term of reference channel }
  \item{D}{ the variance estimation of the random effect, vector of length L }
  \item{oll}{ the log likelihood value without penalties }
  
}
\references{
Ma, W., Chen, L.S., Han, S.W., Lin, C., Yan, P., Paulovich, A., CPTAC consortium, Zhong, H., and Wang, P. (2015) \emph{A Penalized Multivariate Linear Mixed Effects Model for Integrative Proteo-genomics Analysis Based on iTRAQ Data}

}
\author{
Weiping Ma, Sung Won Han, Hua Zhong, and Pei Wang

Maintainer: Pei Wang <pei.wang@mssm.edu>
}

\examples{

library(mvtnorm);
library(corpcor);
library(MASS);
library(remMap);
library(pamr);
library(Matrix);

# call simulation data for X matrix          
data(NNX)

# call simulation data for Y matrix, complete data
data(NNYcomp)

# call simulation data for Y matrix, missing data
data(NNYobs)

# call the index of experiment
data(NNX.index.experiment)

# call the index of reference
data(NNX.index.reference)

### total CV patition number
V = 5;


# Grid the parameter values
L1 = c(11:20);
L2 = seq(0,90,by =10);
Lam = expand.grid(L1,L2);
# Job number
ii = 320 ;

# Obtain solutions with complete data
ProMap.comp <- ProMap(X.data=NNX, Y.data=NNYcomp, num.channel=4, index.experiment=NNX.index.experiment, index.reference=NNX.index.reference , lambda.1=Lam[ceiling(ii/4),1], lambda.2=Lam[ceiling(ii/4),2], nu.1=0.01, nu.2=1, V=V, CV.fold=floor((ii-1)\%\%(V)), max.iteration =100, tolerance=10^-4)


# Obtain solutions with missing data
ProMap.obs <- ProMap(X.data=NNX, Y.data=NNYobs, num.channel=4, index.experiment=NNX.index.experiment, index.reference=NNX.index.reference , lambda.1=Lam[ceiling(ii/4),1], lambda.2=Lam[ceiling(ii/4),2], nu.1=0.01, nu.2=1, V=V, CV.fold=floor((ii-1)\%\%(V)), max.iteration =100, tolerance=10^-4)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
