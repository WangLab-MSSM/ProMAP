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
ProMap(Xarray = Xarray, Yarray = Yarray, lambda.1 = lambda.1, lambda.2 = lambda.2, nu.1 = nu.1, nu.2 = nu.2, CV.i = CV.i, CV.all = CV.all, max.teration = max.teration, seed.CV = 1234)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xarray}{
%%     ~~Describe \code{Xarray} here~~
}
  \item{Yarray}{
%%     ~~Describe \code{Yarray} here~~
}
  \item{lambda.1}{
\code{lambda_1}%%     ~~Describe \code{lambda.1} here~~
}
  \item{lambda.2}{
\code{lambda_2}%%     ~~Describe \code{lambda.2} here~~
}
  \item{nu.1}{
\code{nu_1}
}
  \item{nu.2}{
\code{nu_2}
}
  \item{CV.i}{
current CV fold
}
  \item{CV.all}{
total CV folds
}
  \item{max.teration}{
the maximum number of iteration
}
  \item{seed.CV}{
random seed to generate CV folds
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
   \item{nllik}{  }
  \item{alpha.hat}{    }
  \item{alpha.var}{  }
  \item{sigma1.hat}{  }
  \item{sigma0.hat}{  }
  \item{D}{  }
  \item{oll}{  }
  
  
}
\references{
Ma, W., Chen, L.S., Han, S.W., Zhong, H., Wang, P., and CPTAC consortium. (2015) \emph{A Regularized Linear Mixed Effect Model (Mixed remMap) with application on the iTRAQ experiment data}

}
\author{
Weiping Ma, Sung Won Han, Hua Zhong, and Pei Wang

Maintainer: Pei Wang <pei.wang@mssm.edu>
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
# Extract Data

library(mvtnorm);
library(corpcor);
library(MASS);
library(remMap);
library(Matrix);
library(impute);
library(ProMap);
          
data(cna77)
data(raw_pho770.3)

# Change data format
X_array = array(as.matrix(cna77),dim = c(36,4,dim(cna77)[2]));
Y_array = array(as.matrix(raw_pho770.3),dim = c(36,4,dim(raw_pho770.3)[2]));

L1 = c(11:20);
L2 = seq(0,90,by =10);
Lam = expand.grid(L1,L2);

# Job number
ii = 320 ;

ProMap.result <- ProMap(Xarray=X_array, Yarray=Y_array, lambda.1=Lam[ceiling(ii/4),1], lambda.2=Lam[ceiling(ii/4),2], nu.1=0.01, nu.2=1, CV.i=floor((ii-1)\%\%(4)), CV.all = 4, max.teration =10)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
