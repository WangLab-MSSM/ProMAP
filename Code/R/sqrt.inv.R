sqrt.inv <-
function(PDMatr)
{  
  e <- eigen(PDMatr);
  S <- e$vectors;
  ###  S %*% diag(e$values) %*% t(S) == PDMatr;
  B <- S %*% diag(sqrt(e$values)) %*% t(S);
  B.inv = my.solve(B); 
  return(B.inv);
}
