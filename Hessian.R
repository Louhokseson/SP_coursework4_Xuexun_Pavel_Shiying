
numerical_hessian = function(f,theta,...)
{
  th0 <- theta
  eps <- 1e-7 ## finite difference interval
  number_par <-length(th0)
  fh <- matrix(0,length(th0),length(th0))


  for (i in 1:number_par) ## loop over rows in hessain matrix
  { 
    for (j in i:number_par) ## loop over upper triangle of hessain
    {
      Add_vector <- rep(0, number_par); Add_vector[j] <- eps
      fh[i,j] <- (numerical_gradient(f,theta+Add_vector)[i] - numerical_gradient(f,theta)[i])/eps
      fh[j,i] <- fh[i,j] ## 
    }
  }
  return(fh)
}

numerical_gradient = function(f,theta,...)
{
  fd <- th0 <- theta ## test param value, and approx grad vec nll0 <- nll(th0,t=t80,y=y) ## nll at th0
  f0 <- f(th0,...) ## nll at th0
  eps <- 1e-7 ## finite difference interval
  for (i in 1:length(th0)) 
  { ## loop over parameters
    
    th1 <- th0; th1[i] <- th1[i] + eps ## increase th0[i] by eps nll1 <- nll(th1,t=t80,y=y) ## compute resulting nll
    f1 = f(th1,...)
    fd[i] <- (f1 - f0)/eps ## approximate -dl/dth[i]
    
  }
  return(fd)
}