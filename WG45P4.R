# Group 45, Pavel Bozmarov, Siying Zhu, Xuexun Lu
# https://github.com/Louhokseson/SP_coursework4_Xuexun_Pavel_Shiying.git
#
#                                     ***** OVERVIEW *****
# 
# 
#
# 
#                                            *****  

#########################################################################################################################
Hessian = function(f,theta,...)
{
  #                               ***** DESCRIPTION *****
  #
  # what type of function is this
  # objects of the function
  #
  # How the function works.
  #                                     
  #                                       ***** 
  
  n = length(theta)
  eps <- 1e-7 ## finite difference interval
  Hfd <- matrix(0,n,n) ## finite diference Hessian for (i in 1:length(th0)) { ## loop over parameters
  th0 = theta
  gll0 = gradient(f,th0,...)
  for (i in 1:length(th0)) { ## loop over parameters
    th1 <- th0; th1[i] <- th1[i] + eps ## increase th0[i] by eps gll1 <- gll(th1,t=t80,y=y) ## compute resulting nll
    gll1 = gradient(f,th1,...)
    Hfd[i,] <- (gll1 - gll0)/eps ## approximate second derivs
  }
return(0.5*(t(Hfd)+Hfd))
}

#########################################################################################################################
rb <- function(theta,getg=FALSE,k=10) {
  ## Rosenbrock objective function, suitable for use by ’bfgs’
  z <- theta[1]; x <- theta[2]
  f <- k*(z-x^2)^2 + (1-x)^2 + 1
  if (getg) {
    attr(f,"gradient") <- c(2*k*(z-x^2),
                            -4*k*x*(z-x^2) -2*(1-x))
  }
  f
} ## rb

#########################################################################################################################
gradient = function(f,theta,...)
{
  #                               ***** DESCRIPTION *****
  #
  # what type of function is this
  # objects of the function
  #
  # How the function works.
  #                                     
  #                                       ***** 
  
  function_value = f(theta,...)
  if (length(attr(function_value,'gradient'))>0)
  {
    return(attr(function_value,'gradient'))
  }
  else
  {
    return(numerical_gradient(f,theta,...))    
  }
}

#########################################################################################################################
numerical_gradient = function(f,theta,...)
{
  #                               ***** DESCRIPTION *****
  #
  # what type of function is this
  # objects of the function
  #
  # How the function works.
  #                                     
  #                                       ***** 
  
  fd <- th0 <- theta ## test param value, and approx grad vec nll0 <- nll(th0,t=t80,y=y) ## nll at th0
  f0 <- f(th0,...) ## nll at th0
  eps <- 1e-5 ## finite difference interval
  for (i in 1:length(th0)) 
  { ## loop over parameters
    
    th1 <- th0; th1[i] <- th1[i] + eps ## increase th0[i] by eps nll1 <- nll(th1,t=t80,y=y) ## compute resulting nll
    f1 = f(th1,...)
    fd[i] <- (f1 - f0)/eps ## approximate -dl/dth[i]
    
  }
  return(drop(fd))
}

#########################################################################################################################
our_bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=140){
  
  #                               ***** DESCRIPTION *****
  #
  # what type of function is this
  # objects of the function
  #
  # How the function works.
  #                                     
  #                                       ***** 
  
  number_par <- length(theta) # the number of parameters
  
  current_theta <- theta
  theta_to_return = rep(0,number_par)
  # Initialization of B as a identity matrix and let I be a identity matrix as well
  
  current_B <- I <- diag(number_par)
  
  N_steps = 0
  
  while (N_steps<=maxit)
  {

    current_function_value = f(current_theta,...)
    
    current_gradient = gradient(f,current_theta,...)
   
    #current_step
    QN_step = - current_B %*% current_gradient
    
    updated_theta <- current_theta + QN_step
    updated_function_value = f(updated_theta,...)
    
    updated_gradient = gradient(f,updated_theta,...)
    
    
    #step_check with condition 2, if not satisfied , reduce by half
    
    counter = 0 
    while ( f(updated_theta,...)>=f(current_theta,...) )
    {
      
      QN_step <- QN_step/2
      updated_theta <- current_theta + QN_step
      updated_function_value = f(updated_theta,...)
      
      updated_gradient = gradient(f,updated_theta,...)
      
      counter = counter + 1
      #print("counter")
      #print(counter)
      
    }
    
    while(updated_gradient %*% QN_step < 0.9 * current_gradient %*% QN_step)
    {
      QN_step <- QN_step + 0.1* QN_step
      updated_theta <- current_theta + QN_step
      updated_function_value = f(updated_theta,...)
      updated_gradient = gradient(f,updated_theta,...)
      
      if(updated_function_value >= current_function_value)
      {
        stop("ERROR")
      
      }
    
    }
    
    ########################## if our step is ok, update B ##################################
    # s and y
    s <- updated_theta - current_theta; y <- updated_gradient-current_gradient
    s=drop(s)
    #rho
    rho <- (1/(t(s) %*% y))[1]
    updated_B = (I-rho*(s%*%t(y)))%*%current_B%*%(I-rho*y%*%t(s)) + rho*s%*%t(s)

    
    ########################## move to next theta and B ##################################
    theta_to_return = current_theta
    current_theta = updated_theta
    current_B = updated_B
    
    N_steps = N_steps+1
    
    if( max(abs(current_gradient)) < (abs(current_function_value)+fscale)*tol)
        {
          print("Optimal Solution Found")
          break
        }
   
  }
  
  H = Hessian(f,theta_to_return,...)
  return(list(f=current_function_value,theta=theta_to_return,iter = N_steps, g = current_gradient ,H = H))
  
}

res=our_bfgs(theta=c(-1,2),rb,getg=FALSE)
print(res)
