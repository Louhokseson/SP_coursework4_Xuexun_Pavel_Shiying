pos_def_test = function(A)
{
  eigen_vals = eigen(A)$values
  
  if(min(eigen_vals) >= 0)
  {
    
    return("Positive Definite")
  }
  else
  {
    return("Negative Definite")  
    
  }
  
}

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

numerical_gradient = function(f,theta,...)
{
  fd <- th0 <- theta ## test param value, and approx grad vec nll0 <- nll(th0,t=t80,y=y) ## nll at th0
  f0 <- f(th0,...) ## nll at th0
  eps <- 1e-7 ## finite difference interval
  for (i in 1:length(th0)) 
  { ## loop over parameters
    
    th1 <- th0; th1[i] <- th1[i] + eps ## increase th0[i] by eps nll1 <- nll(th1,t=t80,y=y) ## compute resulting nll
    f1 = f(th1,...)
    fd[i] <- (nll1 - nll0)/eps ## approximate -dl/dth[i]
    
  }
  return(fd)
}


our_bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=140){
  
  #
  #browser()
  
  number_par <- length(theta) # the number of parameters
  
  current_theta <- theta
  
  # Initialization of B as a identity matrix and let I be a identity matrix as well
  
  current_B <- I <- diag(number_par)
  
  N_steps = 0
  
  while (N_steps<=maxit)
  {
    wolfey = "WOLFEY 2 MET"
    
    current_function_value = f(current_theta,...)
    
    if(is.null(attr(current_function_value,"gradient")))
    {
      current_gradient = numerical_gradient(f,current_theta,...)
      
    }
    else
    {
      current_gradient = attr(current_function_value,"gradient")
    }
    #current_step
    QN_step = - current_B %*% current_gradient
    
    updated_theta <- current_theta + QN_step
    updated_function_value = f(updated_theta,...)
    
    if(is.null(attr(updated_function_value,"gradient")))
    {
      updated_gradient = numerical_gradient(f,updated_theta,...)
      
    }
    else
    {
      updated_gradient = attr(updated_function_value,"gradient")
    }
    
    #check if our current step is ok
    # wolfe conditions focus on the second one
    
    
    c_2 <- 0.9
    print(N_steps)
    print(pos_def_test(current_B))
    print(current_B)
    #step_check with condition 2, if not satisfied , reduce by half
    
    counter = 0 
    while ( f(updated_theta,...)>=f(current_theta,...) )
    {
      
      QN_step <- QN_step/2
      updated_theta <- current_theta + QN_step
      updated_function_value = f(updated_theta,...)

      
      if(is.null(attr(updated_function_value,"gradient")))
      {
        updated_gradient = numerical_gradient(f,updated_theta,...)
        
      }
      else
      {
        updated_gradient = attr(updated_function_value,"gradient")
      }
      counter = counter + 1
      #print("counter")
      #print(counter)
      
    }
    
    while(updated_gradient %*% QN_step < c_2 * current_gradient %*% QN_step)
    {
      QN_step <- QN_step + 0.1* QN_step
      updated_theta <- current_theta + QN_step
      updated_function_value = f(updated_theta,...)
      
      
      if(is.null(attr(updated_function_value,"gradient")))
      {
        updated_gradient = numerical_gradient(f,updated_theta,...)
        
      }
      else
      {
        updated_gradient = attr(updated_function_value,"gradient")
      }
      
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
    ########################################################################################
    
    print("updated_B")
    print(updated_B)
    
    ########################## move to next theta and B ##################################
    current_theta = updated_theta
    current_B = updated_B
    
    
    #update step number
    N_steps=N_steps+1
    #print(wolfey)
    print("step")
    print(QN_step)
    print("updated_B")
    print(current_B)
    print("Engaging Next Step")
    print(" ")
    
    if( max(abs(current_gradient)) < (abs(current_function_value)+fscale)*tol)
        {
          print("Optimal Solution Found")
          break
        }
    {
      
      
    }
    }
  
  
  return(list(updated_function_value=updated_function_value,current_theta=current_theta))
  
}

res=our_bfgs(theta=c(-10,100),rb,getg=TRUE)
print(res)
