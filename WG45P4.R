# Group 45, Pavel Bozmarov, Siying Zhu, Xuexun Lu
# https://github.com/Louhokseson/SP_coursework4_Xuexun_Pavel_Shiying.git
#
#                                     ***** OVERVIEW *****
#     
# This program creates a function that implements the BFGS quasi-Newton minimization method.
# The BFGS algorithm is a quasi-Newton method for solving unconstrained optimization problems.
# BFGS determines the descent direction by preconditioning the gradient with curvature information. 
# It does so by gradually improving an approximation to the Hessian matrix of the observable 
# function, obtained only from gradient evaluations.
# 
#
#                                           *****  


#########################################################################################################################
bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=140){
  
  
  #                               ***** DESCRIPTION *****
  #
  # bfgs is a function that implements the BFGS quasi-Newton minimization method. 
  # 
  # Steps:
  #  1) setting initial theta, initial function value ( the function value at the initial theta) and initial gradient.
  #  2) check wheter our initial function value or initial gradient are finate.
  #  3) We loop until we find an optimal solution - reduce our current step by half until we reduce the function
  #                                               - check wolfe second condition
  #                                               - check for convergence
  #                                               - raise errors when needed
  # 
  # Input arguments: 
  #   1) f      - the objective function to minimize. Its first argument is the vector of optimization parameters. 
  #               Its second argument is a logical indicating whether or not gradients of the objective w.r.t. the 
  #               parameters should be computed.Remaining arguments will be passed from bfgs using ‘...’. The scalar 
  #               value returned by f will have a gradient attribute if the second argument to f is TRUE.
  #   
  #   2) theta  - vector of initial values for the optimization parameters.
  #   3) tol    - the convergence tolerance.
  #   4) fscale - a rough estimate of the magnitude of f at the optimum - used in convergence testing.
  #   5) maxit  - the maximum number of BFGS iterations to try before giving up.    
  # 
  # Outputs a list containing the following arguments:
  #   1) f      -  the scalar value of the objective function at the minimum.
  #   2) theta  -  the vector of values of the parameters at the minimum.
  #   3) iter   -  the number of iterations taken to reach the minimum.
  #   4) g      -  the gradient vector at the minimum (so the user can judge closeness to numerical zero).
  #   5) H      -  the approximate Hessian matrix (obtained by finite differencing) at the minimum.
  
  number_par <- length(theta) # the number of parameters in our theta
  
  
  current_theta <- theta # the initial theta
  
  # Initialization of B, the approximate inverse Hessian Matrix, 
  # as an Identity matrix and let I be the identity matrix as well
  current_B <- I <- diag(number_par)
  
  # counter of iterations
  N_steps = 0
  
  current_function_value = f(current_theta,...) # the function value at the initial theta
  
  current_gradient = gradient(f,current_theta,...) # the gradient of the function at the initial theta
  
  # WARNING. We check if  NaN, Inf or -Inf ocurred in in our objective function or in our gradient at yout current theta
  if (is.finite(current_function_value)==FALSE | sum(is.finite(current_gradient)) < number_par)
  {
    stop('Error: the objective or derivatives are not finite at your initial parameters. Please have a check.')
  }
  
  # We loop until out iterations surpass maxit
  while (N_steps<=maxit)
  {
    
    current_function_value = f(current_theta,...) # the function value at the theta of the current step
    
    current_gradient = gradient(f,current_theta,...) # the gradient of the function at theta of the current step
    
    #current_step
    QN_step = - current_B %*% current_gradient
    
    updated_theta <- current_theta + QN_step # our next (updated) theta
    updated_function_value = f(updated_theta,...) # the updated function value - the function value at the updated theta
     
    updated_gradient = gradient(f,updated_theta,...) # the gradient of the function at the updated theta
    
    counter = 0 # counter to represent number in iterations in the below loop. 
    
    # We reduce our step by half until we find a theta at which the function has a smaller value than the function
    # value at the current theta. If we can't find such a theta until looping 100 times, we check for convergence
    while ( f(updated_theta,...)>=f(current_theta,...) )
    {
      
      QN_step <- QN_step/2 # update step
      updated_theta <- current_theta + QN_step # update the current theta
      updated_function_value = f(updated_theta,...) # update the function value
       
      updated_gradient = gradient(f,updated_theta,...) # update the gradient
      
      counter = counter + 1 # update the counter
      
      # if we haven't decreased our function within 100 steps, we stop.
      if (counter == 100 & max(abs(current_gradient)) >= (abs(current_function_value)+fscale)*tol)
      {
        stop('Error:The step fails to reduce the objective but convergence has not occurred')
      }
      
    }
    
    # We check the second Wolfe condition with c2 = 0.9 if not we start enlarging our current step
    # a little bit until we satisfy the Wolfe condition or until we make the observable function bigger.
    while(updated_gradient %*% QN_step < 0.9 * current_gradient %*% QN_step)
    {
      QN_step <- QN_step + 0.1* QN_step # we enlarge our step with 10 percent of the currenst step
      updated_theta <- current_theta + QN_step  # update theta
      updated_function_value = f(updated_theta,...) # update the function value
      updated_gradient = gradient(f,updated_theta,...) # update the gradient
      
      # WARNING. If we our function at the updated step becomes larger than the value of the function
      # at the step that we started the while loop with, then we STOP!
      if(updated_function_value >= current_function_value)
      {
        stop("Error: Step unable to reduce observable function")
        
      }
      
    }
    
    # s - difference between the theta at the current step and the theta at the next step
    # y - difference between the gradient at the current step and the gradient at the next step
    # 
    s <- updated_theta - current_theta; y <- updated_gradient-current_gradient
    s=drop(s) # our updated theta is obtained by matrix vector multiplication, so in R it will be stored as a 
              # matrix and we need to convert it to vector.
    
    
    rho <- (1/(t(s) %*% y))[1]
    
    By = drop(current_B%*%y)
    
    #calculate the approximation of the updated inverse Hessian matrix
    updated_B = current_B - rho * ( By%*%t(s) + s%*%(t(y)%*%current_B) ) + rho^2* ( ((s%*%t(y))%*%By)%*%t(s)) +  rho*s%*%t(s)

    
   # CHECK FOR OPTIMAL SOLUTION AT CURRENT ITERATION  
    if( max(abs(current_gradient)) < (abs(current_function_value)+fscale)*tol)
    {
      print("Optimal Solution Found")
      
      H = Hessian(f,current_theta...) # Calculate approximation of the Hessian at the current theta
      
      return(list(f=current_function_value,theta=current_theta,iter = N_steps, g = current_gradient ,H = H))
    }
    
    # Make our current_theta and current_B to be the updated_theta and the updated_B in order to move to the next
    # iteration appropriately.
    current_theta = updated_theta
    current_B = updated_B
    
    N_steps = N_steps+1 #update iteration
    

    
  }
    
  # WARNING maximum iteration is reached without convergence
  warning("Error: The algorithm reached the maximum iteration but it didn't converge.")

}



#########################################################################################################################
Hessian = function(f,theta,...)
{
  #                               ***** DESCRIPTION *****
  #
  # This function approximates the Hessian of a given function at a given point using finite differencing.
  # 
  # Input arguments:
  #   1) f     - function
  #   2) theta - point
  #
  # Ouput:
  #   H - numerical approximation of the Hessian of f at the point theta.
  #                                     
  #                                       ***** 
  
  n = length(theta) # the length of theta
  
  eps <- 1e-7 ## finite difference interval
  
  th0 = theta # th0 - the point that we want to approximate the Hessian at.
  
  Hfd <- matrix(0,n,n) ## finite diference Hessian for (i in 1:length(th0))  ## loop over parameters
  
  gll0 = gradient(f,th0,...) # taking the gradient of t at point th0
  
  for (i in 1:length(th0)) { ## loop over parameters
    th1 <- th0; th1[i] <- th1[i] + eps ## increase th0[i] by eps
    gll1 = gradient(f,th1,...) # gradient of f at th1
    Hfd[i,] <- (gll1 - gll0)/eps ## approximate second derivs
  }
  return(0.5*(t(Hfd)+Hfd)) # make our Hessian to be symmetric and return it
}

#########################################################################################################################
gradient = function(f,theta,...)
{
  #                               ***** DESCRIPTION *****
  #
  # This function checks whether a function has an attribute 'gradient' inside its implementation. If yes, return this
  # gradient, if no, compute the numerical approximation of the gradient at point theta.
  # 
  # Input arguments:
  #   1) f     - function
  #   2) theta - point
  #
  # Ouput:
  #   gradient or numerical gradient of f.
  #                                     
  #                                       ****  
  
  function_value = f(theta,...) # take the function value at point theta
  
  # if our function_value has attribute gradient, then return it
  if (length(attr(function_value,'gradient'))>0)
  {
    return(attr(function_value,'gradient'))
  }
  else # if not, calculate it using our function numerical_gradient
  {
    return(numerical_gradient(f,theta,...))    
  }
}

#########################################################################################################################
numerical_gradient = function(f,theta,...)
{
  
  
  #                               ***** DESCRIPTION *****
  #
  # This function approximates the gradient of a given function at a given point using finite differencing.
  # 
  # Input arguments:
  #   1) f     - function
  #   2) theta - point
  #
  # Ouput:
  #   fd - numerical approximation of the gradient of f at the point theta.
  #                                     
  #                                       ****  
  
  fd <- th0 <- theta ## fd - our gradient approximation, th0 - the point that we want to approximate the gradient at.
  f0 <- f(th0,...) ## f at th0
  eps <- 1e-5 ## finite difference interval
  for (i in 1:length(th0))  # loop through the elements of theta
  { 
    
    th1 <- th0; th1[i] <- th1[i] + eps ## increase th0[i] by eps 
    f1 = f(th1,...) # f at th1
    fd[i] <- (f1 - f0)/eps ## approximate the i-th component of the gradient
    
  }
  return(drop(fd)) # return the gradient
}

########################################  END OF THE IMPLEMENTATION PART ############################################################################
########################################  * * * * * * * * * * * * * * *  ###########################################################################



#######################################  * * *  TESTING PART  * * *  ##############################################################

=======
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
>>>>>>> 9ce900e155d14c01742dbe46f40779f44752e943
rb <- function(theta,getg=FALSE,k=10) {
  ## Rosenbrock objective function, suitable for use by ’bfgs’
  z <- theta[1]; x <- theta[2]
  f <- k*(z-x^2)^2 + (1-x)^2 + 1
  if (getg) {
    attr(f,"gradient") <- c(2*k*(z-x^2),
                            -4*k*x*(z-x^2) -2*(1-x))
  }
  f
<<<<<<< HEAD
} 

res=bfgs(c(-1,2),rb,getg=TRUE)
print(res)
=======
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
  #browser()
  number_par <- length(theta) # the number of parameters
  
  current_theta <- theta
  theta_to_return = rep(0,number_par)
  # Initialization of B as a identity matrix and let I be a identity matrix as well
  
  current_B <- I <- diag(number_par)
  
  N_steps = 0
  
  current_function_value = f(current_theta,...)
  
  current_gradient = gradient(f,current_theta,...)
  
  # WARNING infinite WARNING
  if (is.infinite(current_function_value)|sum(is.infinite(current_gradient))!=0)
  {
    stop('Error: the objective or derivatives are not finite at your initial parameters. Please have a check.')
  }
  
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
      # WARNING NOT CONVERGED WHEN IT HIT THE MAXIMUM ITERATION WARNING
      if (counter > 100 & max(abs(current_gradient)) > (abs(current_function_value)+fscale)*tol)
      {
        stop('Error:The step fails to reduce the objective but convergence has not occurred')
      }
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
        stop("Error: the update function value is worse than the previous one in Wolfe condition 2")
      
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
  # WARNING maximum iteration is reached. WARNING
  if( max(abs(current_gradient)) > (abs(current_function_value)+fscale)*tol)
  {
    stop("Error: The algorithm reached the maximum iteration but it didn't converge.")
  }
  
  H = Hessian(f,theta_to_return,...)
  return(list(f=current_function_value,theta=theta_to_return,iter = N_steps, g = current_gradient ,H = H))
  
}

res=our_bfgs(theta=c(-1,2),rb,getg=FALSE)
print(res)
>>>>>>> 9ce900e155d14c01742dbe46f40779f44752e943
