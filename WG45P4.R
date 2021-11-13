




#test function 
rb <- function(theta,getg=FALSE,k=10) {
  ## Rosenbrock objective function, suitable for use by 'bfgs'
  z <- theta[1]; x <- theta[2]
  f <- k*(z-x^2)^2 + (1-x)^2 + 1
  if (getg) {
    attr(f,"gradient") <- c(2*k*(z-x^2),
                            -4*k*x*(z-x^2) -2*(1-x))
  }
  f
} ## rb

# Hokseson's practice
myfunction <- function(theta,a,b) {
 a+b+theta
} ## rb

Test_bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100){
  # Here, we need use the three dots to input our variable into the funciton f
  result <- f(theta,...) 
  return (result)
}


Test_bfgs(theta=3,myfunction,a=5,b=6)

fit <- optim(c(-1,2),rb,getg = FALSE,hessian=TRUE)

bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100){
  
  number_par <- length(theta) # the number of parameters
  # BFGS method
  # Initialization of B as a identity matrix
  B <- diag(number_par)
  QN_step = first_derivative(f)
  return (result)
}
