# Hokseson's practice
myfunction <- function(theta,a,b,i=FALSE) {
  
  if(i){
    print(1)
  }
  a+b+theta
} ## rb




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

Test_bfgs <- function(theta,f,...){
  # Here, we need use the three dots to input our variable into the funciton f
  result <- f(theta,...,getg = TRUE) 
  return (result)
}


fit <- Test_bfgs(theta=c(1,2),rb,k=10)

fit <- optim(c(-1,2),rb,getg = FALSE,hessian=TRUE)

bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100){
  browser()
  number_par <- length(theta) # the number of parameters
  current_position <- theta
  current_point <- f(current_position,...,getg = TRUE)
  # BFGS method
  # Initialization of B as a identity matrix and let I be a identity matrix as well
  B <- diag(number_par)
  I <- diag(number_par)
  # call the getg = TRUE when we need to compute the gradient
  QN_step <- -B%*%attr(current_point,"gradient")
  next_position <- current_position + QN_step
  next_point <- f(next_position,...,getg = TRUE)
  
  # wolfe conditions focus on the second one
  c_1 <- 0.2
  c_2 <- 0.9
  while (attr(next_point,"gradient")%*%QN_step < c_2 * attr(current_point,"gradient")%*%QN_step){
    QN_step <- QN_step/2
    next_position <- current_position + QN_step
    next_point <- f(next_position,...,getg = TRUE)
  }
  
  # s and y
  s <- next_position - current_position; y <- attr(next_point,"gradient")-attr(current_point,"gradient")
  # rho
  rho <- s %*% y
  # update the B

  
}
fit_bfgs <- bfgs(theta=c(5,7),rb)
