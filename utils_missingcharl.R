# estimating the parameters for the Charlson missingness model (mm):

# inputs: vector V of n_visits, vector Yobs of 1 for 'code', 0 for no code (0 or NA
# in original data)

# theta = (p,phi) 

estimate_mm <- function(Y_obs,V){
  # initialize
  theta_old <- init_mm(Y_obs,V)
  # optimization parameters
  eps <- 1e-5
  Kmax <- 100
  # loop
  for(kk in 1:Kmax){
    theta_new <- newtom_mm(theta_old,Y_obs,V)
    if(mean((theta_old - theta_new)^2) < eps){
      break
    }
  }
  # return
  return(theta_new)
}

init_mm <- function(Y_obs,V,K=400){
  p_init <- mean(Y_obs[V > K])
  phi_init <- mean(Y_obs[V==1]) / p_init
  return(c(p_init,phi_init))
}

newton_mm <- function(theta,Y_obs,V){
  p <- theta[1]
  phi <- theta[2]
  # sub components of gradient and Hessian
  c1 <- 1 - (1-phi)^V
  c2 <- V*((1-phi)^(V-1))
  c3 <- V*(V-1)*((1-phi)^(V-2))
  # components of gradient
  g1 <- sum((1/p)*(Y_obs==1)) + sum(((-c1)/(1-p*c1))*(Y_obs==0))
  g2 <- sum((c2/c1)*(Y_obs==1)) + sum(((-p*c2)/(1-p*c1))*(Y_obs==0))
  grad <- c(g1,g2)
  # components of Hessian
  H11 <- sum(((-1)/(p^2))*(Y_obs==1)) + sum(((-(c1)^2)/((1-p*c1)^2))*(Y_obs==0))
  H12 <- sum((c2/((1-(p^2)*c1)^2))*(Y_obs==0))
  H22 <- sum(((-c3*c1 - c2^2)/(c1^2))*(Y_obs==1)) + sum((((-c3*((1/p) - c1)) - c2^2)/(((1/p) - c1)^2))*(Y_obs==0))
  hess <- matrix(c(H11,H12,H12,H22),2,2)
  # newton step
  theta_new <- theta - solve(hess,grad)
  return(theta_new)
}

