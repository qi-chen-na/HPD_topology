###simulation
library(alphashape3d)
library(pracma)
##simulate shells in p-dimension uniformly

shell_simulate <- function(n,p,lower_bound = 0.6, upper_bound = 1){
  i <- 0
  simulate <- array(dim = c(n, p))
  while (i < n){
    point <- runif(p)
    point_norm <- sqrt(sum(point ** 2))
    if (lower_bound < point_norm & point_norm < upper_bound){
      simulate[i+1,] <- point
      i <- i+1
    }
  }
  return (simulate)
}

##Simulate points uniformly distributed in p-dimensional ball

ball_simulate <- function(n,p){
  i <- 0
  simulate <- array(dim = c(n,p))
  while (i < n){
    point <- runif(p)
    if (sqrt(sum(point ** 2)) <= 1){
      simulate[i + 1,] <- point
      i <- i + 1
    }
  }
  return (simulate)
}

##Simulate points uniformly distributed in torus in 3d and then extend to p-dimension with orthonormal transformation
torus_simulate <- function(n, p, R = 1, r = 0.6, origin = c(0, 0, 0)){
  points_torus <- rtorus(n, lower_bound, upper_bound, ct= origin)
  simulate_before <- cbind(array(points_torus,dim = c(n,3)),array(0,dim = c(n,p-3)))
  rand_ortho_matrix <- randortho(p, type = 'orthonormal')
  simulate_after <- simulate_before %*% rand_ortho_matrix
  return (simulate_after)
}

shell_poits <- shell_simulate(n,p)
ball_points <- ball_simulate(n,p)
torus_points <- torus_simulate(n,p)

source('funcode.R')


