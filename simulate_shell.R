##simulate shells in p-dimension
n <- 10000
p <- 10
lower_bound <- 0.5
upper_bound <- 1
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
simulate
