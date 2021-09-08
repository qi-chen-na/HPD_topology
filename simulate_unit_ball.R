##Simulate points uniformly distributed in d-dimensional ball
n <- 10000
p <- 10


i <- 0
simulate <- array(dim = c(n,p))
while (i < n){
  point <- runif(p)
  if (sqrt(sum(point ** 2)) <= 1){
    simulate[i + 1,] <- point
    i <- i + 1
  }
}

simulate
