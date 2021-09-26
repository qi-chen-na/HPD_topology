library("rgl")
source('funcode.R')
source('erosion_new.R')



#Simulate points in a 3d torus, points are more denser around the circle in the middle
n <- 100000
R <- 1
#In this simulation we deliberately make the points denser around the center uniformly in the sense that the points are uniformly distributed on the every sphere surface
ball <- function(n,R){
  sigma <- R/2
  theta <- runif(n,min = 0, max = 2*pi)
  phi <- acos(1-2*runif(n,min = 0, max = 1))
  r <- abs(rnorm(n,0,sigma^2))
  x <- r * sin(phi)*cos(theta)
  y <- r * sin(phi)*sin(theta)
  z <- r * cos(phi)
  points <- array(data = NA,dim = c(n,3))
  points[,1] <- x
  points[,2] <- y
  points[,3] <- z
  return (points)
}
ball_points1 <- ball(n,R)
plot3d(ball_points1[,1],ball_points1[,2],ball_points1[,3])
ball_points2 <- ball(n,R)


train.sim <- ball_points1
test.sim <- ball_points2
conv_training <- convert.marginal(train.sim,train.sim)
conv_testing <- convert.marginal(train.sim,test.sim)
vector_tau = seq(0.5,0.01,len=10)
alpha = 0.1
out=cred.est2(test=conv_testing,train.sim=train.sim,train=conv_training,tau=vector_tau,alpha=alpha,p.alpha=0.95)

k <- out$n.tree
HPD <- out$cre.set
all_cells <- out$uni_space
complement <- all_cells[(length(HPD)+1):k]
e <- erosion1(HPD, complement)
graph <- e[[1]]
HPD_vertices <- e[[2]][[1]]
rglplot(induced_subgraph(graph,HPD_vertices), vertex.size = 1)



