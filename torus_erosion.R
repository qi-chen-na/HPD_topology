
library("rgl")
source('funcode.R')
source('erosion_new.R')
##test for torus


#Simulate points in a 3d torus, points are more denser around the circle in the middle
n <- 100000
r <- 0.1
R <- 1
torus_data <- function(n,r,R){
  sigma <- (R-r)/4  #so that the boudnary is at 3 sigma
  theta1 <- runif(n,min = 0, max = 2*pi)
  theta2 <- runif(n, min = 0, max = 2*pi)
  norm_points <- abs(rnorm(n, mean = 0, sd = sigma^2))
  x <- ((r + R)/2 - norm_points*cos(theta2))*cos(theta1)
  y <- ((r + R)/2 - norm_points*cos(theta2))*sin(theta1)
  z <- sin(theta2)*norm_points
  torus_data <- array(data = NA,dim = c(n,3))
  torus_data[,1] <- x
  torus_data[,2] <- y
  torus_data[,3] <- z
  return (torus_data)
}

torus_points1 <- torus_data(n,r,R)
#Interactive plot of the simulated torus
plot3d(torus_points1[,1],torus_points1[,2],torus_points1[,3],xlim = c(-1.3,1.3),ylim = c(-1.3,1.3), zlim = c(-1.3,1.3),xlab = 'x', ylab = 'y',zlab = 'z')
torus_points2 <- torus_data(n,r,R)

#Implement HPD sets estimation
train.sim <- torus_points1
test.sim <- torus_points2
conv_training <- convert.marginal(train.sim,train.sim)
conv_testing <- convert.marginal(train.sim,test.sim)
vector_tau = seq(0.5,0.01,len=10)
alpha = 0.5
out=cred.est2(test=conv_testing,train.sim=train.sim,train=conv_training,tau=vector_tau,alpha=alpha,p.alpha=0.95)

k <- out$n.tree
HPD <- out$cre.set
all_cells <- out$uni_space
complement <- all_cells[(length(HPD)+1):k]
e <- erosion1(HPD, complement)
graph <- e[[1]]
HPD_vertices <- e[[2]][[1]]
rglplot(induced_subgraph(graph,HPD_vertices), vertex.size = 1)

plot(torus_points1)
xx <- c(1,2)
f1=range(torus_points1[,xx[1]]); f2=range(torus_points1[,xx[2]])
for (i in 1:length(HPD)){
  pt=all_cells[[i]][,xx]
  pt[,1]=pt[,1]*diff(f1)+f1[1]
  pt[,2]=pt[,2]*diff(f2)+f2[1]
  xleft <- pt[1,1]; ybottom <- pt[1,2]; xright <- pt[2,1]; ytop <- pt[2,2]
  rect(xleft, ybottom, xright, ytop, density = 11, angle = 45, col = 'red')
}
for (i in HPD_vertices){
  pt=all_cells[[i]][,xx]
  pt[,1]=pt[,1]*diff(f1)+f1[1]
  pt[,2]=pt[,2]*diff(f2)+f2[1]
  xleft <- pt[1,1]; ybottom <- pt[1,2]; xright <- pt[2,1]; ytop <- pt[2,2]
  rect(xleft, ybottom, xright, ytop, density = 11, angle = 45, col = 'green')
}


