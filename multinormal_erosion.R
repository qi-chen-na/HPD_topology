library(MASS)
library("rgl")
source('funcode.R')
source('erosion_new.R')

#Simulate multivariate normal distributed points 
n <- 100000
sigma <- matrix(c(10,3,3,2),2,2)
mu <- c(0,0)
mvrn_points1 <- mvrnorm(n,mu,Sigma = sigma)
# dev.new()
plot(mvrn_points1, xlab = 'x', ylab = 'y')
mvrn_points2 <- mvrnorm(n,mu,Sigma = sigma)

#Implement HPD sets estimation
train.sim <- mvrn_points1
test.sim <- mvrn_points2
conv_training <- convert.marginal(train.sim,train.sim)
conv_testing <- convert.marginal(train.sim,test.sim)
vector_tau = seq(0.5,0.01,len=10)
alpha = 0.1
out=cred.est2(test=conv_testing,train.sim=train.sim,train=conv_training,tau=vector_tau,alpha=alpha,p.alpha=0.95)

k <- out$n.tree
HPD <- out$cre.set
all_cells <- out$uni_space
complement <- all_cells[(length(HPD)+1):k]

##Visualization of the estimated HPD sets
#Plot HPD sets in red
xx <- c(1,2)
f1=range(mvrn_points1[,xx[1]]); f2=range(mvrn_points1[,xx[2]])
for (i in 1:length(HPD)){
  pt=all_cells[[i]][,xx]
  pt[,1]=pt[,1]*diff(f1)+f1[1]
  pt[,2]=pt[,2]*diff(f2)+f2[1]
  xleft <- pt[1,1]; ybottom <- pt[1,2]; xright <- pt[2,1]; ytop <- pt[2,2]
  rect(xleft, ybottom, xright, ytop, density = 20, angle = 45, col = 'red')
}
#Plot complement sets in blue
for (i in (length(HPD)+1):k){
  pt=all_cells[[i]][,xx]
  pt[,1]=pt[,1]*diff(f1)+f1[1]
  pt[,2]=pt[,2]*diff(f2)+f2[1]
  xleft <- pt[1,1]; ybottom <- pt[1,2]; xright <- pt[2,1]; ytop <- pt[2,2]
  rect(xleft, ybottom, xright, ytop, density = 20, angle = 45, col = 'blue')
}

##Apply erosion method on the multivariate normal data
e <- erosion2(HPD, complement)
all_graph <- e[[1]]
HPD_vertices <- e[[2]][[1]]

#plot the HPD set remains
for (i in HPD_vertices){
  pt=all_cells[[i]][,xx]
  pt[,1]=pt[,1]*diff(f1)+f1[1]
  pt[,2]=pt[,2]*diff(f2)+f2[1]
  xleft <- pt[1,1]; ybottom <- pt[1,2]; xright <- pt[2,1]; ytop <- pt[2,2]
  rect(xleft, ybottom, xright, ytop, density = 30, angle = 45, col = 'green')
}

#plot the graph of HPD set remains
plot(induced_subgraph(all_graph, HPD_vertices))

