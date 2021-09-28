library("rgl")
source('funcode.R')
source('erosion_new.R')

##Simulate a "eight-shape" region
n <- 100000
r <- 0.1
R <- 1
eight <- function(n,r,R){
  sigma <- (R-r)/4
  theta1 <- runif(n/2,min = 0, max = 2*pi)
  theta2 <- runif(n/2, min = 0, max = 2*pi)
  n1 <- rnorm(n/2,0,sigma^2)
  n2 <- rnorm(n/2,0,sigma^2)
  left_torus <- array(,dim = c(n/2,2))
  left_torus[,1] <- ((R + r)/2 + n1)*cos(theta1)-0.55
  left_torus[,2] <- ((R + r)/2 + n1)*sin(theta1)
  right_torus <- array(,dim = c(n/2,2))
  right_torus[,1] <- ((R + r)/2 + n2)*cos(theta2)+0.55
  right_torus[,2] <- ((R + r)/2 + n2)*sin(theta2)
  torus_data <- rbind(left_torus,right_torus)
  return(torus_data)
}

eight_points1 <- eight(n,r,R)
#Visualization of the eight-shape object
plot(eight_points1, xlab = 'x', ylab = 'y')
eight_points2 <- eight(n,r,R)

#Implement HPD sets estimation
train.sim <- eight_points1
test.sim <- eight_points2
conv_training <- convert.marginal(train.sim,train.sim)
conv_testing <- convert.marginal(train.sim,test.sim)
vector_tau = seq(0.5,0.01,len=10)
alpha = 0.05
out=cred.est2(test=conv_testing,train.sim=train.sim,train=conv_training,tau=vector_tau,alpha=alpha,p.alpha=0.95)

k <- out$n.tree
HPD <- out$cre.set
all_cells <- out$uni_space
complement <- all_cells[(length(HPD)+1):k]

##Visualization of the estimated HPD sets
#Plot HPD sets in red
# dev.new()
xx <- c(1,2)
f1=range(eight_points1[,xx[1]]); f2=range(eight_points1[,xx[2]])
for (i in 1:length(HPD)){
  pt=all_cells[[i]][,xx]
  pt[,1]=pt[,1]*diff(f1)+f1[1]
  pt[,2]=pt[,2]*diff(f2)+f2[1]
  xleft <- pt[1,1]; ybottom <- pt[1,2]; xright <- pt[2,1]; ytop <- pt[2,2]
  rect(xleft, ybottom, xright, ytop, density = 11, angle = 45, col = 'red')
}
#Plot complement sets in blue
for (i in (length(HPD)+1):k){
  pt=all_cells[[i]][,xx]
  pt[,1]=pt[,1]*diff(f1)+f1[1]
  pt[,2]=pt[,2]*diff(f2)+f2[1]
  xleft <- pt[1,1]; ybottom <- pt[1,2]; xright <- pt[2,1]; ytop <- pt[2,2]
  rect(xleft, ybottom, xright, ytop, density = 11, angle = 45, col = 'blue')
}

##Apply erosion method on the multivariate normal data
e <- erosion2(HPD, complement)
graph <- e[[1]]
HPD_vertices <- e[[2]][[1]]
plot(as.undirected(induced_subgraph(graph,HPD_vertices)), vertex.size = 1)
rglplot(induced_subgraph(graph, HPD_vertices), vertex.size = 1)
#plot the HPD set remains
for (i in HPD_vertices){
  pt=all_cells[[i]][,xx]
  pt[,1]=pt[,1]*diff(f1)+f1[1]
  pt[,2]=pt[,2]*diff(f2)+f2[1]
  xleft <- pt[1,1]; ybottom <- pt[1,2]; xright <- pt[2,1]; ytop <- pt[2,2]
  rect(xleft, ybottom, xright, ytop, density = 11, angle = 45, col = 'green')
}






