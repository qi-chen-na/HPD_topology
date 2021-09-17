###simulation
library(alphashape3d)
library(pracma)
##simulate shells in p-dimension uniformly

shell_simulate <- function(n,p,lower_bound = 0.6, upper_bound = 1){
  i <- 0
  simulate <- array(dim = c(n, p))
  while (i < n){
    point <- runif(p,-1,1)
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
    point <- runif(p,-1,1)
    if (sqrt(sum(point ** 2)) <= 1){
      simulate[i + 1,] <- point
      i <- i + 1
    }
  }
  return (simulate)
}

##Simulate points uniformly distributed in torus in 3d and then extend to p-dimension with orthonormal transformation
torus_simulate <- function(n, p, R = 1, r = 0.6, origin = c(0, 0, 0)){
  points_torus <- rtorus(n, r, R, ct= origin)
  simulate_before <- cbind(array(points_torus,dim = c(n,3)),array(0,dim = c(n,p-3)))
  rand_ortho_matrix <- randortho(p, type = 'orthonormal')
  simulate_after <- simulate_before %*% rand_ortho_matrix
  return (simulate_after)
}


##plot the result of erosion by choosing the first two coords. 

plot_result <- function(all_cells, hpd_vertices, complement_vertices, removed_vertices){
  for(i in removed_vertices){
    pt=all_cells[[i]][,xx]
    pt[,1]=pt[,1]*diff(f1)+f1[1]; pt[,2]=pt[,2]*diff(f2)+f2[1]
    lines(c(pt[1,1],pt[2,1],pt[2,1],pt[1,1],pt[1,1]),c(pt[1,2],pt[1,2],pt[2,2],pt[2,2],pt[1,2]),col='green') 
  }
  for(i in 1:length(out$cre.set)){
    pt=all_cells[[i]][,xx]
    pt[,1]=pt[,1]*diff(f1)+f1[1]
    pt[,2]=pt[,2]*diff(f2)+f2[1]
    lines(c(pt[1,1],pt[2,1],pt[2,1],pt[1,1],pt[1,1]),c(pt[1,2],pt[1,2],pt[2,2],pt[2,2],pt[1,2]),col='blue') 
  }
  for(i in e[[2]]){
    pt=all_cells[[i]][,xx]
    pt[,1]=pt[,1]*diff(f1)+f1[1]; pt[,2]=pt[,2]*diff(f2)+f2[1]
    lines(c(pt[1,1],pt[2,1],pt[2,1],pt[1,1],pt[1,1]),c(pt[1,2],pt[1,2],pt[2,2],pt[2,2],pt[1,2]),col='red') 
  }
}

source('funcode.R')


n <- 10000
p <- 2
train.sim <- ball_simulate(n,p)
test.sim <- ball_simulate(n,p)
conv_training <- convert.marginal(train.sim,train.sim)
conv_testing <- convert.marginal(train.sim,test.sim)
alpha=0.1; #target quantile level
vector_tau = seq(0.5,0.01,len=10)
out=cred.est2(test=conv_testing,train.sim=train.sim,train=conv_training,tau=vector_tau,alpha=alpha,p.alpha=0.95)


source('Erosion.R')

#test erosion algorithm
k <- out$n.tree
HPD <- out$cre.set
all_cells <- out$uni_space
complement <- all_cells[(length(HPD)+1):k]

xx=c(1,2); f1=range(ball_points1[,xx[1]]); f2=range(ball_points1[,xx[2]])
exam.sim <- ball_simulate(n,p)
dev.new()
plot(exam.sim[,xx],col='darkgrey',xlim=range(train.sim[,xx[1]]),ylim=range(train.sim[,xx[2]]))

start_time <- Sys.time()
e <- erosion(HPD, complement)
end_time <- Sys.time()
end_time - start_time


plot_result(all_cells, e[[1]],e[[2]],e[[3]])



e[[1]]
e[[2]]
e[[3]]
all_cells
all_graph
