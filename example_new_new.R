rm(list=ls())
setwd("C:/Users/nicholls.NICHOLLS2389/Documents/collab - Kate/YigeQi")

#99,94,88, 86, 74,73
set.seed(99)
library(alphashape3d)
library(pracma)

source('funcode.R')
source('erosion_new_new.R')

for_grace <- function(train.sim, test.sim,alpha){
  conv_training <- convert.marginal(train.sim,train.sim)
  conv_testing <- convert.marginal(train.sim,test.sim)
  vector_tau = seq(0.5,0.01,len=10)
  out=cred.est2(test=conv_testing,train.sim=train.sim,train=conv_training,tau=vector_tau,alpha=alpha,p.alpha=0.95)
  k <- out$n.tree
  HPD <- out$cre.set
  all_cells <- out$uni_space
  complement <- all_cells[(length(HPD)+1):k]
  e <- erosion1(HPD, complement,all_cells)
  before <- list(c(1:length(HPD)),c((length(HPD)+1):k),c())
  adjacency_matrix <- as_adjacency_matrix(e[[1]], type = 'both')
  return (list(boxes=out$uni_space, before=before, after=e[[2]], adjacency_matrix=as.matrix(adjacency_matrix)))
}

 
#
library(MASS)

norm_mix_ring<-function(p,n.mix,s,n,torus=FALSE) {
  #mixture of p-dim normals (independent components, sigma=s)
  #with means distributed randomly over d-dim sphere
  z=matrix(2*runif(n.mix*p)-1,n.mix,p); 
  if (torus) z[,p]<-0
  mu=t(apply(z,1,function(x) {x/norm(x,"2")}))
  Sig=s*diag(1,p)
  X=matrix(NA,n,p)
  for (k in 1:n) {i=sample(1:n.mix,1,prob=rep(1,n.mix)/n.mix); X[k,]=mvrnorm(1,mu[i,],Sig)}
  #X<-rbind(X,matrix(4*runif(100*p)-2,100,p)) #experiment with adding background points
  return(list(samples=X,mu=mu,Sigma=Sig))
}

#plotting only works for p=2,3 so far only trivial projection by cutting dim>2
p=3; n.mix=1000; s=0.02; 
n=10000
sim<-norm_mix_ring(p,n.mix,s,n,torus=TRUE)
train.sim <- sim$samples 
sim<- norm_mix_ring(p,n.mix,s,n,torus=TRUE)
test.sim<-sim$samples

if (p==3) {
  plot3d(sim$samples,col=2,pch=".",alpha=0.1)
  points3d(sim$mu,pch=16,size=10,alpha=0.2)
} else {
  plot(sim$samples[,1:2],col=2,pch=".")
  points(sim$mu[,1:2],pch=16)
}

output<-for_grace(train.sim, test.sim, alpha=0.25)
attach(output)
box_centres<-t(sapply(boxes,function(bx) {apply(bx,2,mean)}))
n.box<-length(boxes); n.box
after.cut<-after;
after.cut[[2]]<-setdiff(after[[2]],n.box+1)

SHOW3D=TRUE
if (p==3 & SHOW3D==TRUE) {

  open3d()

  plot3d(box_centres[before[[1]],],col=1,size=5,alpha=0.5)
  points3d(box_centres[before[[2]],],col=2,size=5)
   
  for (i in before[[1]]) { for (j in before[[1]]) {
    if (adjacency_matrix[i,j]==1) lines3d(box_centres[c(i,j),1],box_centres[c(i,j),2],box_centres[c(i,j),3],alpha=0.3)
  }}
  for (i in before[[2]]) { for (j in before[[2]]) {
    if (adjacency_matrix[i,j]==1) lines3d(box_centres[c(i,j),1],box_centres[c(i,j),2],box_centres[c(i,j),3],col=2,alpha=0.3)
  }}

  open3d()
 
  plot3d(box_centres[after.cut[[1]],],col=1,size=5,alpha=0.5)
  points3d(box_centres[after.cut[[2]],],col=2,size=5)
  points3d(box_centres[after.cut[[3]],],col=3,size=5)

  for (i in after.cut[[1]]) { for (j in after.cut[[1]]) {
    if (adjacency_matrix[i,j]==1) lines3d(box_centres[c(i,j),1],box_centres[c(i,j),2],box_centres[c(i,j),3])
  }}
  for (i in after.cut[[2]]) { for (j in after.cut[[2]]) {
    if (adjacency_matrix[i,j]==1) lines3d(box_centres[c(i,j),1],box_centres[c(i,j),2],box_centres[c(i,j),3],col=2)
  }}
}

 #2d plot just chop off extra dims
  box_centres<-box_centres[,c(1,2)]
  windows(10,6)
  par(mfrow=c(1,2))
  plot(box_centres[before[[2]],],pch=16,col=2)
  points(box_centres[before[[1]],],pch=16,col=1)

  for (i in before[[1]]) { for (j in before[[1]]) {
    if (adjacency_matrix[i,j]==1) lines(box_centres[c(i,j),1],box_centres[c(i,j),2])
  }}
  for (i in before[[2]]) { for (j in before[[2]]) {
    if (adjacency_matrix[i,j]==1) lines(box_centres[c(i,j),1],box_centres[c(i,j),2],col=2)
  }}

  plot(box_centres[after.cut[[2]],],col=2,pch=16)
  points(box_centres[after.cut[[3]],],pch=16,col=3)
  points(box_centres[cbind(after.cut[[1]],after.cut[[1]]),],pch=16,col=1) #plot all twice in case there is just one so "points" works

  for (i in after.cut[[1]]) { for (j in after.cut[[1]]) {
    if (adjacency_matrix[i,j]==1) lines(box_centres[c(i,j),1],box_centres[c(i,j),2])
  }}
  #for (i in after.cut[[2]]) { for (j in after.cut[[2]]) {
  #  if (adjacency_matrix[i,j]==1) lines(box_centres[c(i,j),1],box_centres[c(i,j),2],col=2)
  #}}


detach(output)
