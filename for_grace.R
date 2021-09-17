source('funcode.R')
source('erosion_new.R')

for_grace <- function(train.sim, test.sim){
  conv_training <- convert.marginal(train.sim,train.sim)
  conv_testing <- convert.marginal(train.sim,test.sim)
  alpha=0.1; #target quantile level
  vector_tau = seq(0.5,0.01,len=10)
  out=cred.est2(test=conv_testing,train.sim=train.sim,train=conv_training,tau=vector_tau,alpha=alpha,p.alpha=0.95)
  k <- out$n.tree
  HPD <- out$cre.set
  all_cells <- out$uni_space
  complement <- all_cells[(length(HPD)+1):k]
  e <- erosion_new(HPD, complement)
  before <- list(c(1:length(HPD)),c((length(HPD)+1):k),c())
  adjacency_matrix <- as_adjacency_matrix(e[[1]], type = 'both')
  return (list(out$uni_space, before, e[[2]], adjacency_matrix))
}