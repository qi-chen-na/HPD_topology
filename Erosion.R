##This function tells you if two hyper-rectangles are adjacent to each other.
##The inputs are two hyper-rectangles, each stored in a 2*d matrix where d is the dimension of dataset.
##Three outputs;
#The first output is logical, TRUE if adjacent;
#The second is the connecting dimension if the first one is TRUE, otherwise -1;
#The third one is the mid point of that on that connecting surface/hyper-plane.
adjacent_direct <- function(R1, R2){
  i0 <- -1
  i <- 0
  d <- ncol(R1)
  for (j in 1:d){
    if (R1[1,j] = R2){
      i0 <- j
      break
    }
  }
  if (i0 = -1) return (list(FALSE, -1, NA)) 
  else{
    coord <- rep(NA, d)
    for (j in 1:d){
      if (j != i0 & (R2[1,j] > R1[2,j] | R1[1,j] > R2[2,j])) return (list(FALSE, -1, NA)) 
      else if (j != i0 & (R2[1,j] <= R1[2,j])) coord[j] <- (R2[1,j] + R1[2,j])/2
      else if (j != i0 & (R1[1,j] <= R2[2,j])) coord[j] <- (R1[1,j] + R2[2,j])/2
      else if (j = i0 & R1[2,j] == R2[1,j]) coord[j] <- R1[2,j]
      else coord[j] <- R2[2,j]
    }
  }
  return (list(TRUE, i0, coord))
}

##Construct a undirected graph
##The input are a list of hyper-rectangles, where each is 2*d.
simple_graph <- function(all_rec){
  n <- length(all_rec)
  dim <- ncol(all_rec[1])
  sparse_connecting <- Matrix(FALSE, n, n, sparse = TRUE)
  for (i in 1:n){
    R1 <- all_rec[i]
    for (j in i:n){
      R2 <- all_rec[j]
      s <-adjacent_direct(R1, R2) 
      connecting[i,j] <- s[1]
    }
  }
  return (connecting)
}


is.removable <- function(all_graph,hpd_index,complement_index,adjacency,v,removed_index, targed = 'HPD'){
  adjacent_vertices <- adjacency[v]
  if (target = 'HPD') compare_index <- complement_index
  else compare <- hpd_index
  k <- setdiff(adjacent_vertices, union(removed_index, compare_index))
  if (length(k)=0) return (FALSE)
      #This means it's not at the boundary
  else{
      sub <- subgraph(all_graph, k)
      if (is.chordal(sub)) return (TRUE)
      else return (FALSE)
  }
}


##Erode the hyper_rectangles iteratively.

erosion <- function(hpd, complement){
  l <- length(hpd)
  n <- length(all)
  all_rec <- append(hpd, complement)
  all_matrix <- simple_graph(all)
  all_graph <- graph_from_adjacency_matrix(all_matrix, mode = 'undirected')
  hpd_index <- c(1:l)
  complement_index <- c(l+1:n)
  boundary_hpd_index <- c()
  boundary_complement_index <- c()
  removed_index<- c()
  
  #Create a list of index to store the adjacent vertices, whichis static
  adjacency <- vector(mode = 'list', length = n)
  for (i in 1:n) adjacency[i] <- adjacent_vertices(all_graph,i)
  
  #First erode HPD hyper_rectangles
  for (i in 1:l){
    if (is.removable(all_graph, hpd_index, complement_index,adjacency,i,removed_index, target = 'HPD')){
      removed_index <- c(removed_index,i)
      hpd_index <- setdiff(hpd_index,i)
      boundary_hpd_index <- setdiff(c(boundary_hpd_index, intersect(adjacency[i],hpd_index)),i)
    }
  }
  change <- 1
  while (change != 0){
    l1 <- length(hpd_index)
    candidate <- copy(boundary_hpd_index)
    for (i in candidate){
      if (is.removable(all_graph, hpd_index, complement_index,adjacency,i,removed_index, target = 'HPD')){
        removed_index <- c(removed_index,i)
        hpd_index <- setdiff(hpd_index,i)
        boundary_hpd_index <- setdiff(c(boundary_hpd_index, intersect(adjacency[i],hpd_index)),i)
      }
    l2 <- length(hpd_index)
    change <- l1-l2
    }
  }
  
  #Then erode Non-HPD hyper-rectangles
  for (i in l+1:n){
    if (is.removable(all_graph, hpd_index, complement_index,adjacency,i,removed_index, target = 'other')){
      removed_index <- c(removed_index,i)
      complement_index <- setdiff(complement,i)
      boundary_complement_index <- setdiff(c(boundary_complement_index, intersect(adjacency[i],complement_index)),i)
    }
  }
  change <- 1
  while (change != 0){
    l1 <- length(complement_index)
    candidate <- copy(boundary_complement_index)
    for (i in candidate){
      if (is.removable(all_graph, hpd_index, complement_index,adjacency,i,removed_index, target = 'other')){
        removed_index <- c(removed_index,i)
        complement_index <- setdiff(complement_index,i)
        boundary_complement_index <- setdiff(c(boundary_complement_index, intersect(adjacency[i],complement_index)),i)
      }
      l2 <- length(complement_index)
      change <- l1-l2
    }
  }
  return(hpd_index, complement_index, removed_index)
}


