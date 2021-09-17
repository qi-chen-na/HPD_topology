library(Matrix)
library(igraph)
library(data.table)
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
    if ((R1[2,j] == R2[1,j]) | (R2[2,j] == R1[1,j])){
      i0 <- j
      break
    }
  }
  if (i0 == -1) return (list(FALSE, -1, NA)) 
  else{
    coord <- rep(NA, d)
    for (j in 1:d){
      if (j != i0 & (R2[1,j] > R1[2,j] | R1[1,j] > R2[2,j])) return (list(FALSE, -1, NA)) 
      else if (j != i0 & (R2[1,j] <= R1[2,j])) coord[j] <- (R2[1,j] + R1[2,j])/2
      else if (j != i0 & (R1[1,j] <= R2[2,j])) coord[j] <- (R1[1,j] + R2[2,j])/2
      else if (j == i0 & R1[2,j] == R2[1,j]) coord[j] <- R1[2,j]
      else coord[j] <- R2[2,j]
    }
  }
  return (list(TRUE, i0, coord))
}

##Construct a undirected graph
##The input are a list of hyper-rectangles, where each is 2*d.
simple_graph <- function(all_rec){
  n <- length(all_rec)
  dim <- ncol(all_rec[[1]])
  sparse_connecting <- Matrix(FALSE, n, n, sparse = TRUE)
  for (i in 1:n){
    R1 <- all_rec[[i]]
    for (j in i:n){
      R2 <- all_rec[[j]]
      s <-adjacent_direct(R1, R2) 
      sparse_connecting[i,j] <- s[[1]]
    }
  }
  return (sparse_connecting)
}


erosion <- function(HPD, complement){
  l <- length(HPD)
  all_rec <- append(HPD, complement)
  n <- length(all_rec)
  d <- ncol(all_rec[[1]])
  hpd_vertices <- c(1:l)
  complement_vertices <- c((l+1):(n+1))
  all_matrix <- simple_graph(all_rec)
  
  #Check over all the cells and store the cells and the edge of the space
  edge_indices <- c()
  for (i in 1:length(all_cells)){
    rec <- all_cells[[i]]
    is_edge <- 0
    for (j in rec[2,]){
      if (j == 1) is_edge <- 1
    }
    for (j in rec[1,]){
      if (j == 0) is_edge <- 1
    }
    if (is_edge == 1) edge_indices <- c(edge_indices, i)
  }
  #Connecting complements at the edge to an additional vertex so that all complements are connected
  all_matrix <- cbind(rbind(all_matrix,rep(0,n)),rep(0,(n+1)))
  for (i in edge_indices) all_matrix[i,(n+1)] <- 1
  all_graph <- graph_from_adjacency_matrix(all_matrix, mode = 'undirected')

  boundary_hpd_vertices <- c()
  boundary_complement_vertices <- c()
  removed_vertices<- c()
  
  #Create a list of indicies to store the adjacent vertices, which is static
  adjacency <- vector(mode = 'list', length = (n+1))
  for (i in 1:(n+1)) adjacency[i] <- adjacent_vertices(all_graph,i)
  
  #First erode HPD sets
  #Initiated with a set of HPD sets at the boundary
  for (i in 1:l){
    if (isempty(intersect(adjacency[[i]],complement_vertices))) next
    else boundary_hpd_vertices <- c(boundary_hpd_vertices, i)
  }
  change <- -1
  while (change != 0){
    before_hpd <- hpd_vertices
    before_boundary <- copy(boundary_hpd_vertices)
    for (i in before_boundary){
      adjacent_vertices <- adjacency[[i]]
      compare_vertices <- union(complement_vertices, removed_vertices)
      hpd_neighbour <- setdiff(adjacent_vertices, compare_vertices)
      other_neighbour <- setdiff(adjacent_vertices, hpd_neighbour)
      if (isempty(hpd_neighbour)) next
      else if (is_connected(induced_subgraph(all_graph, hpd_neighbour)) & is_connected(induced_subgraph(all_graph, other_neighbour))){
        removed_vertices <- c(removed_vertices,i)
        hpd_vertices <- setdiff(hpd_vertices,i)
        boundary_hpd_vertices <- setdiff(union(boundary_hpd_vertices, hpd_neighbour),i)
      }
      else next
    }
    change <- length(setdiff(before_hpd,hpd_vertices))
    change
    length(removed_vertices)
    length(hpd_vertices)
  }
  
  #Then erode complement sets
  #Initiated with a set of complement sets at the boundary
  for (i in (l+1):n){
    if (isempty(intersect(adjacency[[i]],union(hpd_vertices,removed_vertices)))) next
    else boundary_complement_vertices <- c(boundary_complement_vertices, i)
  }
  change <- -1
  while (change != 0){
    before_complement <- complement_vertices
    before_boundary <- copy(boundary_complement_vertices)
    for (i in before_boundary){
      adjacent_vertices <- adjacency[[i]]
      compare_vertices <- union(hpd_vertices, removed_vertices)
      complement_neighbour <- setdiff(adjacent_vertices, compare_vertices)
      other_neighbour <- setdiff(adjacent_vertices, complement_neighbour)
      if (isempty(complement_neighbour)) next
      else if (is_connected(induced_subgraph(all_graph, complement_neighbour)) & is_connected(induced_subgraph(all_graph, other_neighbour))){
        removed_vertices <- c(removed_vertices,i)
        complement_vertices <- setdiff(complement_vertices,i)
        boundary_complement_vertices <- setdiff(union(boundary_complement_vertices, complement_neighbour),i)
      }
      else next
    }
    change <- length(setdiff(before_complement,complement_vertices))
  }
  return (list(hpd_vertices, complement_vertices, removed_vertices))
}

  
  
  
plot(induced_subgraph(all_graph,intersect(adjacency[[35]],hpd_vertices)))
                                                  
