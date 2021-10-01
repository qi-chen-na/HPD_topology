library(Matrix)
library(igraph)
library(pracma)
library(data.table)
########################################################################################
###Adjacency
##Test if two cells are weakly adjacent
#Inputs are two cells stored in 2*d array
adjacent_weak <- function(R1, R2){
  i0 <- -1
  i <- 0
  d <- ncol(R1)
  #check if there're two interval meet at the same dimension
  for (j in 1:d){
    if ((R1[2,j] == R2[1,j]) | (R2[2,j] == R1[1,j])){
      i0 <- j
      break
    }
  }
  if (i0 == -1) return (FALSE) 
  else{
    for (j in 1:d){
      if (j != i0 & (R2[1,j] > R1[2,j] | R1[1,j] > R2[2,j])) return (FALSE)
    }
  }
  return (TRUE)
}
##Test if two cells are strongly adjacent
#Inputs are two cells stored in 2*d array
adjacent_strong <- function(R1, R2){
  i0 <- -1
  i <- 0
  d <- ncol(R1)
  for (j in 1:d){
    if ((R1[2,j] == R2[1,j]) | (R2[2,j] == R1[1,j])){
      i0 <- j
      break
    }
  }
  if (i0 == -1) return (FALSE) 
  else{
    for (j in 1:d){
      #Strongly adjacent would require all other intervals overlap, so equality at some dim would return a false.
      if (j != i0 & (R2[1,j] >= R1[2,j] | R1[1,j] >= R2[2,j])) return (FALSE)
    }
  }
  return (TRUE)
}
########################################################################################




erosion1 <- function(HPD, complement){
  l <- length(HPD)
  all_rec <- append(HPD, complement)
  n <- length(all_rec)
  d <- ncol(all_rec[[1]])
  hpd_vertices <- c(1:l)
  complement_vertices <- c((l+1):(n+1))
  
  ##Construct a graph by adjacency relationship
  #It's a directed graph where the direction always goes from lower index to larger index
  #The type of adjacency is defined by the color of the edge
  all_graph <- make_empty_graph((n+1))
  for (i in 1:n){
    R1 <- all_rec[[i]]
    for (j in i:n){
      R2 <- all_rec[[j]]
      if (adjacent_strong(R1, R2)) all_graph <- all_graph + edges(c(i,j), color = 'red') #if strong adjacent, color red
      else{
        if (adjacent_weak(R1, R2)) all_graph <- all_graph + edges(c(i,j), color = 'blue') #if weak adjacent, color blue
      }
    }
  }
  
  #Define cells at the edge of space by "if the cell share a hyperplane of dim d-1 to with the boundary of the space"
  #Check over all the cells and store the edge indices.
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
  
  #Introduce a ghost cell labelled (n+1), which is strongly adjacent to every vertices at the edge of the space
  for (i in edge_indices){
    all_graph <- all_graph + edges(c(i, (n+1)), color = 'red') #The ghost vertex is strongly adjacent to every vertex at the edge of the space
  }
  
  boundary_hpd_vertices <- c()
  boundary_complement_vertices <- c()
  removed_vertices<- c()

  ##First erode HPD sets
  #Initiated with a set of HPD sets at the boundary
  for (i in 1:l){
    all_neighbours <- adjacent_vertices(all_graph,i, mode = 'all')[[1]]
    c_neighbours <- intersect(all_neighbours, complement_vertices)
    if (isempty(c_neighbours)) next
    else {
      for (index in c_neighbours){
        if (E(all_graph, c(min(i,index),max(i,index)))$color == 'red'){
          boundary_hpd_vertices <- c(boundary_hpd_vertices, i)
          break
        }
      }
    }
  }
  change <- -1
  while (change != 0){
    before_hpd <- hpd_vertices
    before_boundary <- copy(boundary_hpd_vertices)
    for (i in before_boundary){
      all_neighbours <- adjacent_vertices(all_graph, i, mode = 'all')[[1]]
      compare_vertices <- union(complement_vertices, removed_vertices)
      hpd_neighbours <- setdiff(all_neighbours, compare_vertices)
      other_neighbours <- setdiff(all_neighbours, hpd_neighbours)
      if (isempty(hpd_neighbours)) next
      else if (is_connected(induced_subgraph(all_graph, hpd_neighbours)) & is_connected(induced_subgraph(all_graph, other_neighbours))){
        removed_vertices <- c(removed_vertices,i)
        hpd_vertices <- setdiff(hpd_vertices,i)
        boundary_hpd_vertices <- setdiff(boundary_hpd_vertices,i)
        for (j in hpd_neighbours){
          if (E(all_graph, c(min(i,j),max(i,j)))$color == 'red') boundary_hpd_vertices <- union(boundary_hpd_vertices, j)
          else{
            h <- intersect(adjacent_vertices(all_graph, j, mode = 'all')[[1]], union(complement_vertices, removed_vertices))
            for (k in h){
              if (E(all_graph, c(min(k,j), max(k,j)))$color == 'red') {
                boundary_hpd_vertices <- union(boundary_hpd_vertices, j)
                break
              }
            }
          }
        }
      }
      else next
    }
    change <- length(setdiff(before_hpd,hpd_vertices))
  }
  return (list(all_graph, list(hpd_vertices, complement_vertices, removed_vertices)))
}







erosion2 <- function(HPD, complement){
  l <- length(HPD)
  all_rec <- append(HPD, complement)
  n <- length(all_rec)
  d <- ncol(all_rec[[1]])
  hpd_vertices <- c(1:l)
  complement_vertices <- c((l+1):(n+1))
  
  ##Construct a graph by adjacency relationship
  #It's a directed graph where the direction always goes from lower index to larger index
  #The type of adjacency is defined by the color of the edge
  all_graph <- make_empty_graph((n+1))
  for (i in 1:n){
    R1 <- all_rec[[i]]
    for (j in i:n){
      R2 <- all_rec[[j]]
      if (adjacent_strong(R1, R2)) all_graph <- all_graph + edges(c(i,j), color = 'red') #if strong adjacent, color red
      else{
        if (adjacent_weak(R1, R2)) all_graph <- all_graph + edges(c(i,j), color = 'blue') #if weak adjacent, color blue
      }
    }
  }
  
  #Define cells at the edge of space by "if the cell share a hyperplane of dim d-1 to with the boundary of the space"
  #Check over all the cells and store the edge indices.
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
  
  #Introduce a ghost cell labelled (n+1), which is strongly adjacent to every vertices at the edge of the space
  for (i in edge_indices){
    all_graph <- all_graph + edges(c(i, (n+1)), color = 'red') #The ghost vertex is strongly adjacent to every vertex at the edge of the space
  }
  
  boundary_hpd_vertices <- c()
  boundary_complement_vertices <- c()
  removed_vertices<- c()
  
  ##First erode HPD sets
  #Initiated with a set of HPD sets at the boundary
  for (i in 1:l){
    all_neighbours <- adjacent_vertices(all_graph,i, mode = 'all')[[1]]
    c_neighbours <- intersect(all_neighbours, complement_vertices)
    if (isempty(c_neighbours)) next
    else {
      for (index in c_neighbours){
        if (E(all_graph, c(min(i,index),max(i,index)))$color == 'red'){
          boundary_hpd_vertices <- c(boundary_hpd_vertices, i)
          break
        }
      }
    }
  }
  change <- -1
  while (change != 0){
    before_hpd <- hpd_vertices
    before_boundary <- copy(boundary_hpd_vertices)
    for (i in before_boundary){
      all_neighbours <- adjacent_vertices(all_graph, i, mode = 'all')[[1]]
      compare_vertices <- union(complement_vertices, removed_vertices)
      hpd_neighbours <- setdiff(all_neighbours, compare_vertices)
      other_neighbours <- setdiff(all_neighbours, hpd_neighbours)
      hpd_neighbours_graph <- induced_subgraph(all_graph, hpd_neighbours)
      hpd_neighbours_graph <- delete_edges(hpd_neighbours_graph,E(hpd_neighbours_graph)[E(hpd_neighbours_graph)$color == 'blue'])
      other_neighbours_graph <- induced_subgraph(all_graph, other_neighbours)
      other_neighbours_graph <- delete_edges(other_neighbours_graph,E(other_neighbours_graph)[E(other_neighbours_graph)$color == 'blue'])
      if (isempty(hpd_neighbours)) next
      else if (is_connected(hpd_neighbours_graph) & is_connected(other_neighbours_graph)){
        removed_vertices <- c(removed_vertices,i)
        hpd_vertices <- setdiff(hpd_vertices,i)
        boundary_hpd_vertices <- setdiff(boundary_hpd_vertices,i)
        for (j in hpd_neighbours){
          if (E(all_graph, c(min(i,j),max(i,j)))$color == 'red') boundary_hpd_vertices <- union(boundary_hpd_vertices, j)
          else{
            h <- intersect(adjacent_vertices(all_graph, j, mode = 'all')[[1]], union(complement_vertices, removed_vertices))
            for (k in h){
              if (E(all_graph, c(min(k,j), max(k,j)))$color == 'red') {
                boundary_hpd_vertices <- union(boundary_hpd_vertices, j)
                break
              }
            }
          }
        }
      }
      else next
    }
    change <- length(setdiff(before_hpd,hpd_vertices))
  }
  return (list(all_graph, list(hpd_vertices, complement_vertices, removed_vertices)))
}


