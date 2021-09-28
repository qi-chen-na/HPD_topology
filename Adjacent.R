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
#Construct a graph embedded in space.
#Input is a list of hyper-rectangles, each entry is stored in 2*d matrix.
##Two outputs;
#The first is a list of center of a hyper-rectangle;
#The second one is a matrix, (i,j) is NA if not adjacent, a coordinate in dimension d if adjacent.
space_graph <- function(HPD){
  n <- length(HPD)
  dim <- ncol(HPD[1])
  centre <- vector(mode = 'list', length = n)
  connecting <- array(dim = c(n, n, d))
  for (i in 1:n) centre[i] <- (HPD[i][1,]+HPD[i][2,])/2
  for (i in 1:n){
    for (j in i:n){
      R1 <- HPD[i]
      R2 <- HPD[j]
      s <-adjacent_direct(R1, R2) 
      if (s[1]) connecting[i, j, ] <- s[3]
      else connecting[i, j, ] <- NA
    }
  }
  return (list(centre, connecting))
}

#m_adjacency in the sense that two cells are adjacent if there intersection is larger than d-1-m.
#e.g. when m = 0, this is strong adjacent.
m_adjacent<- function(R1, R2, m){
  i <- 0
  d <- ncol(R1)
  for (j in 1:d){
    if (R2[1,j] > R1[2,j] | R1[1,j] > R2[2,j] | i > m + 1) return (FALSE)
    else if (R2[1,j] == R1[2,j] | R1[1,j] == R2[2,j]) i <- i + 1
  }
  return (TRUE)
}

