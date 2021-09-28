# HPD_topology
########################################################################################
adjacent.R has the code for test adjacency.
###
adjacent_direct <- function(R1, R2)
inputs: Two cells.
outputs:1.TRUE if adjacent.
        2.Connecting dimension if the first one is TRUE, otherwise -1.
        3.Mid point of that on that connecting surface/hyper-plane.
###
###
space_graph <- function(HPD)
inputs: A list of cells
outputs:1.List of centers of cells.
        2.Matrix, (i,j) is NA if not adjacent, a coordinate of length d if adjacent.
###
###
m_adjacent <- function(R1,R2,m)
inputs: 1.Two cells.
        2.m a real number less than d where d is the dimension of the cell.
outputs:TRUE if two cells are m_adjacent.
########################################################################################
########################################################################################
funcode.R
Code from Nicholls & Kate's github.
Doing HPD estimation.
########################################################################################
########################################################################################