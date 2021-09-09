##simulate a 3d-solid torus in p dimension by p*p orthogonal matrix transformation
library(alphashape3d)
library(pracma)
n <- 10000
p <- 10
upper_bound <- 1
lower_bound <- 0.5
origin <- c(0,0,0)

points_torus <- rtorus(n, lower_bound, upper_bound, ct= origin)
simulate_before <- cbind(array(points_torus,dim = c(n,3)),array(0,dim = c(n,p-3)))

rand_ortho_matrix <- randortho(10, type = 'orthonormal')
simulate_after <- simulate_before %*% rand_ortho_matrix

simulate_after
