##Torus in 3D
nstall.packages("alphashape3d")
library(alphashape3d)
n<-1000
r<-3
R<-6
d<-10
T1 <- rtorus(n, r, R)
##generate 3d graph into d dimensions.
A<- matrix(0, nrow=n, ncol = d-3)
D<-cbind(T1,A)
library(pracma)
P<-randortho(10, type ="orthonormal")
New_data<-D%*%P#generate a randomised nxd matrix using rows as observations that
#'theoretically' has the same topology as the original data in the 3d torus
