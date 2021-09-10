x1<-c(sample(1:10,10,replace=TRUE))
x2<-c(sample(1:10,10,replace=TRUE))
x3<-c(sample(1:10,10,replace=TRUE))
y1<-c(sample(1:10,10,replace=TRUE))
y2<-c(sample(1:10,10,replace=TRUE))
l1<-x2-x1
l2<-x3-x1
l3<-y1-x1
l4<-y2-x1

data<- c(x1,x2,x3,y1,y2)
mat1<- matrix(data,nrow=5,ncol=10,byrow=TRUE)
test1<-mat1%*%t(mat1)

start_time <- Sys.time()
if(test.function(test1)=="in the same 3 dimensional subspace")
{r<-c(l1,l2,l3)
R<-matrix(r,nrow=3,ncol=10,byrow=TRUE)
n1<-c(0,0,0,0)
n2<-R%*%l1
n3<-R%*%l2
n4<-R%*%l3
n5<-R%*%l4}# should I orthogonalise?
end_time <- Sys.time()
end_time - start_time
test.function <- function(a) {
if(abs(det(a))<1e-10){
  print("in the same 3 dimensional subspace")
}else{print("does not intersect")}}
