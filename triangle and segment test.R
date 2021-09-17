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


plucker<-function(a,b){
  p0<-a[1]*b[2]-a[2]*b[1]
  p1<-a[1]*b[3]-a[3]*b[1]
  p2<-a[1]-b[1]
  p3<-a[2]*b[3]-a[3]*b[2]
  p4<-a[3]-b[3]
  p5<-b[2]-a[2]
  pc<-c(p0,p1,p2,p3,p4,p5)
  return(pc)
}#pc for plucker coordinate
side<-function(a,b){
  s<-a[1]*b[5] + a[2]*b[6] + a[3]*b[4] + a[4]*b[3] + a[5]*b[1] + a[6]*b[2]
}#calculate the permuted inner product of 2 plucker coordinates of 2 lines

operator1<-function(x1,x2,x3,y1,y2){
  L<-plucker(y1,y2)
  e1<-plucker(x1,x2)
  e2<-plucker(x2,x3)
  e3<-plucker(x3,x1)
  s1<-side(L,e1)
  s2<-side(L,e2)
  s3<-side(L,e3)
  if (s1<0 && s2<0 && s3<0){print("proper intersection")} else if (s1>0&&s2>0&&s3>0){print("proper intersection")} else{print("coplaner or did not intersect")}
}
x1<-c(0,0,0)
x2<-c(0,5,0)
x3<-c(5,5,0)
y1<-c(11,10,-1)
y2<-c(11,10,1)
start_time <- Sys.time()
operator1(x1,x2,x3,y1,y2)
end_time <- Sys.time()
end_time - start_time
