heron<-function(x,y,z){
  v1<-y-x
  v2<-z-y
  v3<-x-z
  l1<-sqrt(sum(v1^2))
  l2<-sqrt(sum(v2^2))
  l3<-sqrt(sum(v3^2))
  s<-(l1+l2+l3)/2
  area<-sqrt(s*(s-l1)*(s-l2)*(s-l3))
  return(area) 
}#function to calculate triangle area, input are the coordinates of the vertexes 

interior<-function(x,y,z,p){
  A0<-heron(x,y,z)
  A1<-heron(x,y,p)
  A2<-heron(y,z,p)
  A3<-heron(x,z,p)
  if (abs(A0-A1-A2-A3)<1*e^(-10)){print("point inside triangle")} else{print("point outside triangle")}
}#input are coordinates of the points in d dimensions
