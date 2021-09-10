##Compute barycentric coordinates (u, v, w) for
## point p with respect to triangle (a, b, c)
f<-function(a, b, c, p){
  v0 = b-a;
  v1 = c-a;
  v2 = p-a;
  d00 = Dot(v0, v0);
  d01 = Dot(v0, v1);
  d11 = Dot(v1, v1);
  d20 = Dot(v2, v0);
  d21 = Dot(v2, vl);
  denom = d00 * dll-d01 * d01;
  v = (dll * d20 - d01 * d21) / denom;
  w = (d00 * d21 - d01 * d20) / denom;
  u=1.0-v-w;
}
if(0<v&w&u<1){
  print("TRUE")
}else{
  print("FALSE")
}
