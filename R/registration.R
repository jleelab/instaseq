registration<-function(x, y, d.x, d.y){
  mat<-.Call("ThinPlateRegistration", x, y, d.x, d.y)
  return(mat)
}
