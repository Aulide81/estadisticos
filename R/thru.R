thru<- function(x,inicial,final){
  (inicial<=x & x<=final) & !is.na(x)
}
`%thru%`<- function(x,y) {
  thru(x,y[1],y[2])
}
