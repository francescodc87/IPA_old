### selecting things not to remove in each iteration according to RT
"checking.RT" <- function(RT,i, RT.win){
  out <- c(i,which(abs(RT-RT[i]) > RT.win))
  out
}