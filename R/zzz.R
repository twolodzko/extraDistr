
.onUnload <- function (libpath) {
  library.dynam.unload("extraDistr", libpath)
}
