function (i, j, k, h, M) {
  findpos <- dget("./multicorr/findpos.R")
  temp <- findpos(i, j, k, h)
  fpho <- M[temp]
  return(fpho)
}
