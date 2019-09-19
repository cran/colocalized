#' Create a colocalization matrix
#'
#' @param bed1 A .bed table
#' @param bed2 A .bed table
#' @param n The cut-off distance
#' @return A colocalization matrix whose dimensions are the number of rows in
#'  \code{bed1} by the number of rows in \code{bed2}. Entry i,j is a 1 if the ith
#'  sequence in \code{bed1} is within the cut-off distance of the jth sequence
#'  in \code{bed2}, and 0 otherwise.
#' @examples
#' chr<-"chrY"
#' onedim_dist(nfkb1[which(nfkb1[,1]==chr),],nfkb2[which(nfkb2[,1]==chr),],150)
#' @export

onedim_dist<-function(bed1,bed2,n){
  otz<-matrix(nrow=length(bed1[,1]),ncol=length(bed2[,1]))
  for (i in 1:length(bed1[,1])){
    for (j in 1:length(bed2[,1])){
      otz[i,j]=as.numeric(min(abs(bed1[i,2]-bed2[j,3]),abs((bed1[i,3]-bed2[j,2])))<n)
    }
  }
  if (length(bed1[,1])!=length(bed2[,1])){
    otz
  }

  else{
    if (bed1[1,7]==bed2[1,7]){
      for (i in 1:length(otz[,1])){
        otz[i,i]=0
      }
      otz
    }
    otz

  }
}
