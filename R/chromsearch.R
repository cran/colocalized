
#' Search one chromsome
#'
#' Search a single chromosome for clusters of TF binding sequences.
#' Example produces a null result, test the same complex on "chr9"
#' for a positive reading.
#'
#' @param choose List of .bed tables
#' @param chrom Chromosome to be searched given as e.g. "chr19"
#' @param n Cut-off distance between colocalized sequences
#' @return A table containing the addresses (as one dimensional intervals) of the members of every cluster, with some annotation data.
#' @examples
#' complex<-list(nfkb1,nfkb2,relb)
#' chromsearch(complex,150,"chrY")
#' @export


chromsearch<-function(choose,n,chrom){
  tin<-choose
  ChromOutput<-colocalized(tin,chrom,n)
  raw_matrix<-Reduce(cbind,ChromOutput)[c(1,2,3,4,7,10,11,12,15,18,19,20,23),]
  if (is.null(raw_matrix)){
    OutputTable<-NULL
    OutputTable
  }

  else {
    raw_matrix<-sapply(1:length(raw_matrix[1,]), function(i) unlist(raw_matrix[,i]))
    UniqueColumns<-function(mat){
      i<-1
      RecursionBody<-function(mat,i)
        if (i>=length(mat[1,])){
          mat
        }
      else {
        find_duplicate<-sapply(1:length(mat[1,]), function(j) length(which(mat[,i] %in% mat[,j])))
        find_duplicate[i]<-0
        RecursionBody(mat[,-which(find_duplicate==(length(mat[,i])))],i+1)
      }
      RecursionBody(mat,i)
    }

    #OutputTable<-t(UniqueColumns(raw_matrix))
    #colnames(OutputTable)<-c("Chromosome",rep(c("IntervalStart","IntervalEnd","Sequence","Name"),(length(OutputTable[1,])-1)/4))

    leave<-raw_matrix[,which(sapply(1:length(raw_matrix[1,]),function(i) !(is.na(raw_matrix[2,i]))))]

    flip<-raw_matrix[,which(sapply(1:length(raw_matrix[1,]),function(i) is.na(raw_matrix[2,i])))]

    flip<-flip[1,which(!is.na(flip[1,]))]

    if (length(flip)>0){
      flipped<-cbind(leave,sapply(0:(length(flip)/24-1), function(i) flip[c(1,2,3,4,7,10,11,12,15,18,19,20,23)+24*i]))

      finaltrim<-sapply(1:length(unique(flipped[1,])), function(i) UniqueColumns(flipped[,which(flipped[1,]==unique(flipped[1,])[i])]))
      OutputTable<-t(Reduce(cbind,finaltrim))

      OutputTable

    }

    else{
    t(UniqueColumns(raw_matrix))
    }
  }
}


