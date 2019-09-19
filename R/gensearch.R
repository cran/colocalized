#'@import foreach
#'@import doParallel
NULL

#'NFKB1 bed file
#'
#'@format A dataframe with 230505 rows and 8 columns
#'@source \url{https://ccg.epfl.ch/pwmtools/pwmscan.php}
"nfkb1"


#'NFKB2 bed file
#'
#'@format A dataframe with 1901 rows and 8 columns
#'@source \url{https://ccg.epfl.ch/pwmtools/pwmscan.php}
"nfkb2"

#'RELB bed file
#'
#'@format A dataframe with 1448 rows and 8 columns
#'@source \url{https://ccg.epfl.ch/pwmtools/pwmscan.php}
"relb"



#'Colocalized cluster search.
#'
#'Searches for clusters of colocalized transcription factor (TF) binding sequences.
#'\code{colocalized(choose,chr,n)} searches for instances where the sequences from each
#'table element in choose are colocalized to within a cut-off distance.
#'
#'@inheritParams gensearch
#'@param chr Chromosome
#'@param cores Number of cores for parallel processing. Leaving this blank causes the program to use default (series) processing
#'@return Table of clusters found in \code{chr}
#'@examples
#'complex<-list(nfkb1,nfkb2,relb)
#'colocalized(complex,"chrY",150)
#'@export

colocalized<-function(choose,chr,n,cores){

  if(missing(cores)){
    colocalized_sequential(choose,chr,n)
  }
  else{

  TF.binding.sites<-lapply(1:length(choose), function(i) choose[[i]][which(choose[[i]][,1]==chr),])


  for (i in 1:length(TF.binding.sites)){
    TF.binding.sites[[i]][,6]<-as.character(TF.binding.sites[[i]][,6])
    TF.binding.sites[[i]][,1]<-as.character(TF.binding.sites[[i]][,1])
    TF.binding.sites[[i]][,4]<-as.character(TF.binding.sites[[i]][,4])
    TF.binding.sites[[i]][,7]<-as.character(TF.binding.sites[[i]][,7])
  }

  pairs<-utils::combn(1:length(TF.binding.sites),2)





  if (length(TF.binding.sites)==2){
    colocalized.special<-onedim_dist(TF.binding.sites[[1]],TF.binding.sites[[2]],n)

    locations.special<-sapply(1:length(which(colocalized.special==1,arr.ind=TRUE)[,1]), function(i) TF.binding.sites[[1]][which(colocalized.special==1,arr.ind = TRUE)[i,1],2])

    collect.special<-sapply(1:length(locations.special), function(i) sapply(1:length(TF.binding.sites), function(j) if(as.character(TF.binding.sites[[j]][1,4])!=as.character(TF.binding.sites[[1]][1,4])|j==1) TF.binding.sites[[j]][which(abs(TF.binding.sites[[j]][,3]-locations.special[i])==min(abs(TF.binding.sites[[j]][,3]-locations.special[i]))),] else TF.binding.sites[[j]][which(abs(TF.binding.sites[[j]][,3]-locations.special[i])==sort(abs(TF.binding.sites[[j]][,3]-locations.special[i]))[j]),]))
    collect.special
  }




  else{

    doParallel::registerDoParallel(cores=cores)

    colocalized.pairs.matrix.list<-foreach::foreach(i=1:length(pairs[1,])) %dopar% onedim_dist(TF.binding.sites[[pairs[1,i]]],TF.binding.sites[[pairs[2,i]]],n)

    inv.pairs.matrix.list<-lapply(1:length(colocalized.pairs.matrix.list), function(i) t(colocalized.pairs.matrix.list[[i]]))

    pairs.matrix.list<-append(colocalized.pairs.matrix.list,inv.pairs.matrix.list)

    search.order<-utils::combn(1:length(pairs.matrix.list),(length(TF.binding.sites)-1))


    all.pairs<-cbind(pairs,sapply(1:length(pairs[1,]),function(i) rev(pairs[,i])))

    dimensionality<-sapply(1:length(search.order[1,]), function(i) prod(sapply(1:(length(all.pairs[,search.order[,i]][1,])-1),function(j) as.numeric(all.pairs[,search.order[,i]][2,j]==all.pairs[,search.order[,i]][1,(j+1)]))))



    unqs<-as.numeric(sapply(1:length(search.order[1,]), function(i) length(unique(as.vector(all.pairs[,search.order[,i]])))==length(TF.binding.sites)))

    valid.sequence<-dimensionality*unqs

    coordinates<-foreach::foreach(i=1:length(pairs.matrix.list)) %dopar% {which(pairs.matrix.list[[i]]==1, arr.ind=TRUE)}

    mor<-list()
    mig<-list()

    clustering<-function(x,n,i){
      mor<-list()
      mor[[1]]<-x[[n[1,i]]]
      for (j in 2:length(n[,i])){
        mor[[j]]<-x[[n[j,i]]][which(x[[n[j,i]]][,1] %in% mor[[j-1]][,2]),]
      }
      mor[[length(n[,1])]]
    }


    mig<-foreach::foreach(i=which(valid.sequence==1)) %dopar% clustering(coordinates,search.order,i)

    locations<-mig

    doParallel::stopImplicitCluster()
    address<-list()
    dist<-list()
    collect<-list()
    for (i in 1:length(locations)){
      if (length(dim(locations[[i]]))>1){
        address[[i]]<-TF.binding.sites[[all.pairs[,search.order[,which(valid.sequence==1)[i]][length(search.order[,which(valid.sequence==1)[1]])]][2]]][locations[[i]][,2],3]


        collect[[i]]<-list()
        for (k in 1:length(address[[i]])){
          collect[[i]][[k]]<-sapply(1:length(TF.binding.sites), function(j) unlist(TF.binding.sites[[j]][which(abs(TF.binding.sites[[j]][,3]-address[[i]][k])==min(abs(TF.binding.sites[[j]][,3]-address[[i]][k]))),]))

          same<-list()

          for (j in 1:length(unique(unlist(collect[[i]][[k]][7,])))){

            same[[j]]<-which(unlist(collect[[i]][[k]][7,])==unique(unlist(collect[[i]][[k]][7,]))[j])

            if (length(same[[j]])>1){

              for (y in 2:length(same[[j]])){
                collect[[i]][[k]][,same[[j]][y]]<-as.character(TF.binding.sites[[same[[j]][1]]][which(abs(TF.binding.sites[[same[[j]][1]]][,3]-address[[i]][k])==sort(abs(TF.binding.sites[[same[[j]][1]]][,3]-address[[i]][k]))[y]),])
              }
            }
          }
        }



        #Reduce(cbind,lapply(1:length(collect[[i]]), function(j) Reduce(rbind,collect[[i]][[j]])))[c(1,2,3,4,7,10,11,12,15,18,19,20,23),]

        #if (length(collect[[i]][1,1])==0){
        # collect[[i]]<-NULL
        #}

      }
    }

    address<-purrr::compact(address)
    collect<-purrr::compact(collect)


    collect<-lapply(1:length(collect), function(i) lapply(1:length(collect[[i]]), function(j) Reduce(rbind,collect[[i]][[j]])))
    collect<-lapply(1:length(collect), function(i) Reduce(cbind,collect[[i]]))
    collect
  }
 }
}

#'Colocalized full search.
#'
#'Wrapper for \code{colocalized} that searches every chromosome shared
#'between the given .bed files.
#'
#'@inheritParams gensearch
#'@return List of lists of each cluster found in each chromosome.
#'@export

ColocalizedFullSearch<-function(choose,n,cores){



  tan<-lapply(1:length(choose), function(i) as.character(unique(choose[[i]][,1])))

  tip<-Reduce(intersect,tan)

  if(missing(cores)){
    Full.Search.List<-lapply(1:length(tip), function(i) colocalized_sequential(choose,tip[i],n))
  }
  else{
    Full.Search.List<-lapply(1:length(tip), function(i) colocalized(choose,tip[i],n,cores))

  }
}
#' Whole genome search.
#'
#'Search the whole genome for clusters of colocalized TF binding sequences.
#' @param cores Number of cores for parallel processing. Leaving this blank causes the program to use default (series) processing
#' @param n The cut-off distance
#' @param choose List of .bed tables
#' @return Table containing the addresses (as one dimensional intervals) of the members of every cluster, with some annotation data.
#' @export





gensearch<-function(choose,n,cores){
  OutputTable<-ColocalizedFullSearch(choose,n,cores)

  size<-sapply(1:length(OutputTable), function(i) length(OutputTable[[i]]))

  column_collapse<-lapply(1:length(OutputTable[which(size>0)]), function(i) Reduce(cbind,OutputTable[which(size>0)][[i]])[c(1,2,3,4,7,10,11,12,15,18,19,20,23),])

  raw_matrix<-Reduce(cbind,column_collapse)
  #raw_matrix<-sapply(1:length(raw_matrix[1,]), function(i) unlist(raw_matrix[,i]))

  UniqueColumns<-function(mat){
    if (class(mat)=="character"){
      mat
    }

    else{
      RecursionBody<-function(mat,i){
        if (i>length(mat[1,])){
          mat
        }
        else {
          find_duplicate<-sapply(1:length(mat[1,]), function(j) length(which(mat[,i] %in% mat[,j])))
          find_duplicate[i]<-0
          RecursionBody(mat[,which(find_duplicate!=(length(mat[,i])))],i+1)
        }
      }
      RecursionBody(mat,1)
    }
  }
  leave<-raw_matrix[,which(sapply(1:length(raw_matrix[1,]),function(i) !(is.na(raw_matrix[2,i]))))]


  flip<-raw_matrix[,which(sapply(1:length(raw_matrix[1,]),function(i) is.na(raw_matrix[2,i])))]

  flip<-flip[1,which(!is.na(flip[1,]))]

  if (length(flip)>0){
  flipped<-cbind(leave,sapply(0:(length(flip)/24-1), function(i) flip[c(1,2,3,4,7,10,11,12,15,18,19,20,23)+24*i]))

  finaltrim<-sapply(1:length(unique(flipped[1,])), function(i) UniqueColumns(flipped[,which(flipped[1,]==unique(flipped[1,])[i])]))

  OutputTable<-t(Reduce(cbind,finaltrim))
  }

  else{
    OutputTable<-t(UniqueColumns(raw_matrix))
  }
  colnames(OutputTable)<-c("Chromosome",rep(c("IntervalStart","IntervalEnd","Sequence","Name"),(length(OutputTable[1,])-1)/4))
  OutputTable
}


#'Sequential cluster search
#'
#'Search one chromosome for clusters using default non-parallel processing.
#'@param choose List of .bed tables
#'@param chr Chromosome to be searched given as e.g. "chr19"
#'@param n Cut-off distance between colocalized sequences
#'@return Table containing the addresses (as one dimensional intervals) of the members of every cluster, with some annotation data.
#'@examples
#'complex<-list(nfkb1,nfkb2,relb)
#'colocalized_sequential(complex,"chrY",150)
#'@export

colocalized_sequential<-function(choose,chr,n){

  TF.binding.sites<-lapply(1:length(choose), function(i) choose[[i]][which(choose[[i]][,1]==chr),])


  for (i in 1:length(TF.binding.sites)){
    TF.binding.sites[[i]][,6]<-as.character(TF.binding.sites[[i]][,6])
    TF.binding.sites[[i]][,1]<-as.character(TF.binding.sites[[i]][,1])
    TF.binding.sites[[i]][,4]<-as.character(TF.binding.sites[[i]][,4])
    TF.binding.sites[[i]][,7]<-as.character(TF.binding.sites[[i]][,7])
  }

  pairs<-utils::combn(1:length(TF.binding.sites),2)





  if (length(TF.binding.sites)==2){
    colocalized.special<-onedim_dist(TF.binding.sites[[1]],TF.binding.sites[[2]],n)

    locations.special<-sapply(1:length(which(colocalized.special==1,arr.ind=TRUE)[,1]), function(i) TF.binding.sites[[1]][which(colocalized.special==1,arr.ind = TRUE)[i,1],2])

    collect.special<-sapply(1:length(locations.special), function(i) sapply(1:length(TF.binding.sites), function(j) if(as.character(TF.binding.sites[[j]][1,4])!=as.character(TF.binding.sites[[1]][1,4])|j==1) TF.binding.sites[[j]][which(abs(TF.binding.sites[[j]][,3]-locations.special[i])==min(abs(TF.binding.sites[[j]][,3]-locations.special[i]))),] else TF.binding.sites[[j]][which(abs(TF.binding.sites[[j]][,3]-locations.special[i])==sort(abs(TF.binding.sites[[j]][,3]-locations.special[i]))[j]),]))
    collect.special
  }




  else{

    #doParallel::registerDoParallel(cores=2)

    colocalized.pairs.matrix.list<-lapply(1:length(pairs[1,]), function(i) onedim_dist(TF.binding.sites[[pairs[1,i]]],TF.binding.sites[[pairs[2,i]]],n))

    inv.pairs.matrix.list<-lapply(1:length(colocalized.pairs.matrix.list), function(i) t(colocalized.pairs.matrix.list[[i]]))

    pairs.matrix.list<-append(colocalized.pairs.matrix.list,inv.pairs.matrix.list)

    search.order<-utils::combn(1:length(pairs.matrix.list),(length(TF.binding.sites)-1))


    all.pairs<-cbind(pairs,sapply(1:length(pairs[1,]),function(i) rev(pairs[,i])))

    dimensionality<-sapply(1:length(search.order[1,]), function(i) prod(sapply(1:(length(all.pairs[,search.order[,i]][1,])-1),function(j) as.numeric(all.pairs[,search.order[,i]][2,j]==all.pairs[,search.order[,i]][1,(j+1)]))))



    unqs<-as.numeric(sapply(1:length(search.order[1,]), function(i) length(unique(as.vector(all.pairs[,search.order[,i]])))==length(TF.binding.sites)))

    valid.sequence<-dimensionality*unqs

    #coordinates<-foreach::foreach(i=1:length(pairs.matrix.list)) %dopar% {which(pairs.matrix.list[[i]]==1, arr.ind=TRUE)}
    coordinates<-lapply(1:length(pairs.matrix.list), function(i) which(pairs.matrix.list[[i]]==1, arr.ind=TRUE))
    mor<-list()
    mig<-list()

    clustering<-function(x,n,pairs,sites,i){
      mor<-list()
      if(prod(sapply(1:length(n[,i]), function(j) dim(x[[n[j,i]]])[1]))==0){
        NULL
      }
      else{
        mor[[1]]<-x[[n[1,i]]]
        for (j in 2:length(n[,i])){

          if (as.character(sites[[pairs[,n[j,i]][2]]][1,7])!=as.character(sites[[pairs[,n[j-1,i]][1]]][1,7])){
            mor[[j]]<-x[[n[j,i]]][which(x[[n[j,i]]][,1] %in% mor[[j-1]][,2]),]
          }

          else{
            mor[[j]]<-x[[n[j,i]]][which(x[[n[j,i]]][,1] %in% mor[[j-1]][,2]),]
            mor[[j]]<-mor[[j]][-which(sapply(1:length(mor[[j]][,1]), function(z) as.logical(sum(which(mor[[j-1]][,1]==mor[[j]][z,2])==which(mor[[j-1]][,2]==mor[[j]][z,1]))))),]
          }

        }
        mor[[length(n[,1])]]
      }
    }


    mig<-lapply(which(valid.sequence==1), function(i) clustering(coordinates,search.order,all.pairs,TF.binding.sites,i))

    if (length(unlist(mig))==0){
      NULL
    }

    else{
    locations<-mig

    #doParallel::stopImplicitCluster()
    address<-list()
    dist<-list()
    collect<-list()
    for (i in 1:length(locations)){
      if (length(dim(locations[[i]]))>1){
        address[[i]]<-TF.binding.sites[[all.pairs[,search.order[,which(valid.sequence==1)[i]][length(search.order[,which(valid.sequence==1)[1]])]][2]]][locations[[i]][,2],3]


        collect[[i]]<-list()
        for (k in 1:length(address[[i]])){
          collect[[i]][[k]]<-sapply(1:length(TF.binding.sites), function(j) unlist(TF.binding.sites[[j]][which(abs(TF.binding.sites[[j]][,3]-address[[i]][k])==min(abs(TF.binding.sites[[j]][,3]-address[[i]][k]))),]))

          same<-list()

          for (j in 1:length(unique(unlist(collect[[i]][[k]][7,])))){

            same[[j]]<-which(unlist(collect[[i]][[k]][7,])==unique(unlist(collect[[i]][[k]][7,]))[j])

            if (length(same[[j]])>1){

              for (y in 2:length(same[[j]])){
                collect[[i]][[k]][,same[[j]][y]]<-as.character(TF.binding.sites[[same[[j]][1]]][which(abs(TF.binding.sites[[same[[j]][1]]][,3]-address[[i]][k])==sort(abs(TF.binding.sites[[same[[j]][1]]][,3]-address[[i]][k]))[y])[1],])
              }
            }
          }
        }



        #Reduce(cbind,lapply(1:length(collect[[i]]), function(j) Reduce(rbind,collect[[i]][[j]])))[c(1,2,3,4,7,10,11,12,15,18,19,20,23),]

        #if (length(collect[[i]][1,1])==0){
         # collect[[i]]<-NULL
        #}

      }
    }

    address<-purrr::compact(address)
    collect<-purrr::compact(collect)


    collect<-lapply(1:length(collect), function(i) lapply(1:length(collect[[i]]), function(j) Reduce(rbind,collect[[i]][[j]])))
    collect<-lapply(1:length(collect), function(i) Reduce(cbind,collect[[i]]))
    collect
    }
  }
}





