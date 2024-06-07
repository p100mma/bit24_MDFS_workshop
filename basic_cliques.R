	if (!require(matrixStats)) install.packages('matrixStats')	
	library(matrixStats)
	clique_partitioning<- function(weighted_mat) {
	cliques<-list()

	while(any(!is.na(weighted_mat))) {
	  new<-cluster<-which.max(colMaxs(weighted_mat,na.rm=TRUE))

	  while (length(new)>0) {
	    weighted_mat[new,]<-0*weighted_mat[new,]

	   potential<-which(rowAlls(weighted_mat[,cluster,drop=FALSE]>0,na.rm=TRUE))# & rowAnys(F[,cluster,drop=FALSE]>0,na.rm=TRUE)) #unnecessary because all-NA rows cannot happen
	   new.index<-which.max(rowMeans(weighted_mat[potential,cluster,drop=FALSE],na.rm=TRUE))
	   new<-potential[new.index]

	   cluster<-c(cluster,new)
	  }

	  weighted_mat[,cluster]<-NA

	  cliques[[length(cliques)+1]]<-cluster

	#  print(paste(length(cliques), length(cluster)))
	}
	return(cliques)
	}


#function to convert list of cliques to clique membership vector, to simplify comparison

	clusterList2memVec<- function(clusterList,minSize=3){
		n_obj<- length(unlist(clusterList))
		memVec<- rep(NA, n_obj)
		for (j in seq_along(clusterList))
			{
		if(length(clusterList[[j]])>=minSize) 	memVec[ clusterList[[j]] ] = j else memVec[ clusterList[[j]] ]=0
			}
		stopifnot(all(!is.na(memVec)))
		return(memVec)
	}
		
