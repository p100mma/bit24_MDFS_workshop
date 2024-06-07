load('data_clean.RData')
library(MDFS)

disease.name<-'mental_illness_type_unspecified'

#We investigate interactions, so we are interested in 2D relevant variables
fs2D <- MDFS(taxa, disease[,disease.name], dimensions=2, divisions=1, discretizations=1, range=0, p.adjust.method='fdr')
taxa.2D<-taxa[,fs2D$relevant.variables]

#Use ComputeInterestingTuples function to compute mutual information between taxa
similarities<-pairwiseMIsimilarity(taxa.2D,range=0)$S
rownames(similarities)<-colnames(similarities)<-colnames(taxa.2D)

#Do the similarities depend on the decision variable?
similarities0<-pairwiseMIsimilarity(taxa.2D[disease[,disease.name]==0,],range=0)$S
rownames(similarities0)<-colnames(similarities0)<-colnames(taxa.2D)

similarities1<-pairwiseMIsimilarity(taxa.2D[disease[,disease.name]==1,],range=0)$S
rownames(similarities1)<-colnames(similarities1)<-colnames(taxa.2D)

#Show them graphically
heatmap(similarities,symm=T)
heatmap(similarities0,symm=T)
heatmap(similarities1,symm=T)

#rescaling helps to highlight differences
heatmap(sqrt(similarities),symm=T)
heatmap(sqrt(similarities0),symm=T)
heatmap(sqrt(similarities1),symm=T)

#find interesting pairs: for each decision class network,
#                        calculate differences between plain connections
#                        then find biggest deviations

deviation_MIs<- function(sim, sim0, sim1) {
  ltr<- lower.tri(sim)
  s0d<- -(sim-sim0)
  s1d<- -(sim-sim1)
  pmax(s0d,s1d)-> deviations
   data.frame(deviation=as.vector(deviations[ltr]), pair1=as.vector(col(deviations)[ltr]), 
                                               pair2=as.vector(row(deviations)[ltr]))-> devdf
   devdf[order(-devdf$deviation),]->devdf
   devdf$pnr<-1:nrow(devdf)
   devdf
}
##sort pairs of taxa according to difference between in class MI and overall MI
##now we can browse pairs one by one and look for interesting interactions
deviation_MIs(similarities,similarities0,similarities1)-> pair_df

#auxillary function for interactive exploration of results
#use output of previous function to plot pair nr "pnr" from the list
plot_pair<-function(taxa_table,pair_table,pnr,disease.vec){
  par(mfrow=c(1,3))
  p1<-pair_table[pair_table$pnr==pnr, 2]
  p2<-pair_table[pair_table$pnr==pnr, 3]
  xl= c( min(taxa_table[,p1]),max(taxa_table[,p1]))
  yl= c( min(taxa_table[,p2]),max(taxa_table[,p2]))
  plot(taxa_table[,p1], taxa_table[,p2],col=disease.vec + 14, pch=16, xlim=xl,ylim=yl)
  plot(taxa_table[disease.vec==0,p1], taxa_table[disease.vec==0,p2],col=disease.vec[disease.vec==0] + 14,
       xlim=xl,ylim=yl,pch=16)
  plot(taxa_table[disease.vec==1,p1], taxa_table[disease.vec==1,p2],col=disease.vec[disease.vec==1] + 14,
       xlim=xl,ylim=yl, pch=16)
  par(mfrow=c(1,1))
  
}
#other diseases are also worth checking out
#be sure to re calculate taxa.2D also!
plot_pair(taxa.2D,pair_df,1,disease[,disease.name])



#Alternatively, one can use that approach on the whole set of taxa
#In an unsupervised manner, without consideration of any decision variable

pairwiseMIsimilarityDiscrete(data=taxa.binary, divisions=1)-> mi_data


heatmap(sqrt(mi_data$S), symm=TRUE, main="all taxa")
heatmap(log(mi_data$S + 1e-8), symm=TRUE, main="all taxa")

#look from network analysis perspective:
#graphs based on MI might not be very sparse when 
#one does not account for weight of the connections
## number of connections / # possible connections:
Con_Density<-sum(mi_data$S[lower.tri(mi_data$S)]>0)/( ncol(mi_data$S)*(ncol(mi_data$S)-1)/2  )
## average degree vs number of nodes in total
mean(colSums(mi_data$S>0))
ncol(mi_data$S)
## on average, each node is connected to more than half of all other nodes
ncol(mi_data$S)/4
##degree distributions
par(mfrow=c(1,2))
hist(colSums(mi_data$S), main="weighted degree distr.") 
hist(colSums(mi_data$S>0), main="UNweighted degree distr.") 
par(mfrow=c(1,1))
# comparison of number of connections per taxa vs total of their weights
## we have a correlation, although rather non linear
## meaning there are still nodes with considerable number of connections 
## which are not that strong
## so weights matter!
plot(colSums(mi_data$S>0), colSums(mi_data$S), pch=16,
     main="degree vs weighted degree",
     xlab="# of connections", ylab="total strength")

#example graph based analysis based on MI network
source('basic_cliques.R')

clique_partitioning(mi_data$S)-> clique_list
clusterList2memVec(clique_list)-> clq_membership

#visualize clique decomposition of MI similarity graph by heatmap:
##1. add taxa names
S<- mi_data$S
colnames(S)<-rownames(S)<-colnames(taxa.binary)
##2. reorder according to clique membership
S[ order(clq_membership), order(clq_membership) ]-> S_reordered
clq_membership[ order(clq_membership) ]-> mem_reordered
##3. create color labels for cliques
colors()[ mem_reordered + 17 ] -> colorVec
##first three -- color signifies weight
heatmap( S_reordered, symm=TRUE, Rowv=NA, Colv=NA, ColSideColors= colorVec, RowSideColors=colorVec)
heatmap( sqrt(S_reordered), symm=TRUE, Rowv=NA, Colv=NA, ColSideColors= colorVec, RowSideColors=colorVec)
heatmap( log(S_reordered+1e-8), symm=TRUE, Rowv=NA, Colv=NA, ColSideColors= colorVec, RowSideColors=colorVec)
##here - light/dark -> unconnected/connected (unsignif. MI/signif. MI)
heatmap( (S_reordered>0)*1., symm=TRUE, Rowv=NA, Colv=NA, ColSideColors= colorVec, RowSideColors=colorVec)
