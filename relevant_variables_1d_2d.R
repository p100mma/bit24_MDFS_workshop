library(MDFS)
library(randomForest)

#Load preprocessed data
load("data_clean.RData")


#Perform the analysis and display the results
analyse.1D.2D <- function(disease.name, p.adjust.method='holm') {

 #!D and 2D statistical tests 
 fs1D <- MDFS(taxa, disease[,disease.name], dimensions=1, divisions=1, discretizations=1, range=0, p.adjust.method=p.adjust.method)
 fs2D <- MDFS(taxa, disease[,disease.name], dimensions=2, divisions=1, discretizations=1, range=0, p.adjust.method=p.adjust.method)

 #Number of relevant variables identified by the tests
 print(paste(length(fs1D$relevant.variables),'1D relevant taxa'))
 print(paste(length(fs2D$relevant.variables),'2D relevant taxa'))

 #P-P plots - always strongly recommended!
 plot(fs1D$p.value[order(fs1D$p.value)], col='red', ylab='p-value', main='p-value for 1D and 2D')
 points(fs2D$p.value[order(fs2D$p.value)],col='blue')
 legend(0,1,c('1D','2D'),col=c('red','blue'),pch=1,bty='n')

 #The difference between 1D and 2D - better visible in logarithmic scale
 plot(fs1D$p.value, fs2D$p.value, 
      col=1+
       (1:ncol(taxa))%in%fs1D$relevant.variables+
       2*(1:ncol(taxa))%in%fs2D$relevant.variables,
      pch=19, log='xy', xlab='p-value 1D', ylab='p-value 2D',main='1D vs 2D')

 #Which features are identified only by 1D and only by 2D test?
 common.features <- intersect(fs1D$relevant.variables,fs2D$relevant.variables)
 features.1D.only <- setdiff(fs1D$relevant.variables,fs2D$relevant.variables)
 features.2D.only <- setdiff(fs2D$relevant.variables,fs1D$relevant.variables)

 #How important are these features?
 ranks.1D.only <- rank(fs1D$p.value)[features.1D.only]
 ranks.2D.only <- rank(fs2D$p.value)[features.2D.only]

 print('Ranks of taxa relevant in 1D only:')
 print(round(ranks.1D.only[order(ranks.1D.only)]))
 
 print('Ranks of taxa relevant in 2D only:')
 print(round(ranks.2D.only[order(ranks.2D.only)]))
 
 #We may need MDFS results
 return(list(fs1D, fs2D))
} 

##################################################################################### 

#First example
FS12 <- analyse.1D.2D('mental_illness_type_unspecified')

#One strong interaction is visible. Now find out, which one
features.2D.only <- setdiff(FS12[[2]]$relevant.variables,FS12[[1]]$relevant.variables)
best.2D <- features.2D.only[order(FS12[[2]]$p.value[features.2D.only])][1]

#Use the ComputeInterestingTuples function to identify the partner
best.tuples <- ComputeInterestingTuples(taxa, disease[,'mental_illness_type_unspecified'], range=0, interesting.vars = best.2D)

best.tuple <- unlist(best.tuples[which.max(best.tuples$IG),c(2,3)])

print(FS12[[1]]$p.value[best.tuple])
print(FS12[[2]]$p.value[best.tuple])

#How the interaction look like?
plot(taxa[,best.tuple[1]], taxa[,best.tuple[2]], col=1+disease[,'mental_illness_type_unspecified'], 
     pch=19, xlab=colnames(taxa)[best.tuple[1]], ylab=colnames(taxa)[best.tuple[2]],cex.lab=0.66)

#What are the bacteria?
colnames(taxa)[best.tuple]

#Very strange result - is it real? Is it an artifact? Who knows?...

###################################################################################

#Second example
FS12 <- analyse.1D.2D('subset_healthy')

#Look at PP-plots. 
#There are many variables that contribute a weak signal in 2D.
#So, maybe more liberal multiple test correction is better.
FS12 <- analyse.1D.2D('subset_healthy',p.adjust.method='fdr')

#Which set of taxa is better for classification?
print('1D:')
print(randomForest(taxa[,FS12[[1]]$relevant.variables], as.factor(disease[,'subset_healthy'])))
print('2D:')
print(randomForest(taxa[,FS12[[2]]$relevant.variables], as.factor(disease[,'subset_healthy'])))





