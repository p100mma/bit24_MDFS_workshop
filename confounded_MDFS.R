library(MDFS)

source('confounders.R')

#Load preprocessed data
load("data_clean.RData")

#Take the age into account as a confounder in MDFS 1D
#As previously, encapsulate the analysis into a function
age.influence <- function(disease.name) {
  decision<-disease[,disease.name]
 #Find the best discretization for the age into 2 categories
 discr.age <- opt.discretize(patient.data[,'age.years'],decision,divisions=1) 
 age.discr <- discr.age$x-1

 #You can see the optimal breakpoints
 print(paste('Breaks at:',discr.age$breaks.opt,'years'))

 #The decision variable
 

 #MDFS without confounding variables
 fs.raw <- MDFS(taxa, decision, dimensions=1, divisions=1, discretizations=1, range=0)

 #MDFS conditioned on the age
 fs.age <- MDFS.confounded(taxa, decision, age.discr, dimensions=1, divisions=1, discretizations=1, range=0, n.contrast=0)

 print(paste(length(fs.raw$relevant.variables),'relevant taxa by raw test'))
 print(paste(length(fs.age$relevant.variables),'relevant taxa by age-adjusted test'))
 
 #P-P plots - always strongly recommended!
 plot(fs.raw$p.value[order(fs.raw$p.value)], col='red', ylab='p-value', main='p-value with and without age adjustment')
 points(fs.age$p.value[order(fs.age$p.value)], col='blue')
 legend(0, 1, c('raw','age-adjusted'), col=c('red', 'blue'), pch=1, bty='n')
 
 #The difference between p-values - better visible in logarithmic scale
 plot(fs.raw$p.value, fs.age$p.value,
      col=1+
       (1:ncol(taxa))%in%fs.raw$relevant.variables+
       2*(1:ncol(taxa))%in%fs.age$relevant.variables,
      pch=19, log='xy', xlab='p-value raw', ylab='p-value adjusted',main='Raw vs age-adjusted p-value')
 
 #Which features are identified only by raw and only by age-adjusted test?
 common.features <- intersect(fs.raw$relevant.variables,fs.age$relevant.variables)
 features.raw.only <- setdiff(fs.raw$relevant.variables,fs.age$relevant.variables)
 features.age.only <- setdiff(fs.age$relevant.variables,fs.raw$relevant.variables)
 
 #How important are these features?
 ranks.raw.only <- rank(fs.raw$p.value)[features.raw.only]
 ranks.age.only <- rank(fs.age$p.value)[features.age.only]
 
 print('Ranks of taxa relevant in raw test only:')
 print(round(ranks.raw.only[order(ranks.raw.only)]))
 
 print('Ranks of taxa relevant in adjusted test only:')
 print(round(ranks.age.only[order(ranks.age.only)]))
 
 #We may need MDFS results
 return(list(fs.raw, fs.age))
} 

#Test our two examples
FS.mental <- age.influence('mental_illness_type_unspecified')
 
FS.healthy <- age.influence('subset_healthy')



