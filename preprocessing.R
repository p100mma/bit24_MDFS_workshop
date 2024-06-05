metadata.filename <- "metadata.csv"
taxa.filename <- "taxa.tsv"

#Read microbiome data
taxa <- as.data.frame(t(read.delim(taxa.filename, header = TRUE, sep = "\t")))
colnames(taxa) <- taxa[1,]
taxa <- taxa[-1,]

#remove Archaea and Eukaryotes
archaea.eucaryotes <- c(1:9, 1022:1028)
taxa <- taxa[,-archaea.eucaryotes]

for(i in 1:ncol(taxa)) taxa[,i] <- as.numeric(taxa[,i])

#Read metadata
metadata<-read.csv(metadata.filename, header = TRUE, sep = ",")
rownames(metadata)<-metadata[,1]
metadata<-metadata[,-1]

#Disease variables 
disease <- metadata[, c(6:9, 22:29, 33, 35, 36, 48:50, 54, 67:75, 98:110, 113:115,
                        134, 202:206, 212, 216, 218:229, 234:240, 264, 266, 289:293, 298, 300)]

#Binarize disease data (1- disease, 0- no disease)
disease[disease %in% c("Self-diagnosed",
                       "Diagnosed by an alternative medicine practitioner",
                       "Diagnosed by a medical professional (doctor, physician assistant)",
                       "Unspecified",
                       "Surgery only",
                       "Radiation therapy",
                       "Chemotherapy",
                       "Quite worried",
                       "A little worried",
                       "Type II diabetes",
                       "Type I diabetes",
                       "Gestational diabetes")]<-1
disease[is.na(disease)]<-0
disease[disease!=1]<-0
for (i in 1:ncol(disease)) disease[[i]]<-as.numeric(disease[[i]])


#Basic patient data: age, sex, BMI
patient.data <- metadata[c('sex', 'age.years','height.cm', 'weight.kg', 'bmi')]

#Discretize meta data and set probably erroneus values to NA
patient.data$sex[patient.data$sex=="male"] <- 1
patient.data$sex[patient.data$sex=="female"] <- 0
patient.data$sex[!(patient.data$sex %in% c(1, 0))] <- NA

for(i in 1:ncol(patient.data)) patient.data[,i] <- as.numeric(patient.data[[i]])

patient.data$height.cm[patient.data$height.cm<38] <- NA
patient.data$height.cm[patient.data$height.cm>220] <- NA
patient.data$height.cm <- round(patient.data$height.cm, 2)

patient.data$weight.kg[patient.data$weight.kg<=1] <- NA
patient.data$weight.kg[patient.data$weight.kg>250] <- NA
patient.data$weight.kg <- round(patient.data$weight.kg, 2)

patient.data$age.years[patient.data$age.years<1] <- 0
patient.data$age.years[patient.data$age.years>150] <- NA
patient.data$age.years<-round(patient.data$age.years)

#BMI ic calculated
patient.data$bmi <- patient.data$weight.kg/(patient.data$height.cm/100)^2

#Remove rows with NA in basic data
patient.data <- patient.data[rowSums(is.na(patient.data))==0,]


#Remove taxa records with too small total abundance. 
#There are several methods to establish the threshold, but all result in the number close to the number of taxa
taxa <- taxa[rowSums(taxa)>=ncol(taxa),]

#Normalize data
taxa <- taxa/rowSums(taxa)

#subset of samples for which we have taxonomy info available
common.rows <- intersect(rownames(taxa), rownames(patient.data))

taxa <- taxa[common.rows,]
patient.data <- patient.data[common.rows,]
disease <- disease[common.rows,]

#Remove diseases with less than 30 cases in the smaller class
disease <- disease[,colSums(disease==0)>30 & colSums(disease==1)>30]


#Remove taxa with less than 30 non-zero records
taxa <- taxa[,colSums(taxa>0)>30]

#Binarize taxa
taxa.binary<-0*taxa
for(i in 1:ncol(taxa)) taxa.binary[,i]<-(taxa[,i]>median(taxa[,i]))


#Save the data
save(taxa, taxa.binary, patient.data, disease, file='data_clean.RData')






