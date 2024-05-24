#!/usr/bin
library(mclust)
library(RColorBrewer)
library(scales)
input=commandArgs(trailing=TRUE)

functional_directory=input[1]
patientID=input[2]

#print(functional_directory)

selectMclustG <- function(dist){

if(length(dist) < 10000) {
 samplesize <- length(dist)
}else { samplesize = 10000}

newdist <- sample(x=dist,size=samplesize)

Distance = as.vector(newdist) # 'dist' is a distance matrix

colors = brewer.pal(9, 'Set1') # color set
br = ceiling(max(Distance) - min(Distance)) # number of breaks
xmax = ceiling(max(Distance)) # distance upper bound

dens = densityMclust(Distance,modelNames="E")

bic.diff = max(dens$BIC, na.rm=T) - min(dens$BIC, na.rm=T)
bic = data.frame(dens$BIC[,])
bic$diff = 0
thresh =10 # threshold of 10% BIC loss
g = 1

for (n in 2:nrow(bic)) {
	bic$diff[n] = (bic[n,1] - bic[n-1,1]) * 100 / bic.diff 
	if( bic$diff[n] <= thresh ) {
		g = n-1
		break } }
#Don't let single solutions exist
if(g==1){g <- 3}

dens1 = densityMclust(Distance, modelNames='E', G=g)

n = dens1$n   # number of observations
classes = dens1$classification   # classification of data points
means = dens1$parameters$mean	# distribution means
vars = dens1$parameters$variance	# distribution variances
probs = dens1$parameters$pro	# probability of data points classes

#par(mfrow=c(2,1))

plot(dens, data=Distance, what='BIC')
points(x=g, y=bic[g,1], pch=24, bg='black')
            
plot(dens1, data=Distance, what='density', breaks=br, xlim=c(0, xmax))
points(x=Distance, y=rep(0, n), col=alpha(colors[classes], .3), pch=20)
text(x=means, y=rep(.005, g), col=colors[1:g], labels=as.character(round(means, 1)))

return(dens1)
}

filelocation <- list.files(path = functional_directory,pattern = ".calcvarreadprob.txt",full.names = T)

print(filelocation)

filelocation_shortened <- list.files(path = functional_directory,pattern = ".calcvarreadprob.txt",full.names = F)
filelocation_shortened <- gsub(pattern = ".by.singlesample.vcf.gz.temp.PairTree.txt.calcvarreadprob.txt",replacement = "", x = filelocation_shortened )

#Colors for plotting
colors = brewer.pal(9, 'Set1')

#Lit of locations
plat <- toString(unlist(filelocation_shortened))
write.table(x=plat,file=paste(patientID,"temporary.json",sep="."),sep="\t",row.name=F,col.name=F)


list_of_files <- lapply(filelocation,read.table,sep="\t",header=F,stringsAsFactors=FALSE)

subclonal_variants <- list()

pdf(file=paste(patientID,"_mclustVAF_selection.pdf",sep=""),width=12,height=10)
for (variable in 1:length(list_of_files)) {
samplename <- filelocation_shortened[variable]
t1 <- list_of_files[[variable]]
t1$V1 <- as.character(t1$V1)
t1$VAF <- t1$V2/(t1$V3)
t1 <- t1[t1$VAF!="NaN",]
t1$VAF <- as.numeric(t1$VAF)

par(mfrow=c(3,1))
denser <- selectMclustG(dist = t1$VAF)


#Find the heterozygous clusters and homozygous cluster
hetcluster_index <- ceiling(length(names(table(denser$classification)))/2)
hetcluster <- names(table(denser$classification)[hetcluster_index])
homocluster <- denser$G

res <- cbind(denser$data,denser$classification,denser$uncertainty,rep(0,length(denser$data)),rep(0,length(denser$data)))
res <- as.data.frame(res)

#Remove Factors as struture
#res$V1 <- as.character(res$V1)
res$V1 <- as.numeric(as.character(res$V1))
res$V2 <- as.numeric(as.character(res$V2))
res$V3 <- as.numeric(as.character(res$V4))
res$V4 <- as.numeric(as.character(res$V4))
res$V5 <- as.numeric(as.character(res$V5))

for(i in 1:length(res[,2])){
  clust <- res[i,2]
  res[i,4] <- mean(res[which(res[,2]==clust),1])
  res[i,5] <- median(res[which(res[,2]==clust),1])
}

#Set some basic rules for mclust density function. 
##If there is 
if (denser$G<=2){
  #print(paste("The maximum VAF is:"))
  exclusionVAF <- max(res[res$V2==1, ]$V1)
}else if (denser$G==3 && length(unique(res$V2))==3){
  exclusionVAF <- min(res[res$V2==2, ]$V1)
}else if (denser$G==3 && length(unique(res$V2))==2){
  exclusionVAF <- max(res[res$V2==1, ]$V1)
}else if (denser$G>=5){
  exclusionVAF <- min(res[res$V2==hetcluster,]$V1)
}else {exclusionVAF <- unique(res[res$V2==hetcluster,]$V4)}

#Deal with mclusts failure to assign proper classifications
if(exclusionVAF == -Inf) {exclusionVAF <- max(res[res$V2==2, ]$V1)}

#str(exclusionVAF)
#Set a maximum limit of VAF's for non hypermutated samples
if (exclusionVAF >= 0.30) {exclusionVAF <- 0.30} 

#Set minimum limit of VAF to avoid underfiltering (mainly for hypermutated)
if(exclusionVAF <= 0.20) {exclusionVAF <- 0.20}

VariantList_subclonalRemoved <- as.character(t1[t1$VAF >= exclusionVAF, ]$V1)
#subclonal_variants <- rbind(subclonal_variants,VariantList_subclonalRemoved)


if(length(VariantList_subclonalRemoved) >= 50000){
  denser_internal <- densityMclust(data = sample(t1$VAF,size = 20000), G=4:9,modelnames="E",plot = F)  
  res <- cbind(denser_internal$data,denser_internal$classification,denser_internal$uncertainty,rep(0,length(denser_internal$data)),rep(0,length(denser_internal$data)))
  res <- as.data.frame(res)
  res$V1 <- as.numeric(as.character(res$V1))
  res$V2 <- as.numeric(as.character(res$V2))
  res$V3 <- as.numeric(as.character(res$V3))
  res$V4 <- as.numeric(as.character(res$V4))
  res$V5 <- as.numeric(as.character(res$V5))
  exclusionVAF <- min(res[res$V2==denser_internal$G-1,]$V1)
  if(exclusionVAF <= 0.30) {exclusionVAF <- 0.30}
  #if(exclusionVAF >= 0.50) {exclusionVAF <- 0.50}
}

exclusionVAF <- round(x=exclusionVAF,digits = 3)


#For double checking the number of variants in hypermutated spots
VariantList_subclonalRemoved <- as.character(t1[t1$VAF >= exclusionVAF, ]$V1)

#str(VariantList_subclonalRemoved)
#This is how many variants each sample will consider as clonal (above threshold)
print(length(VariantList_subclonalRemoved))

subclonal_variants <- append(subclonal_variants,VariantList_subclonalRemoved)

#Set out the 
#head(colors,n=denser$G)

#par(mfrow=c(3,1))
plot.Mclust(denser,what = "classification",xlim=c(0,1),col=head(colors,n=denser$G) ,xlab = paste("Het cluster=",hetcluster,"Homo-cluster=",homocluster," & exclusion VAF=", exclusionVAF,"From samples:",samplename, "Variants above threshold=", length(VariantList_subclonalRemoved),sep = " "))
#plot.Mclust(denser,what = "BIC")
#plot.Mclust(denser,what = "density",xlab = paste("Sample Name: ",samplename,sep=""))

print(paste("Exclusion VAF is: ",exclusionVAF," From samples: ",samplename,sep=""))
}
dev.off()

#VariantList_subclonalRemoved <- as.character(VariantList_subclonalRemoved)
subclonal_variants <- as.character(subclonal_variants)
subclonal_variants <- gsub(" ", "", subclonal_variants, fixed = TRUE)
print(length(subclonal_variants))


total <-  Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "V1", all.x = FALSE),list_of_files)
total$V1 <- gsub(" ", "", total$V1, fixed = TRUE)


#Alt-count & total_depth per site & variant_read_probability
colNumber_altCount <- total[,grep(pattern = "V2",x=colnames(total))]
colNumber_totalCount <- total[,grep(pattern = "V3",x=colnames(total))]
colNumber_var_read_prob <- total[,grep(pattern = "V4",x=colnames(total))]

#Replace 0 with 0.01 read probability to avoid divide by zero problems
colNumber_var_read_prob[colNumber_var_read_prob==0]=0.01

#Collapse into appropriate format
if(is.null(dim(colNumber_altCount))==FALSE){
colNumber_altCount <- data.frame(apply(colNumber_altCount, 1, paste, collapse=","), stringsAsFactors=FALSE)
colNumber_totalCount <- data.frame(apply(colNumber_totalCount, 1, paste, collapse=","), stringsAsFactors=FALSE)
colNumber_var_read_prob <- data.frame(apply(colNumber_var_read_prob, 1, paste, collapse=","), stringsAsFactors=FALSE)
}

newCompleted_varRead <- cbind(as.character(total$V1),colNumber_altCount,colNumber_totalCount,colNumber_var_read_prob)
newCompleted_varRead <- data.frame(newCompleted_varRead,stringsAsFactors=FALSE)
colnames(newCompleted_varRead) <- c("VariantName","alt_reads","total_read","var_read_prob")
newCompleted_varRead$VariantName <- as.character(newCompleted_varRead$VariantName)

#Keep only variants from the completed table that are clonal in one or more samples.
NewCompleted <- newCompleted_varRead[newCompleted_varRead$VariantName %in% unique(subclonal_variants), ]

CommonVariants <- read.table(file="/home/ahgillmo/master_scripts_slurm/pairtree_scripts/Tumor_cellLine_Xenograft_common_mutations.tsv",stringsAsFactors = F)
NewCompleted <- NewCompleted[!NewCompleted$VariantName %in% CommonVariants$V1, ]

write.table(x=NewCompleted, file =paste(patientID,".ssm",sep=""),quote = F,col.names = F,row.names = F)
 
