library(ape)
library(picante)
library(geiger)
library(smatr)
library(phytools)

setwd("~/Documents/data/phyloCorrection/")
excludeUncertainSpecies = FALSE
excludeUncertainSpeciesString = "-all-species"
if(excludeUncertainSpecies)
{
  excludeUncertainSpeciesString = "-uncertain-excluded"
}

uncorrected.data<-read.csv(paste0("uncorrectedData",excludeUncertainSpeciesString,".csv"),row.names=1)
results <- matrix(nrow=ncol(uncorrected.data),ncol=8)
colnames(results) <- c("slope","lower limit","upper limit","P_value-tested-slope","elevation","r2","P_value-no-correlation", "P_value-no-scaling")


#average numerical entries for both Meles meles specimens
data.to.average <- uncorrected.data[which(uncorrected.data$image.names.species.name=="Meles meles"),]
average <- data.to.average[1,]
average[4:13] <- colMeans(rbind(average[4:13], data.to.average[2,4:13]), na.rm=TRUE)

#remove Meles meles lines and add averaged line
uncorrected.data <- uncorrected.data[-which(uncorrected.data$image.names.species.name=="Meles meles"),]
uncorrected.data <- rbind(uncorrected.data, average)

#average numerical entries for both Lutra lutra specimens
data.to.average <- uncorrected.data[which(uncorrected.data$image.names.species.name=="Lutra lutra"),]
average <- data.to.average[1,]
average[4:13] <- colMeans(rbind(average[4:13], data.to.average[2,4:13]), na.rm=TRUE)

#remove Lutra lutra lines and add averaged line
uncorrected.data <- uncorrected.data[-which(uncorrected.data$image.names.species.name=="Lutra lutra"),]
uncorrected.data <- rbind(uncorrected.data, average)

#average numerical entries for both Semnopithecus specimens
data.to.average <- uncorrected.data[which(uncorrected.data$image.names.species.name=="Semnopithecus entellus"),]
average <- data.to.average[1,]
average[4:13] <- colMeans(rbind(average[4:13], data.to.average[2,4:13]), na.rm=TRUE)

#remove Semnopithecus lines and add averaged line
uncorrected.data <- uncorrected.data[-which(uncorrected.data$image.names.species.name=="Semnopithecus entellus"),]
uncorrected.data <- rbind(uncorrected.data, average)


#average numerical entries for both Phocoena phocoena specimens
data.to.average <- uncorrected.data[which(uncorrected.data$image.names.species.name=="Phocoena phocoena"),]
average <- data.to.average[1,]
average[4:13] <- colMeans(rbind(average[4:13], data.to.average[2,4:13]), na.rm=TRUE)

#remove Phocoena phocoena lines and add averaged line
uncorrected.data <- uncorrected.data[-which(uncorrected.data$image.names.species.name=="Phocoena phocoena"),]
uncorrected.data <- rbind(uncorrected.data, average)

#average numerical entries for all Mustela putorius specimens
data.to.average <- uncorrected.data[which(uncorrected.data$image.names.species.name=="Mustela putorius"),]
average <- data.to.average[1,]
average[4:13] <- colMeans(rbind(average[4:13], data.to.average[2,4:13], data.to.average[3,4:13]), na.rm=TRUE)

#remove Mustela putorius lines and add averaged line
uncorrected.data <- uncorrected.data[-which(uncorrected.data$image.names.species.name=="Mustela putorius"),]
uncorrected.data <- rbind(uncorrected.data, average)

#set up scaling exponents
isometricExponents = c(rep(2.0/3.0,4),1,rep(2.0/3.0,4),1,rep(2.0/3.0,4),1,rep(1,5),rep(1.0/3.0,4),1)

tree <- read.tree("osteon-scaling-phylogenetic-tree.txt")
mass<-data.frame(uncorrected.data$mass.data)

for(current.Y.variable in c(3:23)) {
  test.column<-as.matrix(as.numeric(uncorrected.data[,current.Y.variable]))
  binomials <- make.names(uncorrected.data$image.names.species.name, FALSE)
  binomials <- gsub('\\.','_',binomials)
  rownames(mass)<-binomials
  rownames(test.column)<-binomials
  
  test.column<-replace(test.column, is.infinite(log(test.column)),NA)
  regression.data<-cbind(log(mass), log(test.column))
  regression.data<-regression.data[complete.cases(regression.data),]
  regression.data<-as.matrix(regression.data)
  
  missing<-name.check(tree,regression.data)
  pruned.tree<-tree
  
  if(missing[1]!="OK")
  {
    missing<-name.check(tree,regression.data)$tree_not_data
    pruned.tree<-drop.tip(tree, missing,root.edge = 1)
  }
  
  #Phylogenetic independent contrasts for testcol and mass against the tree
  contrasts.data <- apply(as.matrix(regression.data[,2]), 2, pic, pruned.tree)
  contrasts.mass <- apply(as.matrix(regression.data[,1]), 2, pic, pruned.tree)
  
  combined<-cbind(contrasts.mass, contrasts.data)
  contrasts.data<-data.frame(combined)
  
  # option -1 forces fitted line through origin
  contrasts.analysis <- sma(contrasts.data[,2]~contrasts.data[,1]-1,slope.test=isometricExponents[current.Y.variable],data=contrasts.data)
  
  results[current.Y.variable,"slope"] <- contrasts.analysis$coef[[1]]["slope","coef(SMA)"]
  results[current.Y.variable,"r2"] <- contrasts.analysis$r2[[1]]
  results[current.Y.variable,"P_value-tested-slope"] <- contrasts.analysis$slopetest[[1]]$p
  results[current.Y.variable,"elevation"] <- contrasts.analysis$coef[[1]]["elevation","coef(SMA)"]
  results[current.Y.variable,"lower limit"] <- contrasts.analysis$coef[[1]]["slope","lower limit"]
  results[current.Y.variable,"upper limit"] <- contrasts.analysis$coef[[1]]["slope","upper limit"]
  results[current.Y.variable,"P_value-no-correlation"] <- contrasts.analysis$pval[[1]]
  
  contrasts.analysis <- sma(contrasts.data[,2]~contrasts.data[,1]-1,slope.test=0.0,data=contrasts.data)
  results[current.Y.variable,"P_value-no-scaling"] <- contrasts.analysis$slopetest[[1]]$p
  
  
  png(file = paste0("osteon-scaling-",colnames(uncorrected.data[current.Y.variable]),excludeUncertainSpeciesString,"-phy-indep-contrasts.png"), width = 600, height = 400)
  plot(contrasts.data[,2]~contrasts.mass, main=colnames(uncorrected.data[current.Y.variable]))
  dev.off()
}

rownames(results)=rep("", nrow(results))
rownames(results)[1]="mass"
rownames(results)[2]="binomial"
rownames(results)[3]="file name"
measurements = c("mean","minimum","maximum","median")
quantities = c("osteon area","canal area","infill area","infill ratio","infill distance")
for(q in c(1:length(quantities)))
{
  for(m in c(1:length(measurements)))
  {
    rownames(results)[m+length(measurements)*(q-1)+3] = paste(measurements[m],quantities[q])          
  }
}
 
results
write.csv(format(results,digits=3, scientific = FALSE), file=paste0("osteon-scaling-phylogenetic-regression",excludeUncertainSpeciesString,".csv"),quote=FALSE)

