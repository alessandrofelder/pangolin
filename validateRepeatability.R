library(plyr)

#set up validation data
setwd("/media/rvc_projects/Research_Storage/Doube_Michael/Felder/images/histomorphometry/Quekett/Quekett Collection/Autofluorescence/Analysis/ProcessedData/Validation-Intact-Osteons/")
#figure out OS dependent path separator symbol                          
separator <- '/'
if(Sys.info()[['sysname']]=="Windows")
{separator <- '\\' }

#get only required subdirectories, which start with "run"
sub.directories = list.dirs('.', recursive=FALSE)
sub.directories = sub.directories[-which(substr(sub.directories,1,5)!="./run")]

hippo.data <- data.frame()
paca.data <- data.frame()
puma.data <- data.frame()

for (current.sub.directory in 1:(length(sub.directories)))
{ 
  current.files <-list.files(path=sub.directories[current.sub.directory], pattern="*.csv$", full.names=T, recursive=FALSE)

  temp.hippo.data <- try(read.csv(current.files[1], header=T,stringsAsFactors=F))
  if(!inherits(temp.hippo.data,"try-error"))
  {
    temp.hippo.data <- as.data.frame(temp.hippo.data$Area)
    colnames(temp.hippo.data) <- paste0("run",current.sub.directory) 
    hippo.data <- rbind.fill(hippo.data, as.data.frame(temp.hippo.data))
  } 
  
  temp.paca.data <- try(read.csv(current.files[2], header=T,stringsAsFactors=F))
  if(!inherits(temp.paca.data,"try-error"))
  {
    temp.paca.data <- as.data.frame(temp.paca.data$Area)
    colnames(temp.paca.data) <- paste0("run",current.sub.directory) 
    paca.data <- rbind.fill(paca.data, temp.paca.data)
  } 
  
  temp.puma.data <- try(read.csv(current.files[3], header=T,stringsAsFactors=F))
  if(!inherits(temp.puma.data,"try-error"))
  {
    temp.puma.data <- as.data.frame(temp.puma.data$Area)
    colnames(temp.puma.data) <- paste0("run",current.sub.directory) 
    puma.data <- rbind.fill(puma.data, temp.puma.data)
  } 
}

#validate puma
puma.run1 <- puma.data[-which(is.na(puma.data[,1])),1]
puma.run2 <- puma.data[-which(is.na(puma.data[,2])),2]
puma.run3 <- puma.data[-which(is.na(puma.data[,3])),3]
puma.combined.data <- c(puma.run1, puma.run2, puma.run3)
puma.factors <- factor(rep(1:3, c(length(puma.run1),length(puma.run2),length(puma.run3))),labels=c("run 1","run 2","run 3"))

boxplot(puma.combined.data~puma.factors,main="Bd 258")
kruskal.test(puma.combined.data,puma.factors)

shapiro.test(puma.run1)
shapiro.test(puma.run2)
shapiro.test(puma.run3)

puma.fit = lm(puma.combined.data~puma.factors)
anova(puma.fit)


#validate paca
paca.run1 <- paca.data[-which(is.na(paca.data[,1])),1]
paca.run2 <- paca.data[-which(is.na(paca.data[,2])),2]
paca.run3 <- paca.data[-which(is.na(paca.data[,3])),3]
paca.combined.data <- c(paca.run1, paca.run2, paca.run3)
paca.factors <- factor(rep(1:3, c(length(paca.run1),length(paca.run2),length(paca.run3))),labels=c("run 1","run 2","run 3"))

boxplot(paca.combined.data~paca.factors,main="Bd 37")
kruskal.test(paca.combined.data,paca.factors)

shapiro.test(paca.run1)
shapiro.test(paca.run2)
shapiro.test(paca.run3)

paca.fit = lm(paca.combined.data~paca.factors)
anova(paca.fit)


#validate hippo
hippo.run1 <- hippo.data[-which(is.na(hippo.data[,1])),1]
hippo.run2 <- hippo.data[-which(is.na(hippo.data[,2])),2]
hippo.run3 <- hippo.data[-which(is.na(hippo.data[,3])),3]
hippo.combined.data <- c(hippo.run1, hippo.run2, hippo.run3)
hippo.factors <- factor(rep(1:3, c(length(hippo.run1),length(hippo.run2),length(hippo.run3))),labels=c("run 1","run 2","run 3"))

boxplot(hippo.combined.data~hippo.factors, main="Bd 167")
kruskal.test(hippo.combined.data,hippo.factors)

shapiro.test(hippo.run1)
shapiro.test(hippo.run2)
shapiro.test(hippo.run3)

hippo.fit = lm(hippo.combined.data~hippo.factors)
anova(hippo.fit)
