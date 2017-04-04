library(plyr)
library(smatr)

#figure out OS dependent path separator symbol                          
separator <- '/'
if(Sys.info()[['sysname']]=="Windows")
{separator <- '\\' }

setUpDataFromCSVfile <- function(files, var.name) 
{
  #find all csv files in current folder
  n.files <- length(files)
  
  #set up containers for data
  quantities.of.interest = c("file.name","mean", "min", "max", "median","length")
  n.columns.of.interest <- length(quantities.of.interest)
  data.of.interest <- as.data.frame(matrix(ncol=n.columns.of.interest, nrow=n.files)) 
  
  #set first column to filename
  data.of.interest[names(data.of.interest)=="V1"] <- files
  current.entry.name <-  tolower(paste(as.name(var.name), "file.name",sep="."))
  names(data.of.interest)[names(data.of.interest)=="V1"] <- current.entry.name
  
  #read input data
  for(current.file in 1:n.files)
  {
    temp.data <- try(read.csv(current.files[current.file], header=T,stringsAsFactors=F))
    if(!inherits(temp.data,"try-error"))
    {
      for(current.column in 2:(n.columns.of.interest-1))
      {
        current.quantity <- quantities.of.interest[current.column]
        current.entry.name <- tolower(paste(as.name(var.name), as.name(current.quantity),sep="."))
        names(data.of.interest)[names(data.of.interest)==paste0("V", as.name(current.column))] <- current.entry.name
        current.quantity.function <- get(quantities.of.interest[current.column])
        data.of.interest[,current.entry.name][current.file] <- current.quantity.function(temp.data$Area, na.rm = TRUE) 
      }
      
      #treat length differently atm, since length takes different number of arguments
      current.entry.name <-  tolower(paste(as.name(var.name), "length",sep="."))
      names(data.of.interest)[names(data.of.interest)==paste0("V", as.name(n.columns.of.interest))] <- current.entry.name
      data.of.interest[,current.entry.name][current.file] <- length(temp.data$Area) 
    }
  }
  return(data.of.interest)
}

#set up all data: vars designates Osteons/Canals/etc
setwd("/media/rvc_projects/Research_Storage/Doube_Michael/Felder/images/histomorphometry/Quekett/Quekett Collection/Autofluorescence/Analysis/ProcessedData/run-through-twice-more/")
input.paths <- c(
  "IntactOsteons",
  "HaversianCanals",
  "InfilledAreas"
)  

vars.of.interest <- vector()
for (current.input.path in 1:(length(input.paths)))
{ 
  current.var.name <- tail(as.vector(strsplit(input.paths[current.input.path],separator))[[1]],n=1)
  current.var.name <- strsplit(current.var.name, "(?<=[a-z])(?=[A-Z])", perl = TRUE)
  current.var.name <- paste(tolower(as.vector(current.var.name[[1]])), collapse='.')
  
  current.files <-list.files(path=input.paths[current.input.path], pattern="*.csv$", full.names=T, recursive=FALSE)
  current.data <- setUpDataFromCSVfile(current.files, current.var.name)
  
  #save a copy of data frame with appropriate name
  data.frame.name <- tolower(paste0(current.var.name, ".data"))
  assign(data.frame.name, current.data)
  vars.of.interest <- c(vars.of.interest, data.frame.name)
}

#derived variables
infilling.ratio.data <- infilled.areas.data[,seq(2,6)]/intact.osteons.data[,seq(2,6)]
infilling.ratio.data <- cbind(as.data.frame(intact.osteons.data$intact.osteons.file.name,stringsAsFactors = FALSE), infilling.ratio.data)
colnames(infilling.ratio.data) <- gsub("infilled.areas", "infilling.ratio", colnames(infilling.ratio.data))
vars.of.interest <- c(vars.of.interest, "infilling.ratio.data")

canal.cement.line.data <- as.data.frame(sqrt(intact.osteons.data[,seq(2,6)]/pi)-sqrt(haversian.canals.data[,seq(2,6)]/pi))
canal.cement.line.data <- cbind(as.data.frame(intact.osteons.data$intact.osteons.file.name,stringsAsFactors = FALSE), canal.cement.line.data)
colnames(canal.cement.line.data) <- gsub("intact.osteon", "canal.cement.line", colnames(canal.cement.line.data)) 
vars.of.interest <- c(vars.of.interest, "canal.cement.line.data")

out <- intact.osteons.data
out <- cbind(out, haversian.canals.data[,2:6])
out <- cbind(out, infilled.areas.data[,2:6])
out <- cbind(out, infilling.ratio.data[,2:6])
out <- cbind(out, canal.cement.line.data[,2:6])
