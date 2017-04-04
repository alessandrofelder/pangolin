library(smatr)

#set working directory
setwd("C://Dropbox//RVC//R stuff")

#read and clean data
comparative.osteoblast.data <- read.csv("fakeComparativeCellData.csv",header=TRUE, stringsAsFactors = TRUE)
comparative.osteoblast.data$specimenID <- as.character(comparative.osteoblast.data$specimenID)



options(digits = 2, scipen =3)

getYUnit <- function(var.index) {
  y.units = c(rep("[mm2]",3),"","[mm]","")
  return(y.units[var.index])
}

variables.of.interest <- c("doubling.time", "basal.alkphos", "mineral.alkphos", "percent.mineralised.area")
n.variables.of.interest <- length(variables.of.interest)

csv.output <- c("","b","b-","b+","a","a-","a+","R2","p")
for(current.variable.number in 1:n.variables.of.interest)
{
  #perform scaling analysis
  current.variable <- eval(parse(text = paste0("comparative.osteoblast.data$",variables.of.interest[current.variable.number])))
  fit <- sma(current.variable~comparative.osteoblast.data$specimen.body.mass,Robust=T)
  
  #plot results and errors
  plot(fit, xlab = "body mass [kg]", ylab = variables.of.interest[current.variable.number])
  plot(fit, which="residual")
  plot(fit, which="qq")
  
  #output stats
  print(variables.of.interest[current.variable.number])
  summary(fit) 
  
  #save stats
  csv.output.column <- unlist(c(variables.of.interest[current.variable.number],fit$coef[[1]]$"coef(SMA)"[2],fit$coef[[1]]$"lower limit"[2], fit$coef[[1]]$"upper limit"[2],
                             fit$coef[[1]]$"coef(SMA)"[1],fit$coef[[1]]$"lower limit"[1], fit$coef[[1]]$"upper limit"[1],
                             fit$r2,fit$pval))
  
  csv.output.column.rounded <- csv.output.column
  csv.output.column.rounded[2:9] <- signif(as.numeric(csv.output.column.rounded[2:9]), 5)   
  csv.output <- cbind(csv.output, csv.output.column.rounded)
}
colnames(csv.output) <- c(rep("", length(colnames(csv.output))))
write.csv(csv.output,"osteoblast-statistics.csv", row.names = FALSE)   

