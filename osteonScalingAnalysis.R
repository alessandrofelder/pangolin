library(plyr)
library(smatr)
library(ggplot2)
library(ggrepel)
library(scales)

#figure out OS dependent path separator symbol                          
separator <- '/'
if(Sys.info()[['sysname']]=="Windows")
{separator <- '\\' }


#shapiro wilk test for normality
swTestP <- function(x, na.rm = FALSE)
{
  if(length(x)>2){
    st <- shapiro.test(log(x))
    p <- st$p.value
  }
  else{
      p<-1e-12
    }
  return(p)
}

#helper functions
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#overall limits: 0.0037, 0.27
#overall log limits: -5.6, -1.3
all.intact.osteon.areas <- c()

setUpDataFromCSVfile <- function(files, var.name) 
{
  #find all csv files in current folder
  n.files <- length(files)
  
  #set up containers for data
  quantities.of.interest = c("file.name","mean", "min", "max", "median", "swTestP")
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
      if(length(temp.data$Area)<minimumN) next()
      if(substrRight(current.files[current.file],10)=="Intact.csv")
      {
        all.intact.osteon.areas <<- c(all.intact.osteon.areas, temp.data$Area)  #"<<-" : scopes in R are weird!
      }
      
      for(current.column in 2:(n.columns.of.interest))
      {
        current.quantity <- quantities.of.interest[current.column]
        current.entry.name <- tolower(paste(as.name(var.name), as.name(current.quantity),sep="."))
        names(data.of.interest)[names(data.of.interest)==paste0("V", as.name(current.column))] <- current.entry.name
        current.quantity.function <- get(quantities.of.interest[current.column])
        data.of.interest[,current.entry.name][current.file] <- current.quantity.function(temp.data$Area, na.rm = TRUE) 
      }
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


#get binomials and masses
file.to.latin.map <- read.csv("/media/rvc_projects/Research_Storage/Doube_Michael/Felder/images/histomorphometry/sizeData/imageCorrespondenceMap.csv",header=TRUE,stringsAsFactors=FALSE)

excludeUncertainSpecies = FALSE #to exclude uncertain species, set this to true
excludeFullyAquaticSpecies = FALSE
excludeHumeri = FALSE
excludeFemora = FALSE
speciesMeasured=0
minimumN=0; #minimum number of osteons present for a specimen to be included
excludeUncertainSpeciesString = "-all-species"
if(excludeUncertainSpecies)
{
  excludeUncertainSpeciesString = "-uncertain-excluded"
}

if(excludeFullyAquaticSpecies)
{
  excludeUncertainSpeciesString = paste0(excludeUncertainSpeciesString,"-no-aquatic")
}

if(excludeFemora)
{
  excludeUncertainSpeciesString = paste0(excludeUncertainSpeciesString,"-only-humeri")
}
if(excludeHumeri)
{
  excludeUncertainSpeciesString = paste0(excludeUncertainSpeciesString,"-only-femora")
}
if(minimumN>0)
{
  excludeUncertainSpeciesString = paste0(excludeUncertainSpeciesString,"-minimumN-",minimumN)  
}

vars.of.interest <- vector()
for (current.input.path in 1:(length(input.paths)))
{ 
  current.var.name <- tail(as.vector(strsplit(input.paths[current.input.path],separator))[[1]],n=1)
  current.var.name <- strsplit(current.var.name, "(?<=[a-z])(?=[A-Z])", perl = TRUE)
  current.var.name <- paste(tolower(as.vector(current.var.name[[1]])), collapse='.')
  
  current.files <-list.files(path=input.paths[current.input.path], pattern="*.csv$", full.names=T, recursive=FALSE)
  current.data <- setUpDataFromCSVfile(current.files, current.var.name)
  if(excludeUncertainSpecies)
  {
    current.data <- current.data[-which(file.to.latin.map$bool.species.certain=="species uncertain"),]
  }
  if(excludeFullyAquaticSpecies)
  {
    current.data <- current.data[-which(file.to.latin.map$lifestyle=="aquatic"),]
  }
  if(excludeFemora)
  {
    current.data <- current.data[-which(file.to.latin.map$bone=="femur"),]
  }
  if(excludeHumeri)
  {
    current.data <- current.data[-which(file.to.latin.map$bone=="humerus"),]
  }
  #save a copy of data frame with appropriate name
  data.frame.name <- tolower(paste0(current.var.name, ".data"))
  assign(data.frame.name, current.data)
  vars.of.interest <- c(vars.of.interest, data.frame.name)
}


#derived variables
infilling.ratio.data <- infilled.areas.data[,seq(2,5)]/intact.osteons.data[,seq(2,5)]
infilling.ratio.data <- cbind(as.data.frame(intact.osteons.data$intact.osteons.file.name,stringsAsFactors = FALSE), infilling.ratio.data)
colnames(infilling.ratio.data) <- gsub("infilled.areas", "infilling.ratio", colnames(infilling.ratio.data))
vars.of.interest <- c(vars.of.interest, "infilling.ratio.data")

canal.cement.line.data <- as.data.frame(sqrt(intact.osteons.data[,seq(2,5)]/pi)-sqrt(haversian.canals.data[,seq(2,5)]/pi))
canal.cement.line.data <- cbind(as.data.frame(intact.osteons.data$intact.osteons.file.name,stringsAsFactors = FALSE), canal.cement.line.data)
colnames(canal.cement.line.data) <- gsub("intact.osteon", "canal.cement.line", colnames(canal.cement.line.data)) 
vars.of.interest <- c(vars.of.interest, "canal.cement.line.data")

surface.to.volume.ratio.data <- haversian.canals.data[,seq(2,5)]/infilled.areas.data[,seq(2,5)]
surface.to.volume.ratio.data <- cbind(as.data.frame(intact.osteons.data$intact.osteons.file.name,stringsAsFactors = FALSE), surface.to.volume.ratio.data)
colnames(surface.to.volume.ratio.data) <- gsub("haversian.canals", "surface.to.volume.ratio", colnames(infilling.ratio.data))
vars.of.interest <- c(vars.of.interest, "surface.to.volume.ratio.data")

image.names <- as.data.frame(matrix(ncol=2, nrow=length(intact.osteons.data$intact.osteons.file.name)))
names(image.names)[names(image.names)=="V1"] <- "image.name"
names(image.names)[names(image.names)=="V2"] <- "species.name"
image.names$image.name = gsub("_Intact.csv", ".tif", intact.osteons.data$intact.osteons.file.name)

latin.to.mass.map <- read.table("/media/rvc_projects/Research_Storage/Doube_Michael/Felder/images/histomorphometry/sizeData/PanTHERIA_1-0_WR05_Aug2008_alldata_mastodon_and_palaeoloxodon_extended_paca_corrected.csv",header=TRUE, row.names=NULL,stringsAsFactors=FALSE, sep=",")

mass.data <- vector(mode = "numeric", length = length(image.names$image.name))
abbreviated.names <- vector(mode = "character", length = length(image.names$image.name))

getLatinNameIndex <- function(image.name)
{
  index <- grep(basename(image.name),file.to.latin.map$image.name)
  return(index)
}

getMassFromLatinName <- function(latin.name)
{
  mass.index <- grep(latin.name,latin.to.mass.map$MSW05_Binomial)
  if(length(mass.index)!=1)
  {
    print(paste("species not in this database:",latin.name,"(",i,")"))
    mass <- NA
  }
  else
  {
    mass <- latin.to.mass.map$X5.1_AdultBodyMass_g[mass.index]/1000.0#in kg
  }
}

getSnoutVentLengthFromLatinName <- function(latin.name)
{
  mass.index <- grep(latin.name,latin.to.mass.map$MSW05_Binomial)
  if(length(mass.index)!=1)
  {
    print(paste("species not in this database:",latin.name))
    svlength <- NA
  }
  else
  {
    svlength <- latin.to.mass.map$X13.1_AdultHeadBodyLen_mm[mass.index]/1000.0#in kg
  }
}

for(i in 1:length(intact.osteons.data$intact.osteons.file.name))
{
  latin.name.index <- getLatinNameIndex(image.names$image.name[i])
  
  if(file.to.latin.map$species.name[latin.name.index]=="NULL")
  {
    print(paste("unclear species or missing latin name",i))
    image.names$species.name[i] <- NA
    abbreviated.names[i]<- NA
  } else {
    image.names$species.name[i] <- file.to.latin.map$species.name[latin.name.index]
    #print(image.names$species.name[i])
    abbreviated.names[i] <- abbreviate(image.names$species.name[i])
  }  
  mass.data[i] <- getMassFromLatinName(image.names$species.name[i])
  #mass.data[i] <- getSnoutVentLengthFromLatinName(image.names$species.name[i])
  #mass.data[i] <- utf8ToInt(substring(image.names$species.name[i],1,1))
}

#add 0.05 kg as average mass for Chiroptera species consistent with CSA of Quekett specimen
mass.data[which(image.names$species.name=="Chiroptera")] <- 0.05
out <- intact.osteons.data
out <- cbind(mass.data,image.names$species.name, out)
out <- cbind(out, haversian.canals.data[,2:5])
out <- cbind(out, infilled.areas.data[,2:5])
out <- cbind(out, infilling.ratio.data[,2:5])
out <- cbind(out, canal.cement.line.data[,2:5])
out <- cbind(out, surface.to.volume.ratio.data[,2:5])

write.csv(out,paste0("uncorrectedData",excludeUncertainSpeciesString,".csv"))

#run scaling analysis and plot results for med/min/max
options(digits = 2, scipen =3)

yLabels = c("osteon area ", "canal area ", "infill area ", "infill ratio ", "infill distance ", "border ratio ")
getYLabel <- function(var.index) {
#  yl <- paste(tolower(vars.of.interest[var.index]), collapse=' ')
#  yl <- gsub("haversian","Haversian",yl) 
#  yl <- gsub("canal\\.cement","canal border-cement",yl)
#  measure.string=""
#  if(var.index < 3){measure.string="\\.area"}
#  if(var.index == 5){measure.string="\\.distance"}
#  yl <- gsub("\\.data",measure.string,yl) 
#  yl <- gsub("\\."," ",yl) 

    return(yLabels[var.index])
}

getYUnit <- function(var.index) {
  y.units = c(rep("[mm2]",3),""," [mm]","")
  return(y.units[var.index])
}
hackyCounter<-0
plotScalingAnalysis <- function(fit.to.plot, input.data, input.labels, y.label, plot.title,file.name){
  if(hackyCounter==0||hackyCounter==4)
  {
    x.label <- ""
  }
  else
  {
    x.label <- "adult body mass [kg]"
  }
  hackyCounter<<-hackyCounter+1
  
  fit.coefficients <- fit.to.plot$coef[[1]]
  p.value.tolerance <- 0.005
  r2.value.tolerance <- 0.3
  adapted.y.label <- y.label;
  data.to.plot <- input.data[-which(is.na(input.data[,2])),]
  data.labels <- input.labels[-which(is.na(input.data[,2]))]
  xmin <- min(data.to.plot[,1])
  xmax <- max(data.to.plot[,1])
  axminpb <- 10^fit.coefficients$`coef(SMA)`[1]*xmin^fit.coefficients$`coef(SMA)`[2]
  axmaxpb <- 10^fit.coefficients$`coef(SMA)`[1]*xmax^fit.coefficients$`coef(SMA)`[2]
  if(substrRight(adapted.y.label, 4)=="mm2]")
  {
    adapted.y.label<-substr(adapted.y.label, 1, nchar(adapted.y.label)-3)
    adapted.y.label<-substitute(label*m^2*"]", list(label=adapted.y.label))
  }
  if(fit.to.plot$pval<p.value.tolerance && fit.to.plot$r2>r2.value.tolerance){
    ggplot(data.to.plot, aes(x=data.to.plot[,1],y=data.to.plot[,2]))+ geom_segment(aes(x = xmin, y = axminpb, xend = xmax, yend = axmaxpb)) + geom_point()  + labs(x = x.label, y = adapted.y.label) + theme_bw(base_size = 20) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_text_repel(label=data.labels,size=3,fontface = "italic",force=5)  + scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) 
  }else{
    ggplot(data.to.plot, aes(x=data.to.plot[,1],y=data.to.plot[,2])) + geom_point()  + labs(x = x.label, y = adapted.y.label) + theme_bw(base_size = 20) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_text_repel(label=data.labels,size=3,fontface = "italic",force=7) + scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x))
  }
 ggsave(file=file.name,scale=1,dpi=600)
}

n.vars.of.interest <- length(vars.of.interest)
csv.output <- as.data.frame(matrix(ncol = 4, nrow = 12))
rownames(csv.output) <- c("b", "b^{-}", "b^{+}","a", "a^{-}", "a^{+}","R^2", "p_{\\textit{uncorrelated}}", "p_{\\textit{isometry}}","p_{\\textit{shapiro-wilk}}", "\\rho^2","p_{\\textit{non-param.}}")

isometric_exponents <- c(rep(2.0/3.0,3),0.0,1.0/3.0,-1.0/3.0)
for(current.var.to.plot in 1:(n.vars.of.interest))
{
  var.to.plot <- eval(as.name(vars.of.interest[current.var.to.plot]))
  par(mfrow = c(1, 1))
  png.file.name <- gsub("\\.","-",vars.of.interest[current.var.to.plot])
  png.file.name <- paste0(png.file.name,excludeUncertainSpeciesString,"-mean.png")
  fit.mean   <- sma(var.to.plot[,2]~mass.data,log="xy",Robust=T, slope.test=isometric_exponents[current.var.to.plot])
  fit.mean$spearman.rho.squared <- (cor.test(log(var.to.plot[,2]),log(mass.data),method = "spearman")$estimate)^2
  fit.mean$spearman.cor.test.p <- (cor.test(log(var.to.plot[,2]),log(mass.data),method = "spearman")$p.value)
  fit.mean$shapiro.p <- shapiro.test(log(var.to.plot[,2]))$p.value
  plotScalingAnalysis(fit.mean,as.data.frame(cbind(mass.data,var.to.plot[,2])),input.labels = abbreviated.names,y.label = paste0("mean ",getYLabel(current.var.to.plot),getYUnit(current.var.to.plot)),plot.title="",file.name = png.file.name)
  
  png.file.name <- gsub("\\.","-",vars.of.interest[current.var.to.plot])
  png.file.name <- paste0(png.file.name,excludeUncertainSpeciesString,"-min.png")
  var.to.plot[,3] <- replace(var.to.plot[,3], is.infinite(var.to.plot[,3]),NA) #avoid Inf values in sma call
  fit.min <- sma(var.to.plot[,3]~mass.data,log="xy",Robust=T, slope.test=isometric_exponents[current.var.to.plot])
  fit.min$spearman.rho.squared <- (cor.test(log(var.to.plot[,3]),log(mass.data),method = "spearman")$estimate)^2
  fit.min$spearman.cor.test.p <- (cor.test(log(var.to.plot[,3]),log(mass.data),method = "spearman")$p.value)
  fit.min$shapiro.p <- shapiro.test(log(var.to.plot[,3]))$p.value
  plotScalingAnalysis(fit.min,as.data.frame(cbind(mass.data,var.to.plot[,3])),input.labels = abbreviated.names,y.label = paste0("minimum ",getYLabel(current.var.to.plot),getYUnit(current.var.to.plot)),plot.title="",file.name = png.file.name)
  
  png.file.name <- gsub("\\.","-",vars.of.interest[current.var.to.plot])
  png.file.name <- paste0(png.file.name,excludeUncertainSpeciesString,"-max.png")
  var.to.plot[,4] <- replace(var.to.plot[,4], is.infinite(var.to.plot[,4]),NA) #avoid Inf values in sma call
  fit.max    <- sma(var.to.plot[,4]~mass.data,log="xy",Robust=T, slope.test=isometric_exponents[current.var.to.plot])
  fit.max$spearman.rho.squared <- (cor.test(log(var.to.plot[,4]),log(mass.data),method = "spearman")$estimate)^2
  fit.max$spearman.cor.test.p <- (cor.test(log(var.to.plot[,4]),log(mass.data),method = "spearman")$p.value)
  fit.max$shapiro.p <- shapiro.test(log(var.to.plot[,4]))$p.value
  plotScalingAnalysis(fit.max,as.data.frame(cbind(mass.data,var.to.plot[,4])),input.labels = abbreviated.names,y.label = paste0("maximum ",getYLabel(current.var.to.plot),getYUnit(current.var.to.plot)),plot.title="",file.name = png.file.name)

  png.file.name <- gsub("\\.","-",vars.of.interest[current.var.to.plot])
  png.file.name <- paste0(png.file.name,excludeUncertainSpeciesString,"-median.png")
  fit.median <- sma(var.to.plot[,5]~mass.data,log="xy",Robust=T, slope.test=isometric_exponents[current.var.to.plot])
  fit.median$spearman.rho.squared <- (cor.test(log(var.to.plot[,5]),log(mass.data),method = "spearman")$estimate)^2
  fit.median$spearman.cor.test.p <- (cor.test(log(var.to.plot[,5]),log(mass.data),method = "spearman")$p.value)
  fit.median$shapiro.p <- shapiro.test(log(var.to.plot[,5]))$p.value
  plotScalingAnalysis(fit.median,as.data.frame(cbind(mass.data,var.to.plot[,5])),input.labels = abbreviated.names,y.label = paste0("median ",getYLabel(current.var.to.plot),getYUnit(current.var.to.plot)),plot.title="",file.name = png.file.name)
  

  #redirect text and graphic output
  fit.names = c("mean","min","max","median")
  n.fit.names = length(fit.names)
  current.range <- (1:n.fit.names)
  colnames(csv.output)[current.range] <- paste0(vars.of.interest[current.var.to.plot],".", fit.names)
  
  dev.off()
  sink(file = paste0(as.name(vars.of.interest[current.var.to.plot]),excludeUncertainSpeciesString,"-summary.txt"))
  png.file.name <- gsub("\\.","-",vars.of.interest[current.var.to.plot])
  png(file = paste0(png.file.name,excludeUncertainSpeciesString,"-errors.png"), width = 1000, height = 2000)
  
  par(mfrow = c(5, 2))
  
  for(current.fit.to.output in 1:n.fit.names)
  { 
    
    fit <- eval(as.name(paste0("fit.",fit.names[current.fit.to.output])))
    
    #assumption checking plots
    print(plot(main=fit.names[current.fit.to.output],fit, which="residual"))
    print(plot(main=fit.names[current.fit.to.output],fit, which="qq")) 
    #
    
    #output raw summary 
    summary(fit) 
    cat("\n################################################################################# \n\n")  
    
    #prepare formatted values
    csv.output.row <- unlist(c(fit$coef[[1]]$"coef(SMA)"[2],fit$coef[[1]]$"lower limit"[2], fit$coef[[1]]$"upper limit"[2],
                               10.0^fit$coef[[1]]$"coef(SMA)"[1],10.0^fit$coef[[1]]$"lower limit"[1], 10.0^fit$coef[[1]]$"upper limit"[1],
                               fit$r2,fit$pval, fit$slopetest[[1]]$p, fit$invariance.test.p, fit$shapiro.p, fit$spearman.rho.squared, fit$spearman.cor.test.p))
    
    csv.output.row.rounded <- signif(csv.output.row, 2)   
    csv.output[current.range[1]+current.fit.to.output-1] <- csv.output.row.rounded
  }
  dev.off()
  sink(NULL)
  #write formatted output
  write.csv(csv.output,paste0(as.name(vars.of.interest[current.var.to.plot]),excludeUncertainSpeciesString,"-table.csv"), quote = FALSE)
  
}

