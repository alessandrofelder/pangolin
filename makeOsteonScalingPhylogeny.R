library(ape)
library(picante)
library(geiger)
library(phytools)

setwd("/home/alessandro/Documents/data/phyloCorrection")

#Thanks to Liam Revell for writing this function (http://blog.phytools.org/2012/11/addendum-to-adding-single-tip-to-tree.html)
bind.tip<-function(tree,tip.label,edge.length=NULL, where=NULL,position=0){
  if(is.null(where)) where<-length(tree$tip)+1
  if(is.null(edge.length)&&is.ultrametric(tree)){
    H<-nodeHeights(tree)
    if(where==(length(tree$tip)+1))
      edge.length<-max(H)
    else 
      edge.length<-max(H)-H[tree$edge[,2]==where,2]+position
  }
  tip<-list(edge=matrix(c(2,1),1,2),
            tip.label=tip.label,
            edge.length=edge.length,
            Nnode=1)
  class(tip)<-"phylo"
  obj<-bind.tree(tree,tip,where=where,position=position)
  return(obj)
}

#Thanks to 10k project for providing a template for this function (http://10ktrees.fas.harvard.edu/exercises/10kTrees_Exercise3.pdf)
insertSpecies <- function(original.tree, species.name, sister.species.name, branch.length)
{ 
  speciesToBeInserted = species.name
  insertSpeciesNextTo = sister.species.name
  branchLengthSplitAgo = branch.length
  
  treeNewStr = paste("(", speciesToBeInserted, ":", branchLengthSplitAgo, ");", sep="")
  treeNew <- read.tree(text = treeNewStr)
  # Read tree block
  treeBlock = original.tree
  treeBlockNew = list()
  # For each tree in the tree block
  for (i in 1:length(treeBlock)) { 
    pos = 0
    for (j in 1: length(treeBlock$tip.label)) { if ( insertSpeciesNextTo == treeBlock$tip.label[j]) { pos = j
    }
      
    }
    if (pos == 0) { stop ("Error: Could not find species ",insertSpeciesNextTo, " in tree ",i,", exit\n")
    }
    #Merge both trees and write tree to new list of trees
    treeBlockNew = bind.tree(treeBlock,
                             treeNew,where = pos, position = branchLengthSplitAgo)
  } 
  return(treeBlockNew)
}

#phylogeny from Bininda-Emonds et al. 2007
tree<-read.nexus("./Bininda-emonds_2007_mammals.nex")
tree<-tree$mammalST_bestDates
#update binomials
tree$tip.label[which(tree$tip.label=="Arvicola_terrestris")]="Arvicola_amphibius"
tree$tip.label[which(tree$tip.label=="Agouti_paca")]="Cuniculus_paca"

#adapt binomials without changing the essence of the phylogeny
tree$tip.label[which(tree$tip.label=="Elephas_maximus")]="Palaeoloxodon_namadicus"
tree$tip.label[which(tree$tip.label=="Sphiggurus_mexicanus")]="Sphiggurus_melanurus" #watch out for the relationship of Erethizontidae, Cunuculidae and Hystricidae. Does not matter for our analysis, as porcupine specimen has no secondary osteons

#replace random bat species with "Chiroptera sp" as we have little idea of what the species was. Exact species will not make any difference to phyloCorrection as it is the only bat species in the sample
tree$tip.label[which(tree$tip.label=="Pteropus_conspicillatus")]="Chiroptera"

uncorrected.data<-read.csv("./uncorrectedData-all-species.csv",row.names=1)
binomials <- make.names(uncorrected.data$image.names.species.name,FALSE)
binomials <- gsub('\\.','_',binomials)
dummy.data <- as.matrix(rep(1,length(binomials)))
rownames(dummy.data) <- binomials

missing<-name.check(tree,dummy.data)$tree_not_data
pruned.tree<-drop.tip(tree, missing,root.edge = 1)
pruned.tree<-multi2di(pruned.tree)
pruned.tree <- insertSpecies(pruned.tree, "Mammut_americanum", "Palaeoloxodon_namadicus", 28.3)


png(file = "osteon-scaling-phylogenetic-tree.png", width = 1800, height = 1200)
plot(ladderize(pruned.tree))
dev.off()

write.tree(pruned.tree, "./osteon-scaling-phylogenetic-tree-uncertain-excluded.txt")


