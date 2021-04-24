##################################################################
#This script can determine the tree toplogy of each gene family, 
#whose bootstrap>60 is considered as a resolved gene tree.
#Then the resolved and unresolved gene tree will be visualized.
##################################################################

library("ggplot2")
library("ggtree")
library("ape")
library("methods")
library("phangorn")

scan_ancestor_for_topo<-function(rawtree,tdf,cutoff,genomes)
{
  tag=0
  node.C12 <- tdf[which(tdf$label %in% genomes[which(genomes[,2] %in% c("C1","C2")),1]),"node"]
  node.C13 <- tdf[which(tdf$label %in% genomes[which(genomes[,2] %in% c("C1","C3")),1]),"node"]
  node.C23 <- tdf[which(tdf$label %in% genomes[which(genomes[,2] %in% c("C2","C3")),1]),"node"]
  if(31 %in% node.C13 || 31 %in% node.C12 || 31 %in% node.C23) ##root node was traced, outgroup clustered into P1
  {
    return(tag)
  }else{
    LAS_P12<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2] %in% c("C1","C2")),1]),"node"])
    LAS_P13<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2] %in% c("C1","C3")),1]),"node"])
    LAS_P23<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2] %in% c("C2","C3")),1]),"node"])
    ###if two groups cluster together, 
    ###MCAS shared by A and outgroup should be the same with that shared by B and outgroup
    if(LAS_P12==LAS_P13 && LAS_P12==LAS_P23 && LAS_P13==LAS_P23)
    {
      tag = 0 
      return(tag)
    }
    if(LAS_P12==LAS_P13)##P1 as out
    {  ###from up to bottom, get decendants of each population
      LAS_P1<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2]=="C1"),1]),"node"])
      ban1=as.numeric(tdf[which(tdf[,2]==LAS_P1),4])
      bcan23=as.numeric(tdf[which(tdf[,2]==LAS_P23),4])
      if(is.na(bcan23) || is.na(ban1) || bcan23<cutoff || ban1<cutoff) #####in case that no bootstrap value was attached
      {tag=0} else if(bcan23>=cutoff && ban1>=cutoff)
      {tag=3}
      LAS_P2<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2]=="C2"),1]),"node"])
      LAS_P3<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2]=="C3"),1]),"node"])
      if(((LAS_P2 == LAS_P23)|| (LAS_P3 == LAS_P23)) && tag!=3){tag=0} #mixed topology but not the mix between clades entirely
      return(tag)
    }
    else if(LAS_P12==LAS_P23) ##P2 as out
    {
      LAS_P2<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2]=="C2"),1]),"node"])
      ban2=as.numeric(tdf[which(tdf[,2]==LAS_P2),4])
      bcan13=as.numeric(tdf[which(tdf[,2]==LAS_P13),4])
      if(is.na(bcan13) || is.na(ban2) || bcan13<cutoff || ban2<cutoff) #####in case that no bootstrap value was attached
      {tag=0} else if(bcan13>=cutoff && ban2>=cutoff)
      {tag=2}
      LAS_P1<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2]=="C1"),1]),"node"])
      LAS_P3<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2]=="C3"),1]),"node"])
      if(((LAS_P1 == LAS_P13)||(LAS_P3 == LAS_P13)) && tag!=2){tag=0}
      return(tag)
    }
    else if(LAS_P13==LAS_P23) ##P3 as out
    {
      LAS_P3<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2]=="C3"),1]),"node"])
      ban3=as.numeric(tdf[which(tdf$node==LAS_P3),"label"])
      bcan12=as.numeric(tdf[which(tdf$node==LAS_P12),"label"])
      if(is.na(bcan12) || is.na(ban3) || bcan12<cutoff || ban3<cutoff) #####in case that no bootstrap value was attached
      {tag=0} else if(bcan12>=cutoff && ban3>=cutoff)
      {tag=1}
      LAS_P1<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2]=="C1"),1]),"node"])
      LAS_P2<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2]=="C2"),1]),"node"])
      if(((LAS_P1 == LAS_P12) || (LAS_P2 == LAS_P12)) && tag!=1){tag=0}
      return(tag)
    }
    else
    {
      tag = 0 
      return(tag)
    }
  }
  
}
setwd("/path/of/all/rooted/gene/trees")
# read *nwk file and lables of the tree
tree_dir <- "/path/of/all/rooted/gene/trees"
TREE <- list.files(tree_dir,pattern = "*.nwk", full.names = T)

P12_tree<-list()
P23_tree<-list()
P13_tree<-list()
PNA_tree<-list()

P12_fam<-c()
P23_fam<-c()
P13_fam<-c()
PNA_fam<-c()

#conf <- read.table("conflict.txt")[,1]
#TREE <- paste(conf,".rooted.nwk",sep = "")
for(i in seq.int(TREE))#seq.int(TREE[subTREE]))#
{
  LABLE<-sub("rooted.nwk","labels",TREE[i],perl = T)    # get tip labels 
  genomes <- read.table(LABLE,sep="\t",stringsAsFactors = F)
  
  ##shorten name
  shorten <- function(x)
  {
    y<-unlist(strsplit(x,"_"))
    paste(y[1],y[2],y[3],sep="_")
  }
  genomes[,1] <- sapply(genomes[,1],shorten)
  
  colnames(genomes)<-c("gene","pop","hbt")
  
  title = sub(".rooted.nwk","",basename(TREE[i]))
  rawtree <- read.tree(TREE[i])
  rawtree$tip.label <- sapply(rawtree$tip.label,shorten)
  tree.df<-as.data.frame(fortify(rawtree))
  rownames(genomes)<-genomes[,1]
  
  
  if(scan_ancestor_for_topo(rawtree,tree.df,60,genomes) == 1)
  {P12_fam<-c(P12_fam,title)
  #p<-ggtree(rawtree,branch.length = "none") %<+% LS + 
  p<-ggtree(rawtree) %<+% genomes +
    geom_point(aes(shape = pop,color = hbt), size=2.5) + 
    geom_point2(aes(label=label,subset=!isTip & as.numeric(label) > 60),size=1.5,shape=20)  ##change the btsp as binary point
    geom_text2(aes(label = label,subset = !isTip),
               nudge_x = -0.003,nudge_y = 0.5,size = 2.5) +
    geom_tiplab(size=1.5 , aes(color = hbt)) +
    theme(legend.position ="left")+
    ggtitle(title)
  P12_tree[[length(P12_tree)+1]]<-p}
  
  if(scan_ancestor_for_topo(rawtree,tree.df,60,genomes) == 2)
  {P13_fam<-c(P13_fam,title)
  #p<-ggtree(rawtree,branch.length = "none") %<+% LS + 
  p<-ggtree(rawtree) %<+% genomes +
    geom_point(aes(shape = pop,color = hbt), size=2.5) + 
    geom_text2(aes(label = label,subset = !isTip),
               nudge_x = -0.003,nudge_y = 0.5,size = 2.5) +
    geom_tiplab(size=1.5 , aes(color = hbt)) +
    theme(legend.position ="left")+
    ggtitle(title)
  P13_tree[[length(P13_tree)+1]]<-p}
  
  if(scan_ancestor_for_topo(rawtree,tree.df,60,genomes) == 3)
  {P23_fam<-c(P23_fam,title)
  #p<-ggtree(rawtree,branch.length = "none") %<+% LS + 
  p<-ggtree(rawtree) %<+% genomes +
    geom_point(aes(shape = pop,color = hbt), size=2.5) + 
    geom_text2(aes(label = label,subset = !isTip),
               nudge_x = -0.003,nudge_y = 0.5,size = 2.5) +
    geom_tiplab(size=1.5 , aes(color = hbt)) +
    theme(legend.position ="left")+
    ggtitle(title)
  P23_tree[[length(P23_tree)+1]]<-p}
  
  if(scan_ancestor_for_topo(rawtree,tree.df,60,genomes) == 0)
  {
    PNA_fam<-c(PNA_fam,title)
    #p<-ggtree(rawtree,branch.length = "none") %<+% LS + 
    p<-ggtree(rawtree) %<+% genomes +
      geom_point(aes(shape = pop,color = hbt), size=2.5) + 
      geom_text2(aes(label = label,subset = !isTip),
                 nudge_x = -0.003,nudge_y = 0.5,size = 2.5) +
      geom_tiplab(size=1.5 , aes(color = hbt)) +
      theme(legend.position ="left")+
      ggtitle(title)
    PNA_tree[[length(PNA_tree)+1]]<-p}
}

write.table(P12_fam,file = "P12.txt",quote = F,col.names = F,row.names = F)
write.table(P13_fam,file = "P13.txt",quote = F,col.names = F,row.names = F)
write.table(P23_fam,file = "P23.txt",quote = F,col.names = F,row.names = F)
write.table(PNA_fam,file = "PNA.txt",quote = F,col.names = F,row.names = F)

pdf("P12.pdf",width = 15,height = 8,onefile=T)
lapply(P12_tree, function(x) x)
dev.off()

pdf("P13.pdf",width = 15,height = 8,onefile=T)
lapply(P13_tree, function(x) x)
dev.off()

pdf("P23.pdf",width = 15,height = 8,onefile=T)
lapply(P23_tree, function(x) x)
dev.off()

pdf("PNA.pdf",width = 15,height = 8,onefile=T)
lapply(PNA_tree, function(x) x)
dev.off()



save(P12_tree,P13_tree,P23_tree,PNA_tree,file="classified_trees.RData")


