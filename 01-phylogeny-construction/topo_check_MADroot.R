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

scan_ancestor_for_topo<-function(rootedtree,rtdf,genomes)
{
  tag=0
    LAS_P12<-getMRCA(rootedtree,rtdf[which(rtdf$label %in% genomes[which(genomes[,2] %in% c("C1","C2")),1]),"node"])
    LAS_P13<-getMRCA(rootedtree,rtdf[which(rtdf$label %in% genomes[which(genomes[,2] %in% c("C1","C3")),1]),"node"])
    LAS_P23<-getMRCA(rootedtree,rtdf[which(rtdf$label %in% genomes[which(genomes[,2] %in% c("C2","C3")),1]),"node"])
    ###if two groups cluster together, 
    ###MCAS shared by A and outgroup should be the same with that shared by B and outgroup
    if(LAS_P12==LAS_P13 && LAS_P12==LAS_P23 && LAS_P13==LAS_P23)
    {
      tag = 0 
      return(tag)
    }
    if(LAS_P12==LAS_P13)##P1 as out
    {  
      tag=3
      return(tag)
    }
    else if(LAS_P12==LAS_P23) ##P2 as out
    {
      tag=2
      return(tag)
    }
    else if(LAS_P13==LAS_P23) ##P3 as out
    {
      tag=1
      return(tag)
    }
    else
    {
      tag = 0 
      return(tag)
    }
  
}

get_btsp <- function(rawtree,tdf,cutoff,genomes)
{
  a <- 1
  b <- 0
  #get the node id and boostrap values for pairwise clades
  LAS_P12<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2] %in% c("C1","C2")),1]),"node"])
  LAS_P13<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2] %in% c("C1","C3")),1]),"node"])
  LAS_P23<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2] %in% c("C2","C3")),1]),"node"])
  bcan12=as.numeric(tdf[which(tdf$node==LAS_P12),"label"])
  bcan13=as.numeric(tdf[which(tdf$node==LAS_P13),"label"])
  bcan23=as.numeric(tdf[which(tdf$node==LAS_P23),"label"])
  
  btsp <- min(c(bcan12,bcan13,bcan23),na.rm = T)
  
  #get the node id and bootstrap value for each clade
  LAS_P1<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2]=="C1"),1]),"node"])
  ban1=as.numeric(tdf[which(tdf$node==LAS_P1),"label"])
  LAS_P2<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2]=="C2"),1]),"node"])
  ban2=as.numeric(tdf[which(tdf$node==LAS_P2),"label"])
  LAS_P3<-getMRCA(rawtree,tdf[which(tdf$label %in% genomes[which(genomes[,2]=="C3"),1]),"node"])
  ban3=as.numeric(tdf[which(tdf$node==LAS_P3),"label"])
  
  btsp_single_clade <- min(c(ban1,ban2,ban3),na.rm = T)
  
  if(btsp > 60 & btsp_single_clade >60)
  {
    return(a)
  }else{
    return(b)
  }
}
##shorten name
shorten <- function(x)
{
  y<-unlist(strsplit(x,"_"))
  paste(y[1],y[2],y[3],sep="_")
}

setwd("E:/by28strains/MAD_results/")
# read *nwk file and lables of the tree
tree_dir <- "E:/by28strains/MAD_results/"
TREE <- list.files(tree_dir,pattern = "*.tre", full.names = T)
TREE_rooted <- list.files(tree_dir,pattern = "*.rooted.nwk", full.names = T)

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
  LABLE<-sub("tre","labels",TREE[i],perl = T)    # get tip labels 
  genomes <- read.table(LABLE,sep="\t",stringsAsFactors = F)
  
  genomes[,1] <- sapply(genomes[,1],shorten)
  
  colnames(genomes)<-c("gene","pop","hbt")
  
  title = sub(".tre","",basename(TREE[i]))
  rawtree <- read.tree(TREE[i])
  rawtree$tip.label <- sapply(rawtree$tip.label,shorten)
  tree.df<-as.data.frame(fortify(rawtree))
  
  rootedtree <- read.tree(TREE_rooted[i])
  rootedtree$tip.label <- sapply(rootedtree$tip.label,shorten)
  rootedtree.df<-as.data.frame(fortify(rootedtree))
  
  rownames(genomes)<-genomes[,1]
  
  topo.type <- scan_ancestor_for_topo(rootedtree,rootedtree.df,genomes)
  conf <- get_btsp(rawtree,tree.df,60,genomes) 
  
  if(topo.type == 1 & conf == 1)
  {P12_fam<-c(P12_fam,title)
  #p<-ggtree(rawtree,branch.length = "none") %<+% LS + 
  p<-ggtree(rawtree) %<+% genomes +
    geom_point(aes(shape = pop,color = hbt), size=2.5) + 
    geom_point2(aes(label=label,subset=!isTip & as.numeric(label) > 60),size=1.5,shape=20) + ##change the btsp as binary point
  geom_text2(aes(label = label,subset = !isTip),
             nudge_x = -0.003,nudge_y = 0.5,size = 2.5) +
    geom_tiplab(size=1.5 , aes(color = hbt)) +
    theme(legend.position ="left")+
    ggtitle(title)
  P12_tree[[length(P12_tree)+1]]<-p}
  
  if(topo.type == 2 & conf == 1)
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
  
  if(topo.type == 3 & conf == 1)
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
  
  if(conf == 0)
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


