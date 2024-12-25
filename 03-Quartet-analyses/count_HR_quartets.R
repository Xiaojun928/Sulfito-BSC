#!/usr/bin/R
library(parallel)
library(stringr)
library(reshape2)
library(dplyr)
library(ggplot2)

idmap <- read.table("gnm_clade_id.txt",sep="\t",stringsAsFactors = F,row.names = 2)
colnames(idmap) <- c('strain','clade')
clades <- unique(idmap$clade)

## 0. list all possible quartets for homologous recombination detection
target_qid <- c()
multistrain_clade <- c('C1','C2','C3','C6','C7')
singleton_clade <- setdiff(clades,multistrain_clade)
sample_w_replacement <- expand.grid(clades,clades,clades,clades)
factor2char <- function(x)
{as.character(unlist(x))}
sample_w_replacement <- sample_w_replacement %>% mutate_if(is.factor,factor2char)
for (n in 1:nrow(sample_w_replacement)) {
  tmp_clades <- as.character(sample_w_replacement[n,])
  uni <- unique(tmp_clades)
  if(length(uni) == 3) #e.g., two strain from C1, one from C2 and the other from C3
  {
    if(names(which.max(table(tmp_clades))) %in% multistrain_clade)
      target_qid <- c(target_qid,n)
  }
  if(length(uni) == 2) 
  {
    cnt_singleton <- length(intersect(uni,singleton_clade))
    if(min(table(tmp_clades)) == 2 & cnt_singleton == 0){ ##e.g. two strain from C1, another two from C2
      target_qid <- c(target_qid,n)
    }
    if(min(table(tmp_clades)) == 2 & cnt_singleton > 0) #e.g., two strain from C1, at least one from singleton clades
	{
		target_qid <- c(target_qid,n)
	}
  }
}
target_quartet <- sample_w_replacement[target_qid,c(4,3,2,1)]

## 1. find quartets with strains from 3 or 2 clades
if (! file.exists("QuartetScores_quartet4HR.RData")) {
  quartet <- read.table('quartet_cnt_QuartetScores.txt',sep="\t",stringsAsFactors = F)
  colnames(quartet) <- c('comb','topo1','topo2','topo3')
  
  process_row <- function(r) {
    tmp <- quartet[r,]
    if(length(unique(unlist(str_extract_all(tmp[1,1], "[0-9]+")))) == 4){
      involved_clades <- idmap[str_extract_all(tmp[1,1], "[0-9]+")[[1]],'clade']
      uni <- unique(involved_clades)
      if(length(uni) == 3 | (length(uni) == 2 & min(table(involved_clades)) == 2)) {
        return(r)
      } else {
        return(NA)
      }
    }
  }
  
  num_cores <- 32
  qid_parallel <- mclapply(1:nrow(quartet), process_row, mc.cores = num_cores)
  
  qid_parallel <- unlist(qid_parallel)
  qid_parallel <- qid_parallel[!is.na(qid_parallel)]
  
  quartet_w_same_clds <- quartet[qid_parallel,]
  save(quartet_w_same_clds,file="QuartetScores_quartet4HR.RData")
}else{

  load("QuartetScores_quartet4HR.RData")
  # 2. count the quartet supporting/against recent recombination
  cnt4HR <- matrix(rep(0,2*nrow(target_quartet)),nrow=nrow(target_quartet),ncol=2)
  colnames(cnt4HR) <- c('against_hgt','support_hgt')
  rownames(cnt4HR) <- apply(target_quartet,1,function(row) paste(row,collapse = ", "))
  quartets_record <- list()
  
 for (r in 1:nrow(quartet_w_same_clds)) 
  {
    tmp <- quartet_w_same_clds[r,]
    involved_clades <- idmap[str_extract_all(tmp[1,1], "[0-9]+")[[1]],'clade']
    uni <- unique(involved_clades)
    refcomb <- combn(involved_clades,2)
    if(length(uni) == 3)
    {
      for (idx in 1:ncol(refcomb)) {
        if(refcomb[,idx][1] == refcomb[,idx][2])
        {
          if(idx >3){idx <- 7-idx}
          rowid <- paste(sort(involved_clades),collapse = ", ")
          sp <- sum(tmp[1,2:4]) - tmp[1,c(idx+1)]
          ag <- tmp[1,c(idx+1)]
          cnt4HR[rowid,'against_hgt'] <- ag + cnt4HR[rowid,'against_hgt']
          cnt4HR[rowid,'support_hgt'] <- sp + cnt4HR[rowid,'support_hgt']
          quartets_record[[rowid]][length(quartets_record[[rowid]]) + 1] <- paste(idmap[str_extract_all(tmp[1,1], "[0-9]+")[[1]],'strain'],collapse=", ")
          break
        }
      }
    }
    if(length(uni) == 2)
    {
      for (idx in 1:ncol(refcomb)) {
        if(refcomb[,idx][1] == refcomb[,idx][2])
        {
          if(idx >3){idx <- 7-idx}
          rowid <- paste(sort(involved_clades),collapse = ", ")
          sp <- sum(tmp[1,2:4]) - tmp[1,c(idx+1)]
          ag <- tmp[1,c(idx+1)]
          cnt4HR[rowid,'against_hgt'] <- ag + cnt4HR[rowid,'against_hgt']
          cnt4HR[rowid,'support_hgt'] <- sp + cnt4HR[rowid,'support_hgt']
          quartets_record[[rowid]][length(quartets_record[[rowid]]) + 1] <- paste(idmap[str_extract_all(tmp[1,1], "[0-9]+")[[1]],'strain'],collapse=", ")
          break
        }
      }
    }    
  } 
}

list_result <- list(cnt4HR,quartets_record)
save(list_result,file="cnt4HR.RData")

cnt4HR <- as.data.frame(cnt4HR)
cnt4HR <- cnt4HR %>% filter(against_hgt > 0 & support_hgt > 0)
prop <- cnt4HR %>% mutate_all(~ . / rowSums(cnt4HR) * 100) %>% mutate(Quartet = rownames(cnt4HR))
long_cnt <- melt(prop,id.vars='Quartet',variable.name = "topo",value.name = "Proportion")
p <- ggplot(long_cnt, aes(x = Quartet, y = Proportion, fill=topo)) +
  scale_fill_manual(name="Topology", values = c("against_hgt" = "gray", "support_hgt" = "#fb8702")) +
  geom_bar(stat = "identity") +  
  theme(axis.text.x = element_blank())

pdf(file="quartet_prop4hgt_w_singleton.pdf",height=8,width=12)
p
dev.off()



