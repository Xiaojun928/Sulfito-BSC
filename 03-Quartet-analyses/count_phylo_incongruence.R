#!/usr/bin/R
library('dplyr')
library(stringr)

idmap <- read.table("gnm_clade_id.txt",sep="\t",stringsAsFactors = F,row.names = 2)
colnames(idmap) <- c('strain','clade')
quartet <- read.table('quartet_cnt_QuartetScores_btwn_clades.txt',sep="\t",stringsAsFactors = F)
colnames(quartet) <- c('comb','topo1','topo2','topo3')

#clades <- c('C1','C2','C3','C6','C7')
clades <- unique(idmap$clade)
target_quartet_comb <- combn(clades,4)

cnt <- matrix(rep(0,3*ncol(target_quartet_comb)),ncol=3,nrow = ncol(target_quartet_comb))
colnames(cnt) <- c('AB|CD','AC|BD','AD|BC')
rownames(cnt) <- sapply(1:ncol(target_quartet_comb), function(i) paste(target_quartet_comb[,i], collapse = ","))
determine_index <- function(ordercomb,refcomb){
  #the resulted quartets are not ordered by C1,C2,C3..Cn,
  #So the order of alternative topologies "AB|CD","AC|BD","AD|BC" would be different from the prior setting, this function is to link the resulted order to the prior setting
  match_list <- list()
  for(m in 1:ncol(ordercomb)) {
    for(n in 1:3) {
      if(all(sort(ordercomb[,m]) == sort(refcomb[,n]))) {
        match_list[[length(match_list)+1]] <- list(ordercomb_index = m, refcomb_index = n)
      }
    }
  }
  match_df <- do.call(rbind, match_list)
  match_df <- as.data.frame(match_df)
  colnames(match_df) <- c("order_index", "ref_index")
  match_df <- match_df %>%
    mutate_at(vars(order_index, ref_index), ~ as.numeric(unlist(.)))
  match_df_reorder<-  match_df %>% mutate(order_index = case_when(
    order_index > 3 ~ 7 - order_index, TRUE ~ order_index))
  return(match_df_reorder)
}

for (r in 1:nrow(quartet)) 
{
  tmp <- quartet[r,]
  involved_clades <- idmap[str_extract_all(tmp[1,1], "[0-9]+")[[1]],'clade']
  #for (i in 1:125){
  for (i in 1:ncol(target_quartet_comb)) {
  {  isct <- intersect(involved_clades,target_quartet_comb[,i])
    if(length(isct) == 4)
    {
      ordercomb <- combn(target_quartet_comb[,i],2)
      refcomb <- combn(involved_clades,2)
      idx_map <- determine_index(ordercomb,refcomb)
      tmpcomb <- paste(target_quartet_comb[,i], collapse = ",")
      for(tpi in 1:3)
      {
        cnt[tmpcomb,idx_map$order_index[tpi]] <- cnt[tmpcomb,idx_map$order_index[tpi]] + tmp[,(idx_map$ref_index[tpi]+1)]
      }
    }
  }
}

#cnt125 <- cnt
#save(cnt125,file="cnt1-125.RData")
save(cnt,file="quartet_cnt.RData")
