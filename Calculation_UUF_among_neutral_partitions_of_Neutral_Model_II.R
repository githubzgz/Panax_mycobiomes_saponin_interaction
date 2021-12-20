#Figure_6 --> Calculation of unweighted-UniFrac distances among neutral partitions in Neutral Model II

library(picante)
library(dplyr)
library(stringr)
library(phyloseq)

my_rarefy_mean <- function(gin_above,gin_below,gin_neutral,
                           qui_above,qui_below,qui_neutral,
                           not_above,not_below,not_neutral,
                           tree,rand=20){
  #gin_above:A vector of OTU identifier that belonging to the above partition of PG.
  Merge_func <- function(x,y){
    df <- merge(x,y,by='OTU',all.x = T,all.y = T)
    return(df)
  }
  min_num <- min(c(length(gin_above),length(gin_below),length(gin_neutral),
                   length(qui_above),length(qui_below),length(qui_neutral),
                   length(not_above),length(not_below),length(not_neutral)
  ))
  uuf_result <- array(NA,dim = c(9,9,rand))
  for(i in 1:rand){
    print(paste(i,'in',rand,sep=''))
    otu <- Reduce(Merge_func,list(data.frame(OTU=gin_above[sample(length(gin_above),min_num)],gin_above=1),
                                  data.frame(OTU=gin_below[sample(length(gin_below),min_num)],gin_below=1),
                                  data.frame(OTU=gin_neutral[sample(length(gin_neutral),min_num)],gin_neutral=1),
                                  
                                  data.frame(OTU=qui_above[sample(length(qui_above),min_num)],qui_above=1),
                                  data.frame(OTU=qui_below[sample(length(qui_below),min_num)],qui_below=1),
                                  data.frame(OTU=qui_neutral[sample(length(qui_neutral),min_num)],qui_neutral=1),
                                  
                                  data.frame(OTU=not_above[sample(length(not_above),min_num)],not_above=1),
                                  data.frame(OTU=not_below[sample(length(not_below),min_num)],not_below=1),
                                  data.frame(OTU=not_neutral[sample(length(not_neutral),min_num)],not_neutral=1),
    ))
    otu[is.na(otu)] <- 0
    rownames(otu) <- otu$OTU
    otu$OTU <- NULL
    phylo_rarefy <- phyloseq(otu_table(otu,taxa_are_rows = T),phy_tree(f_tree))
    uuf_rarefy <- UniFrac(phylo_rarefy,weighted = F)
    uuf_result[,,i] <- uuf_rarefy %>% as.matrix
  }
  result <- matrix(NA,nrow=9,ncol=9)
  for(i in 1:11){
    for(j in (i+1):9){
      result[i,j] <- mean(uuf_result[i,j,])}}
  result[is.na(result)]<-0
  result <- result + t(result)
  result <- data.frame(result)
  rownames(result) <- c('gin_above','gin_below','gin_neutral',
                        'qui_above','qui_below','qui_neutral',
                        'not_above','not_below','not_neutral')
  names(result) <- rownames(result)
  result <- result %>% as.matrix %>% as.dist
  return(result)
}

rs_uuf_rarefy_mean <- my_rarefy_mean(gin_above = gin_rs_above,gin_below = gin_rs_below,gin_neutral = gin_rs_neutral,
                                     qui_above = qui_rs_above,qui_below = qui_rs_below,qui_neutral = qui_rs_neutral,
                                     not_above = not_rs_above,not_below = not_rs_below,not_neutral = not_rs_neutral,
                                     tree = f_tree,rand = 20)