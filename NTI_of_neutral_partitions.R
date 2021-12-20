#Figure 4 --> NTI calculation
library(reshape2)
library(picante)
library(iCAMP)

#We use the three partitions in Neutral Model I of the rhizosphere soil (RS) as an example
#f_rs_neutral_result -> The neutral partitions of each OTU in this compartment;
#Column 1: OTU identifer; Column 2:Neutral partitions (above,below,and neutral)
#f_tree:Fungal phylogenetic tree constructed using ghost-tree


#1.Integrating community and neutral partition data
f_rs_total_weight <- merge(neutral_taxa$f_rs_neutral_result,f_norm[,names(f_norm)%in%str_subset(names(f_norm),'_RS')],
                           by.x='OTU',by.y='row.names')

#2.Calculating the mean relative abundance of each OTU in this compartment
f_rs_total_weight$mean <- rowMeans(f_rs_total_weight[,3:83]) # A total of 81 samples
f_rs_total_weight[,3:83] <- NULL
f_rs_total_weight <- dcast(f_rs_total_weight,OTU~group)
f_rs_total_weight[is.na(f_rs_total_weight)] <- 0
rownames(f_rs_total_weight) <- f_rs_total_weight$OTU
f_rs_total_weight$OTU <- NULL

#3.Integrating community data and phylogenetic data
f_rs_total_weight_match <- match.phylo.comm(f_tree,f_rs_total_weight %>% t)

#4.Calculating NTI
f_rs_nti <- NTI.p(f_rs_total_weight_match$comm,cophenetic(f_rs_total_weight_match$phy),nworker = 7)