# Figure 1 - The variation patterns of Panax mycobiomes

# f_com --> Fungal community profile (row is sample, column is OTU)
# metadata --> Factor table (compartment->plant compartments;plant->plant species;site->sampling site groups;year->plant growth year)

library(vegan)

#1. ¦Â-diversity of Panax mycobiomes
f_bray <- vegdist(f_com)

#2. PERMANOVA
adonis(f_bray ~ compartment + plant + site + year,data = metadata)

#3. PPER for each compartment
#For each compartment, the combined explanatory effects of plant species, site groups, and growth year consists of seven parts:
#a -> Pure effect of plant species;
#b -> Pure effect of site group;
#c -> Combined effect of plant species + site group;

#Now we use the Bray-curtis dissimilarity of mycobiome in rhisozphere soil (RS) as an example (rs_bray);
#Corresponding metadata is rs_metadata.

#Calculating the part a.
test_1 <- adonis(rs_bray ~ plant,data = rs_metadata)$R2[1] # a + c
test_2 <- adonis(rs_bray ~ site,data = rs_metadata)$R2[1] # b + c
test_3 <- 1 - adonis(rs_bray ~ plant + site,data = rs_metadata)$R2[3] # a + b + c

a <- test_3 - test_2
b <- test_3 - test_1

PPER = a/b

