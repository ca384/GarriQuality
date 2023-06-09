# r read in the dosage file and modify the row names to match with that in the phenotypic file
snps_dosage <- read.csv(here::here("data/snps_199_genotype.csv"), head=T)

dim(snps_dosage)
head(snps_dosage[,1:15])
snps_dosage1 <- (snps_dosage[,-1])
head(snps_dosage1[1:5,1:5])

mt = rrBLUP::A.mat(snps_dosage1)
heatmap(mt)


#install.packages("FactoMineR")
#library(FactoMineR)
PCA(mt)

genotype_snps <- snps_dosage$X # get the genotype names from the snp dosage file to make it look like the names on the phenotype file

write.csv(genotype_snps,file=here("genotype_snps.csv")) # Had genotype names in new order (NR17C2aF10P011)

#r import the new data frame
snps_geno_name <- read.csv(here::here("data/genotype_snps.csv"), head=T)  # had genotype names in old order eg (NR17F10C2aP011)

dim(snps_geno_name)
#snps_dosage_N <- sub('IITA.TMS.','IITA-TMS-',geno_snps$X) # To ma sure all the genotypes names are correct . Rewrite the checks

# combine snps_dosage_N with the dosage file

snps_dosage$X <- snps_geno_name$X # with genotype name in correct order
colnames(snps_dosage)[1] = "geno_name" # rename column 1 of the snps file

rownames(snps_dosage) <- snps_dosage$geno_name # Naming the row names
#removing missing data from the snps file
for (j in 1:ncol(snps_dosage)) {
  snps_dosage[, j] <- ifelse(is.na(snps_dosage[, j]), mean(snps_dosage[, j], na.rm = TRUE), snps_dosage[,
                                                                                                        j])
}
