#removing missing data from the snps file
for (j in 1:ncol(snps_dosage)) {
  snps_dosage[, j] <- ifelse(is.na(snps_dosage[, j]), mean(snps_dosage[, j], na.rm = TRUE), snps_dosage[,
                                                                                                        j])
}

table(rownames(snps_dosage) == pyt_functional_1$Genotype)

table(pyt_functional_1$Genotype == rownames(snps_dosage) )

#looking for snps that are not in the phenotype data
# gen = which(rownames(snps_dosage) %in%pyt_functional_1$Genotype ) #
# Subset the pheno data that has snp data
gt = which(pyt_functional_1$Genotype %in% rownames(snps_dosage) ) # identify the row number/genotype name those gentypes from pheno to snps


phen1 = droplevels(pyt_functional_1[gt,]) # subset the pheno data that found in the snp data

head(phen1)
dim(phen1)
dim(pyt_functional_1)
length(rownames(snps_dosage)) # length row name of the snps
length(levels(as.factor(phen1$Genotype))) # the length of the selected pheno data genotype levels

# if the legnths are equal no need subseting the snp data
# but if the length is different we should subset the snp data in relation to the selected pheno genotypes
gsnp = which(rownames(snps_dosage) %in% levels(phen1$Genotype))

snp2 = snps_dosage[rownames(snps_dosage)[gsnp],] # subset the snps based on the genotypes with snp and also phenotype
dim(snp2)
length(levels(phen1$Genotype))
Geno_present <- rownames(snp2) %in% levels(phen1$Genotype)

# A_mtrix
A <- rrBLUP::A.mat(snp2[,-1])
heatmap(A)# heatmap for genotypes present in 194 genotypes

A_pca<-prcomp(A, scale=TRUE)

fviz_eig(A_pca) # scree plot of the PCA

# function to loop all mixed models using the GA matrix for all traits
traitnames<- (c("swelling.index","swelling.power", "WAC", "Bulk.Density"))
library(sommer)
h2c = tibble() # an empty table that will contain the heritability for all the traits
vcomp = tibble() # an empty table that will contains the variance components
names(vcomp) = c("vg","vgl","vgy","vp","ve")
# h2c = as.data.frame(h2c)
# rownames(h2c) = traitnames
# h2c = tibble(h2c)
for(traits in traitnames){
  trs = paste0(traits) # the loop will be able to select  and analyse each trait
  eval(parse(text=paste("model1= mmer(",traits,"~ Year + Location + Year:Location,
             random=~vsr(Genotype, Gu=A)+ Genotype:Year + Genotype:Location,
             rcov=~units,
             data=phen1)")))

  ss = summary(model1)
  ss = data.frame(ss$varcomp)
  vg = ss[1,"VarComp"] # assigning names to the variance component in column 1
  vgy = ss[2,"VarComp"]
  vgl = ss[3,"VarComp"]
  ve = ss[4,"VarComp"]
  vp =  vg + vgl/2 + vgy/3  + ve/6
  tt = cbind(vg,vgl, vgy, vp, ve)

  rownames(tt) = paste(traits) # This will calculate the heritability of each of the traits

  vcomp = rbind(vcomp,tt) # this bind the variance components to the functional properties

  name = paste0(traits)
  h2 = vpredict(model1, ~ V1/(V1 + (V2/3) + (V3/2)+(V4/6)) )
  h2 = as.data.frame(h2)
  rownames(h2) = paste0(traits)
  h2c = rbind(h2c,h2)

}
vcomp
h2c

df_mmer <- cbind(vcomp,h2c) # binding the mmer variance components with heritability
df_mmer <- data.frame(df_mmer) #converting to a data frame
colnames(df_mmer)[6] <- "h2" # rename Estimates to h2
colnames(df_mmer)[7] <- "SE h2"
View(df_mmer)
df_mmer

# creating tables
df_lmer <- cbind(Result.vcomp, T_h2 ) # bind the two data frames

colnames(df_mmer)[6] = "H2Comb"
df_mmer1 = df_mmer[,-7]
rownames(df_mmer1) <- c("swelling.index_mmer", "swelling.power_mmer", "WAC_mmer", "Bulk.Density_mmer")
rownames(df_lmer) <- c("swelling.index_lmer", "swelling.power_lmer", "WAC_lmer", "Bulk.Density_lmer")

lmer_mmer <- rbind(df_lmer, df_mmer1) # bind the two data frames

lmer_mmer

lmer_mmer_sort <- lmer_mmer[sort(rownames(lmer_mmer)),] # sort the lmer_mmer data frame according to the rownames
lmer_mmer_sort

