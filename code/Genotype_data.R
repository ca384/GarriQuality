packages_used <- c("workflowr", "tidyverse", "ggplot2", "here", "stringr" ,
                   "sommer","lme4","metan", "rrBLUP","FactoMineR","corrplot", "ASRgenomics", "adegenet")
ip <- installed.packages()
all_packages_installed <- TRUE
for (package in packages_used){
  if (!(package %in% ip[,"Package"])){
    print(paste("Please install package", package))
    all_packages_installed <- FALSE
  } else {
    library(package, character.only = T) # load required package
  } # end else statement
}#END packages_used
if (!all_packages_installed) stop("Need to install packages")



# r read in the dosage file and modify the row names to match with that in the phenotypic file
snps_dosage <- read.csv(here::here("data/snps_199_genotype.csv"), head=T)


snps_dosage1 <- (snps_dosage[,-1])
snps_dosage1 <- as.data.frame(snps_dosage1)
snps_dosage1<- as.matrix(snps_dosage1)

mt = rrBLUP::A.mat(snps_dosage1)
heatmap(mt)


#pc = PCA(mt)
#ds = dist(mt)
#hc = hclust(ds)
#cl = cutree(tree = hc, k = 4)
#cl = as.factor(cl)

#factoextra::fviz_pca_ind(pc, col.ind = "blue", addEllipses = T, habillage = cl)


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
  snps_dosage[, j] <- ifelse(is.na(snps_dosage[, j]), mean(snps_dosage[, j], na.rm = TRUE), snps_dosage[,j])
}


snps_dosage1 <- snps_dosage[,-1]
snps_dosage1 <- as.data.frame(snps_dosage1)
snps_dosage1<- snps_dosage1[,-1]
snps_dosage1<- as.matrix(snps_dosage1)




snpsFilter  <- qc.filtering(M=snps_dosage1, base=FALSE, ref=NULL,
                            marker.callrate=0.2,ind.callrate=0.2, maf=0.05,heterozygosity=0.95,
                            Fis=1, impute=FALSE,  plots=FALSE,message=FALSE)$M.clean




# compute genomic relationship matrix from filtered SNPS marker
Gmat  <- G.matrix(M = snpsFilter, method = "VanRaden", na.string = NA,message=FALSE)$G

heatmap(Gmat)
# Tuning the G-matrix (Gmat) by bending i.e. adjust it near a positive definite
# by making some of its very small or negative eigenvalues slightly positive (Matrix::nearPD # Bend the Gmat matr
G_bend <- G.tuneup(G=Gmat, bend=TRUE, eig.tol = 1e-06,message=FALSE)$Gb

# Tuning the G-matrix (Gmat) by blending #
G_blend <- G.tuneup(G=G_bend, blend=TRUE, pblend=0.02,message=FALSE)$Gb






##### for garri data
snp_numberd = snpsFilter
rownames(snp_numberd) <- 1:198
pc = PCA(snp_numberd,scale.unit = T)
snpsc_new = scale(snp_numberd, scale = T)
ds = dist(snpsc_new)
hc = hclust(ds)
plot(hc, hang = -1)
cl = cutree(tree = hc, k = 5)
cl = as.factor(cl)
factoextra::fviz_pca_ind(pc,
                         col.ind = "blue", addEllipses = F, habillage = cl,
                         geom = "point")

cl = cutree(tree = hc, k = 5)
cl = as.factor(cl)
factoextra::fviz_pca_ind(pc,
                         col.ind = "blue", addEllipses = F, habillage = cl,
                         geom = "point")

#adegenet::find.clusters()
#data.frame(cl_n)

#factoextra::fviz_pca_ind(pc_n,
                         #col.ind = "blue", addEllipses = T, habillage = cl_n)# biplot of the snps data





## here, n.clust is specified, so that only on K value is used
#grp <- find.clusters(snpsFilter,
                    # max.n=30, n.pca=200,
                     #scale= T,
                     #n.clust=4) # takes about 2 minutes
#names(grp)
#grp$stat
#grp$stat

#x <- glSim(50,4e3, 50, ploidy=2)
#x$gen
#plot(x)
#snpda = as.data.frame(snpsFilter)
#pc1 = dapc(x = snpda, grp = grp$grp)

#scatter(pc1, scree.da=F )
#pc1$ind.coord
#plot(pc1$ind.coord[,1], pc1$ind.coord[,2])

snp_numberd = snpsFilter
rownames(snp_numberd) <- 1:198
pc = PCA(snp_numberd,scale.unit = T)
snpsc_new = scale(snp_numberd, scale = T)
ds = dist(snpsc_new)
hc = hclust(ds)
plot(hc, hang = -1)
cl = cutree(tree = hc, k = 5)
cl = as.data.frame(cl)
cl$cl = as.factor(cl$cl)
summary(cl)
cl = as.factor(cl)
factoextra::fviz_pca_ind(pc,
                         col.ind = "blue", addEllipses = T, habillage = cl)

######
#To see the genotypes that are present in each cluster

k4 = kmeans(x = snpsc_new,centers = 4)
k4$cluster

k5 = kmeans(x = snpsc_new,centers = 5)
k5$cluster

head(ddd_new)
mean_dddNew = ddd_new %>% group_by(Geno) %>%
  summarise_at(.vars = c("AMY",    "BD", "CF",
                         "SC",    "SI",
                         "SP",    "STC",
                         "WAC"),
               .funs = mean, na.rm =T)
unique(ddd_new$Geno)
dim(mean_dddNew)


gen_cls = cbind(rownames(snpsFilter),as.data.frame(k4$cluster))
gencl_snp = cbind(snpsFilter,as.data.frame(k4$cluster))
dim(gencl_snp)
#id = which(cl)
#gencl_phen = cbind(meandddNew_ord,as.data.frame(cl))
#dim(meandddNew_ord)
gencls_order = gen_cls[order(gen_cls$`rownames(snpsFilter)`),]
dim(gencls_order)
meandddNew_ord = mean_dddNew[order(mean_dddNew$Geno),]

id = which(gencls_order$`rownames(snpsFilter)` %in% meandddNew_ord$Geno )
gencls_192 = gencls_order[id,]
dim(gencls_192)
gencl_phen = cbind(gencls_192,meandddNew_ord)
ord_cls_phen= cbind(gencls_192, meandddNew_ord)
head(ord_cls_phen)
colnames(ord_cls_phen)[2] = "Cluster"
mean_clusters =ord_cls_phen %>% group_by(Cluster) %>%
  summarise_at(.vars = c("AMY",  "BD",       "CF",
               "SC",       "SI",       "SP",
               "STC", "WAC"), .funs = mean, na.rm = T)
write.csv(x = mean_clusters, file = "mean_of_the_clusters4.csv")
ord_cls_phen_ordbyclt=  ord_cls_phen[order(ord_cls_phen$Cluster),]
write.csv(x = ord_cls_phen_ordbyclt, file = "genotype_cluster.csv", row.names = T)
gncl = read.csv("genotype_cluster.csv")
head(gen_cls)


gencls_order = gen_cls[order(gen_cls$`rownames(snpsFilter)`),]
dim(gencls_order)
meandddNew_ord = mean_dddNew[order(mean_dddNew$Geno),]

id = which(gencls_order$`rownames(snpsFilter)` %in% meandddNew_ord$Geno )
gencls_192 = gencls_order[id,]
dim(gencls_192)
gencl_phen = cbind(gencls_192,meandddNew_ord)
ord_cls_phen= cbind(gencls_192, meandddNew_ord)
head(ord_cls_phen)
colnames(ord_cls_phen)[2] = "Cluster"
mean_clusters =ord_cls_phen %>% group_by(Cluster) %>%
  summarise_at(.vars = c("AMY",  "BD",       "CF",
               "SC",       "SI",       "SP",
               "STC", "WAC"), .funs = mean, na.rm = T)
write.csv(x = mean_clusters, file = "mean_of_the_clusters4.csv")
ord_cls_phen_ordbyclt=  ord_cls_phen[order(ord_cls_phen$Cluster),]
write.csv(x = ord_cls_phen_ordbyclt, file = "genotype_cluster.csv", row.names = T)
gncl = read.csv("genotype_cluster.csv")
head(gen_cls)




gencls_order = gen_cls[order(gen_cls$`rownames(snpsFilter)`),]
dim(gencls_order)
meandddNew_ord = mean_dddNew[order(mean_dddNew$Geno),]

id = which(gencls_order$`rownames(snpsFilter)` %in% meandddNew_ord$Geno )
gencls_192 = gencls_order[id,]
dim(gencls_192)
gencl_phen = cbind(gencls_192,meandddNew_ord)
ord_cls_phen= cbind(gencls_192, meandddNew_ord)
head(ord_cls_phen)

colnames(ord_cls_phen)[2] = "Cluster"
mean_clusters =ord_cls_phen %>% group_by(Cluster) %>%
  summarise_at(.vars = c("AMY",  "BD",       "CF",
                         "SC",       "SI",       "SP",
                         "STC", "WAC"), .funs = mean, na.rm = T)
write.csv(x = mean_clusters, file = "mean_of_the_clusters4.csv")
ord_cls_phen_ordbyclt=  ord_cls_phen[order(ord_cls_phen$Cluster),]
write.csv(x = ord_cls_phen_ordbyclt, file = "genotype_cluster.csv", row.names = T)
gncl = read.csv("genotype_cluster.csv")
head(gen_cls)










## heritability with cluster
ord_cls_phen
str(ord_cls_phen)
ord_cls_phen$Geno = as.factor(ord_cls_phen$Geno)
ord_cls_phen$Cluster = as.factor(ord_cls_phen$Cluster)
Env = levels(ord_cls_phen$Cluster)
Traits = colnames(ord_cls_phen[-c(1:3)])
hert = tibble()
for(env in Env){
  phcl = ord_cls_phen[ord_cls_phen$Cluster== env,]
  for(Trait in Traits){
    eval(parse(text = paste("modelcl= mmer(",Trait,"~  1,
             random=~ vsr(Geno, Gu=Gmat2),
             rcov=~units,
             data=phcl)")))
    ss = summary(modelcl)
    vc = ss$varcomp

    gv = vc[1,"VarComp"]
    ve = vc[2,"VarComp"]
    h2 = gv/(gv +ve)
    h2= cbind("Cluster" = paste(env), "Trait" = paste(Trait),
                 "Heritability" = round(h2,3))
    hert = rbind(hert,h2)
  }
}

write.csv(x = hert, file = "Heritability_of_clusters.csv", row.names = T)


#PCA plot with snps and phenotype
meandddNew_ord_df = as.data.frame(meandddNew_ord)
rownames(meandddNew_ord_df) = meandddNew_ord_df$Geno
pcph = PCA(X = meandddNew_ord_df[,-1], scale.unit = T)

fviz_pca_biplot(X = pcph,geom.ind = "point" ,
                col.ind = ord_cls_phen$Cluster)



######## for chips data
Data_new

snp_numberd_ch = snps_filter_new
rownames(snp_numberd_ch) <- 1:196
snpsc_ch = scale(snp_numberd_ch, scale = T)
pc_ch = PCA(snp_numberd_ch,scale.unit = T)# pca of chips

ds_ch = dist(snpsc_ch)# distance of snps
hc_ch = hclust(ds_ch)# hierarchy cluster of cassava chips snps data
plot(hc_ch, hang = -1)
cl_ch = cutree(tree = hc_ch, k = 5)
cl_ch = as.factor(cl_ch)
factoextra::fviz_pca_ind(pc_ch,
                         col.ind = "blue", addEllipses = F, habillage = cl_ch)

kmeans()

head(Data_new)
mean_Dat_new = Data_new %>% group_by(Genotype) %>%
  summarise_at(.vars = c("AMY_Ch",     "CF_Ch",
                         "SC_Ch",     "STC_Ch"),
               .funs = mean, na.rm =T) # make the mean of Data_new
unique(Data_new$Genotype)
dim(mean_Dat_new)

k_ch = kmeans(x = snps_filter_new,centers = 4) # k cluster for chips
k_ch$cluster

gen_cls_ch = cbind(rownames(snps_filter_new),as.data.frame(k_ch$cluster))# attach names of genotype to the cluster names
gencl_snp_ch = cbind(snps_filter_new,as.data.frame(k_ch$cluster)) # attach the snp data to the cluster
dim(gencl_snp_ch)
#id = which(cl)
#gencl_phen = cbind(meandddNew_ord,as.data.frame(cl))
#dim(meandddNew_ord)
gencls_ch_order = gen_cls_ch[order(gen_cls_ch$`rownames(snps_filter_new)`),]# arrange the cluster alphabetically to make attachment to genotype or phenotype correct
dim(gencls_ch_order)
mean_Dat_new_ord = mean_Dat_new[order(mean_Dat_new$Genotype),] # arrange in alphabetic order

id_ch = which(gencls_ch_order$`rownames(snps_filter_new)` %in% mean_Dat_new$Genotype ) # phenotype that have snps
gencls_ch = gencls_ch_order[id_ch,] # Order
dim(gencls_ch)
gencl_phen_ch = cbind(gencls_ch,mean_Dat_new_ord)
ord_cls_phen_ch= cbind(gencls_ch,mean_Dat_new_ord)# phenotype with cluster number
head(ord_cls_phen_ch)
colnames(ord_cls_phen_ch)[2] = "Cluster"
mean_clusters_ch =ord_cls_phen_ch %>% group_by(Cluster) %>%
  summarise_at(.vars = c("AMY_Ch","CF_Ch", "SC_Ch","STC_Ch"), .funs = mean, na.rm = T)
write.csv(x = mean_clusters_ch, file = "mean_of_the_clusters_Chips.csv")
ord_cls_phen_ordbyclt_ch= ord_cls_phen_ch[order(ord_cls_phen_ch$Cluster),] # arrange ord_cls_phen_ch by cluster
write.csv(x = ord_cls_phen_ordbyclt_ch, file = "genotype_cluster_chips.csv", row.names = T)
gncl_ch = read.csv("genotype_cluster_chips.csv")
head(gncl_ch)


## heritability with cluster for chips
str(ord_cls_phen_ch)
ord_cls_phen_ch$Cluster = as.factor(ord_cls_phen_ch$Cluster)
Env = levels(ord_cls_phen_ch$Cluster)
Traits = colnames(ord_cls_phen_ch[-c(1:3)])
hert_ch = tibble()
for(env in Env){
  phcl_ch = ord_cls_phen_ch[ord_cls_phen_ch$Cluster== env,]
  for(Trait in Traits){
    eval(parse(text = paste("modelcl_ch= mmer(",Trait,"~  1,
             random=~ vsr(Genotype, Gu=Gmat3),
             rcov=~units,
             data=phcl_ch)")))
    ss = summary(modelcl_ch)
    vc = ss$varcomp

    gv = vc[1,"VarComp"]
    ve = vc[2,"VarComp"]
    h2 = gv/(gv +ve)
    h2= cbind("Cluster" = paste(env), "Trait" = paste(Trait),
              "Heritability" = round(h2,3))
    hert_ch = rbind(hert_ch,h2)
  }
}
hert_ch
write.csv(x = hert_ch, file = "Heritability_of_clusters_ch.csv", row.names = T)


#PCA plot with snps and phenotype

mean_Dat_new_ord_df = as.data.frame(mean_Dat_new_ord)
rownames(mean_Dat_new_ord_df) = mean_Dat_new_ord_df$Genotype
pcph_ch = PCA(X = mean_Dat_new_ord_df[,-1], scale.unit = T)

fviz_pca_biplot(X = pcph_ch,geom.ind = "point" ,
                col.ind = ord_cls_phen_ch$Cluster)


