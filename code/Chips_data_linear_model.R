packages_used <- c("workflowr", "tidyverse", "ggplot2", "here", "stringr" , "sommer","lme4","metan",
                   "rrBLUP","FactoMineR","corrplot", "ASRgenomics", "reshape2","pheatmap", "Hmisc", "agricolae")

#install.packages(c("agricolae"))
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


# read in the .csv for all chips and years and locationchips
oto_2021_data <- read.csv(here::here("data/PYT_chemical_data/chips/All_Otobi_2020_2021_chips_data.csv"))

umu_2021_data <- read.csv(here::here("data/PYT_chemical_data/chips/All_umudikie_2020_2021_chips_data.csv"))

oto_2022_data <- read.csv(here::here("data/PYT_chemical_data/chips/All_Otobi_2021_2022_chips_data.csv"))

umu_2022_data <- read.csv(here::here("data/PYT_chemical_data/chips/All_umudikie_2021_2022_chips_data.csv"))

All_chips_data <- rbind(oto_2021_data, umu_2021_data, oto_2022_data,umu_2022_data)


dim(All_chips_data)
head(All_chips_data)
str(All_chips_data)
All_chips_data$Genotype <- as.factor(All_chips_data$Genotype)
All_chips_data$Year <- as.factor(All_chips_data$Year)
All_chips_data$Location <- as.factor(All_chips_data$Location)
dim(All_chips_data)
str(All_chips_data)

chips_all<- All_chips_data %>%
  dplyr::rename("AMY_Ch" ="amylose_chips","CF_Ch" = "cf_chips", "STC_Ch"="starch_chips", "SC_Ch"= "sugar_chips")

source("/Users/ca384/Documents/ChinedoziRepo/GarriQuality/code/Genotype_data.R")
XXX <- which(chips_all$Genotype %in% rownames(snpsFilter)) # All chips present in snpsfilter
YYY <- which(rownames(snpsFilter) %in% chips_all$Genotype) # All snps present in alls chips dataframe
chips_all =droplevels(chips_all[XXX, ]) # Select the XXX ie phenotypes that have genotype data i
head(snpsFilter[,1:10])
snps_filter_new = snpsFilter[YYY,] # new snps data frame containing phenotype
all(which(chips_all$Genotype %in% rownames(snps_filter_new))) # logical to confirm chips rownames in snps rownames
all(which( rownames(snps_filter_new) %in%chips_all$Genotype)) # logical statement to confirm snps_filter_new in chips_new$Genotype line 27 and 28 should give true
colnames(chips_all)
dd_no_snp = chips_all[-XXX,] # data frame without the snps data
dim(dd_no_snp)
unique(dd_no_snp$Genotype)


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


#yyy <- which(rownames(snpsFilter) %in% All_chips_data$Genotype)
#unique(All_chips_data$Genotype[XXX]) %in% unique(rownames(snpsFilter)[yyy])
length(unique(chips_all$Genotype))

head(chips_all)
# add dmc
dmc_oto_21 <- read.csv("~/Documents/ChinedoziRepo/GarriQuality/data/ALL_PYT_DMC/oto_2021_DMC.csv")
dmc_oto_21$Year = "2021"
dmc_oto_21$Location = "Otobi"
dmc_umu_21 <- read.csv("~/Documents/ChinedoziRepo/GarriQuality/data/ALL_PYT_DMC/umu_2021_DMC.csv")
dmc_umu_21$Year = "2021"
dmc_umu_21$Location ="Umudike"
dmc_oto_22 <- read.csv("~/Documents/ChinedoziRepo/GarriQuality/data/ALL_PYT_DMC/Oto_dmc_2022.csv")
dmc_oto_22$Year = "2022"
dmc_oto_22$Location = "Otobi"
dmc_umu_22 <- read.csv("~/Documents/ChinedoziRepo/GarriQuality/data/ALL_PYT_DMC/Umu_DMC_2022.csv")
dmc_umu_22$Year = "2022"
dmc_umu_22$Location = "Umudike"
ALL_DMC <- rbind(dmc_oto_21, dmc_umu_21, dmc_oto_22, dmc_umu_22)

ALL_DMC<- ALL_DMC %>%
  dplyr::rename("Genotype"= "genotype")
All_dmc1 = ALL_DMC %>% group_by(Location, Year, Genotype) %>%
  summarise_at(.vars = "DMC", .funs = mean)
str(All_dmc1)
chips_all<-chips_all%>%
  full_join(All_dmc1, by = c("Location", "Year", "Genotype") )

chips_all = chips_all[,c("Location", "Year", "Genotype", "AMY_Ch", "CF_Ch", "STC_Ch", "SC_Ch", "DMC")]

# convert each of the samples to dry weight basis

chips_all$Starch_Dm = (chips_all$STC_Ch*(chips_all$DMC/100))
chips_all$AMY_Dm = (chips_all$AMY_Ch*(chips_all$DMC/100))
chips_all$CF_Dm = (chips_all$CF_Ch*(chips_all$DMC/100))
chips_all$SC_Dm = (chips_all$SC_Ch*(chips_all$DMC/100))

 chips_all<- chips_all %>%
 dplyr::rename("STC_Dm" = "Starch_Dm")

chips_all2<- chips_all %>%
  dplyr::rename("STC_ch" = "STC_Dm","AMY_ch" ="AMY_Dm", "CF_ch" = "CF_Dm", "SC_ch"= "SC_Dm")
str(chips_all2)

chips_all2$Genotype = as.factor(chips_all2$Genotype)
chips_all2$Year = as.factor(chips_all2$Year)
chips_all2$Location = as.factor(chips_all2$Location)
# STC_ch is starch contented converted to fresh weight while STC_Ch is on dry basis
#chips_new$CF_T= (chips_new$CF_Ch)^(1/3)
# compute genomic relationship matrix from filtered SNPS marker
Gmat2  <- G.matrix(M = snps_filter_new, method = "VanRaden", na.string = NA,message=FALSE)$G
Gmat3 = A.mat(snps_filter_new)
heatmap(Gmat3)
id1 = which(chips_all2$Genotype %in% rownames(Gmat3))

chips_all2 = droplevels(chips_all2[id1,])
dim(chips_all2)
chips_all2$Genotype %in% rownames(Gmat3)
rownames(Gmat3) %in% chips_all2$Genotype
# Tuning the G-matrix (Gmat) by bending i.e. adjust it near a positive definite
# by making some of its very small or negative eigenvalues slightly positive (Matrix::nearPD # Bend the Gmat matr
G_bend2 <- G.tuneup(G=Gmat2, bend=TRUE, eig.tol = 1e-06,message=FALSE)$Gb

# Tuning the G-matrix (Gmat) by blending #
G_blend2 <- G.tuneup(G=G_bend2, blend=TRUE, pblend=0.02,message=FALSE)$Gb


colnames(chips_all2)
traitnames<- (c( "DMC", "STC_ch", "AMY_ch", "CF_ch", "SC_ch") )
vcomp = tibble() # an empty table that will contains the variance components
filtered_value_new <- tibble() # for boxplot

  for(traits in traitnames){
    trs_chips = paste0(traits) # the loop will be able to select  and analyse each trait
    out_ind_ch <- which(chips_all2[,paste(traits)] %in% boxplot.stats(chips_all2[,paste(traits)])$out) # looking for outliers

    if(length(out_ind_ch) == 0){
      Data_new = chips_all2
    }else{
      Data_new = chips_all2[-out_ind_ch,]
    }
    eval(parse(text=paste(
      "modelA= mmer(",traits,"~ 1,
  random= ~Year  + Location + Year:Location + vsr(Genotype, Gu=Gmat2) +  Genotype:Year + Genotype: Location,
     rcov=~units, tolParInv = 0.1,
  data=Data_new)")))
    summary(modelA)
    ss = summary(modelA)
    ss = data.frame(ss$varcomp)
    vg = ss[4,"VarComp"] # assigning names to the variance component in column 1
    vgy = ss[5,"VarComp"]
    vgl = ss[6,"VarComp"]
    ve = ss[7,"VarComp"]
    vp =  vg + vgy/2 + vgl/2 + ve/4
    h = round(vg/vp,3)
    tt = cbind("Trait" = paste(traits), "Method" = "mmer", "vg" = round(vg,3), "vgy" = round(vgy,3), "vgl" = round(vgl,3), "vph" = round(vp,3), "ve" = round(ve,3), "H2" = h)

    rownames(tt) = paste(traits) # This will calculate the heritability of each of the traits

    ## variance comp and h2 for lmer

    eval(parse(text=paste(
      "modelB= lmer(",traits,"~ (1|Year) + (1|Location) + (1|Year:Location)
      + (1|Genotype) + (1|Genotype:Year)+ (1|Genotype:Location),
  data=Data_new)")))


    summary(modelB)
    ss1 = summary(modelB)
    ss1 = data.frame(ss1$varcor)
    vg1 = ss1[ss1$grp == "Genotype","vcov"] # assigning names to the variance component in column 1
    vgy1 = ss1[ss1$grp == "Genotype:Year","vcov"]
     vgl1 = ss1[ss1$grp == "Genotype:Location","vcov"]
    ve1 = ss1[ss1$grp == "Residual","vcov"]
    vp1 =  vg1  + (vgy1/2)+ (vgl1/2) + (ve1/4)
    h1 = round(vg1/vp1,3)

    tt1 = cbind("Trait" = paste(traits), "Method" ="lmer","vg" = round(vg1,3), "vgy" = round(vgy1,3), "vgl"= round(vgl1,3),"vph" = round(vp1,3), "ve" = round(ve1,3), H2 = h1)

    vccomb = rbind(tt,tt1)
    vcomp = rbind(vcomp,vccomb) # this bind the variance components to the functional properties

    Bx_A <-Data_new[,traits]
    Bx_B <- cbind(Year=as.character(Data_new$Year), Location=as.character(Data_new$Location),Geno = as.character(Data_new$Genotype), trait=paste0(traits), value=Bx_A)
    filtered_value_new= rbind(filtered_value_new, Bx_B)
  }

  # long format
 vcomp_chips <- data.frame(vcomp)

 write.csv(x =vcomp_chips, file = "~/Documents/ChinedoziRepo/GarriQuality/output/Heritabity of chips_fresh_weight_with_GL.csv")

#filtered_value_new$value = as.numeric(filtered_value_new$value)
 # filtered_value_new$value=as.numeric(filtered_value_new$value)
 # ggplot(filtered_value_new, aes(x = Geno, y=value, fill=Location))+
  #  geom_boxplot() +
  #  facet_grid(trait~Year,  scales = "free")


  filtered_value_new$value=as.numeric(filtered_value_new$value)
  filtered_value_new$Year= as.factor( filtered_value_new$Year)
  filtered_value_new$Geno= as.factor( filtered_value_new$Geno)
  filtered_value_new$Location= as.factor( filtered_value_new$Location)
  filtered_value_new$trait= as.factor( filtered_value_new$trait)

  filtered_value_new2<- filtered_value_new[(filtered_value_new$trait %in% c("AMY_ch", "CF_ch", "STC_ch", "SC_ch", "DMC")),]
unique(filtered_value_new2$trait)
  global_size=8 # font size of the axises
  ggplot(filtered_value_new2, aes(x = Location, y=value, fill=Location))+
    geom_boxplot(outlier.shape = NA) + facet_grid(trait~Year,  scales = "free") +
     theme_classic(base_size = global_size) # gives the same fount size to the axes
    # geom_point (data = aggregate(traits ~ Genotype, data =Data_new , mean, use=complete.obs)) +
    # facet_grid(trait~Year,  scales = "free")
  # boxplot

  head(filtered_value_new2)
datfilterd_new_final = dcast(data = filtered_value_new2, formula = Location + Year + Geno ~ trait,
      fun.aggregate = mean, value.var = "value")
str(datfilterd_new_final)
datfilterd_new_final$Location  = as.factor(datfilterd_new_final$Location)
datfilterd_new_final$Year = as.factor(datfilterd_new_final$Year)
datfilterd_new_final$Geno = as.factor(datfilterd_new_final$Geno)
  write.csv(x =vcomp, file = "~/Documents/ChinedoziRepo/GarriQuality/output/Heritab_lmer_mmer_chips_2yrs.csv")
library(reshape2)
  filtered_value_new$value = as.numeric(filtered_value_new$value)
  ddd_chips = dcast(data = filtered_value_new2, formula = Year + Location + Geno ~ trait, fun.aggregate = mean, value.var = "value", na.rm= T) # mean of the genotypes across location and years
  indv_cor_chips <- round(cor(ddd_chips[, -c(1,2,3)],use = "pairwise.complete.obs"),3)
  corrplot( indv_cor_chips, type = "upper",method = "number", number.digits = 3,
            is.corr = T, number.cex =0.5,sig.level=T  ) # correlation plot


######################################################
# Genomic selection of chips
  library(emmeans)
  library(lme4)

  set.seed(123)

  traits<- (c("AMY_ch","CF_ch","STC_ch","SC_ch", "DMC"))

  mean_Data2<- tibble()
  for(Trait in traits){
    eval(parse(text = paste("ansm1 <- lmer(",Trait,"~ (Genotype) + Location +
    (1|Year) ,
      data= Data_new)")))
    em2 = emmeans(object = ansm1, specs = ~ Genotype)
    em2 = as.data.frame(em2)
    m2 = cbind(Trait = paste(Trait), em2[,1:2])
    mean_Data2 = rbind(mean_Data2,m2)

  }

  library(reshape)
  Meanwide3 <-dcast(data = mean_Data2, formula = Genotype ~ Trait, fun.aggregate = mean,
                    value.var = "emmean")
  #mean_Data1 <- Data1%>% group_by(Genotype) %>%
  #summarise_if(is.numeric, mean, na.rm = TRUE)


  Acc_pred_n <- tibble()

  idmark = which(rownames(Gmat3) %in% Meanwide3$Genotype)
  Gmat3 = Gmat3[idmark, idmark]
  idphn = which(Meanwide3$Genotype %in% rownames(Gmat3))
  Meanwide3 = Meanwide3[idphn,]
  all(unique(Meanwide3$Genotype) %in% rownames(Gmat3))

  for(Trait in traits) {
    # Creating a folder that contain 5 subset with 5times with a total of 5*5*8(traits)= 200
    fold5 = caret::createMultiFolds(y = unique(Meanwide3$Genotype), k = 5, times = 5)


    for(i in 1:length(fold5)){
      index = fold5[[i]] # the index of the sample for training set
      #subset the phenotypic data
      train_geno_n <- droplevels(unique(Meanwide3$Genotype)[index])
      train_geno_ind_n <- which(Meanwide3$Genotype %in% train_geno_n)
      train.data_n <- droplevels(Meanwide3 %>%
                                 filter(row_number() %in% train_geno_ind_n)) # subset the training set

      test.data_n <- droplevels(Meanwide3 %>%
                                filter(!row_number() %in% train_geno_ind_n)) # subset the testing set

      test.data2 <- test.data_n
      test.data2[,traits] = NA # change the the test data of the training set to NA value

      mod_dat3<- rbind(train.data_n, test.data2) # combine the the data set for analysis

      #####################
        # str(datfilterd_new_final)
      # head(datfilterd_new_final)
      # colnames(datfilterd_new_final)[3] = "Genotype"
      # summary(mod_dat3)
      eval(parse(text = paste("ansm2 <- mmer(",Trait,"~1,
  random=~ vsr(Genotype, Gu=Gmat3),
      rcov=~units,
      data=mod_dat3)")))
ft = fitted(ansm2)
ft$dataWithFitted
ansm2$U$`u:Genotype`
      pblup2 = as.data.frame(ansm2$U$`u:Genotype`) + ansm2$Beta[,"Estimate"]
      pblup2 = cbind(rownames(pblup2),pblup2)
      colnames(pblup2)[1] = "Genotype"
      testind_n = which(pblup2$Genotype %in% unique(test.data2$Genotype))
      testpred_n = pblup2[testind_n,]
      idm = which(unique(test.data2$Genotype) %in% testpred_n$Genotype)
      testpred_n = testpred_n[order(testpred_n$Genotype),]

      obs_test_n =  test.data_n %>% group_by(Genotype) %>%
        summarise_at(.vars = Trait, .funs = mean)

      obs_test_n= obs_test_n[order(obs_test_n$Genotype),]
      obs_test_n = obs_test_n[idm,]
      r_2 = cbind( Trait, predictability = round(cor(testpred_n[,Trait],obs_test_n[,Trait], use = "pairwise.complete.obs"),3))
      colnames(r_2)[2] <- c("predictability")
      Acc_pred_n = rbind(Acc_pred_n,r_2)


    }
  }
  #}
str(Acc_pred_n)
  Acc_pred_n$predictability= as.numeric(Acc_pred_n$predictability)
  ggplot(Acc_pred_n, mapping = aes(x = Trait, y = predictability, fill = Trait)) +
    geom_boxplot() + theme_bw()

#####################################
#combining garri and chips on a single data frame
  #Data7 is the data frame containing the garri data from line 1264 linear-model_garri_quality
  # Data1 contains the chips data

  garri_chips<-Data_new%>%
                     full_join(Data7, by = c("Location", "Year", "Genotype") )

length(Data7$Genotype %in% Data_new$Genotype)
length(Data_new$Genotype %in% Data1$Genotype)
garri_chips2<- garri_chips %>%
  dplyr::select("Location", "Year", "Genotype","SI" , "SP", "WAC", "BD", "AMY", "CF", "STC","SC","AMY_ch","CF_ch","STC_ch","SC_ch", "DMC" )
colnames(garri_chips2)
summary_stat=summary(garri_chips2)


colnames(garri_chips2)
summary_stat_ch_dm=summary(garri_chips2)
write.csv(x =summary_stat_ch_dm, file = "~/Documents/ChinedoziRepo/GarriQuality/output/summary_stat_chips_dm")


# r1 = cor(garri_chips2[,-c(1:3)], use = "pairwise.complete.obs")
# corrplot.mixed(r1, type ="full"  ,method = "number", number.digits = 2,
#          is.corr = T,  cl.cex = 0.5, tl.cex = 0.5, number.font = 0.01,number.cex = 0.6, addshade = c("negative", "positive", "all"))
# r2= correlation(x=garri_chips2[,-c(1:3,13:16)], y=NULL,
#            method = "pearson",
#            alternative ="two.sided")

write.csv(x =r2, file = "~/Documents/ChinedoziRepo/GarriQuality/output/correlation_chips_matrix.csv")

out_cf <- which(garri_chips2$CF %in% boxplot.stats(garri_chips2$CF)$out)
garri_chips2[out_cf, "CF"] = NA # make outlier NA

library(psych)


pairs.panels(garri_chips2[,-c(1:3)],   # plot distributions and correlations for all the data
             gap = 0,
             pch = ".",
             cex = 1.5,
             lm = TRUE,
             ellipses=FALSE,stars = TRUE, method = "pearson")



head(garri_chips2)

Mean_GC2 <- garri_chips2 %>% group_by(Genotype) %>% summarise_at(.vars = c( "SI", "SP", "WAC", "BD","AMY", "CF","STC","SC","DMC", "AMY_ch","CF_ch","STC_ch","SC_ch"), .funs = mean, na.rm=T)
Mean_GC2 = as.data.frame(Mean_GC2 )
rownames(Mean_GC2) <- Mean_GC2$Genotype
Mean_GC2 = Mean_GC2[,-1]
Mean_GC2_sc = scale(Mean_GC2, scale = T)
ds_ch = dist(Mean_GC2_sc)
view(ds_ch)
dim(ds_ch)
ds1 = as.matrix(ds_ch)
# ds2 = ds1[rowSums(is.na(ds1)) == 0, colSums(is.na(ds1)) == 0, drop = F]


which(is.na(ds1))
x = tibble()
for(i in rownames(ds1) ){
ln = length(which(is.na(ds1[i, ])))
ln1 = cbind(rn = i, ln)
x = rbind(x, ln1)
}
x$ln = as.numeric(x$ln)
x[x$ln >0, ]

xc = tibble()
for(i in colnames(ds1) ){
  cn = length(which(is.na(ds1[, i])))
  cn1 = cbind(colN = i, cn)
  xc = rbind(xc, cn1)
}
xc$cn = as.numeric(xc$cn)
ncl3 = xc[xc$cn >10, ]
idc = which(colnames(ds1) %in% ncl3[,1])
ds2 = ds1[!rownames(ds1)%in% as.character(ncl3[,1]),!rownames(ds1)%in% as.character(ncl3[,1])]
rownames(ds2)
hc = hclust(dist(ds2))
gen_list = as.data.frame(hc$labels)
gen_list$code = 1:195 # add numbers to the genotype name list

rownames(ds2) <- gen_list$code
colnames(ds2) <- gen_list$code
hc = hclust(dist(ds2))
plot(hc, hang = -1, cex = 0.4)

install.packages(c("circlize", "dendoextras"))
install.packages("dendextend")
# install.packages("circlize")
library(circlize)
library(dendextend)
hcd = as.dendrogram(hc)
hcd <- hcd %>%
  color_branches(k = 4) %>%
  color_labels(k = 4)
labels_cex(hcd) <- .5
circlize_dendrogram(hcd)
summary(hcd)
is.list(hcd)
tree = cutree(hcd, k = 4)
gen_list$cluster = tree # contain the cluster and the genotypes in the each cluster.

write.csv(gen_list, "/Users/ca384/Documents/ChinedoziRepo/GarriQuality/output/dendrogram_cluster_list.csv")
write.csv(Mean_GC2, "/Users/ca384/Documents/ChinedoziRepo/GarriQuality/output/trait_mean_for_dendro.csv")

#hc <- hc %>%
 # color_branches(k = 3) %>%
  #color_labels(k = 3)

#circlize_dendrogram(hc)

#plot(hcd)
#plot(hc, hang= -1, cex = 0.3)
#rownames(ds1)
#dim(ds2)
#dim(ds1)
#hclust(ds2)
#is.matrix(ds2)
#summary(ds2)
#view(ds2)
#is.na(ds2[,1])
#garri_chips_dmc2<- garri_chips2 %>%
#  dplyr::select("Location", "Year", "Genotype", "DMC","AMY_ch" , "CF_ch", "STC_ch", "SC_ch" )


Mean_GC2 <- garri_chips_dmc2 %>% group_by(Genotype) %>% summarise_at(.vars = c( "DMC", "AMY_ch","CF_ch","STC_ch","SC_ch"), .funs = mean, na.rm=T)
library (reshape2)

wide_mean_data2 <- dcast(data = mean_Data2, formula = Genotype~ Trait,fun.aggregate = mean,value.var = "emmean" ) # long to wide format


#cor_GC <- cor(Mean_GC[,-1], use = "pairwise.complete.obs")
#corrplot(cor_GC, method = "pie",is.corr = F , )
#corrplot(cor_GC, type = "upper",method = "number", number.digits = 3, is.corr = T)

Mean_GC2 = as.data.frame(wide_mean_data2)
rownames(Mean_GC2) <- Mean_GC2$Genotype
Mean_GC2 = Mean_GC2[,-1]
pc = PCA(Mean_GC2)
dim(pc$ind$coord)
str(Mean_GC2)
 Mean_GC2_sc = scale(Mean_GC2, scale = T)
 pc_ch = PCA(Mean_GC2)
 ds_ch = dist(Mean_GC2)

 ds1 = as.matrix(ds_ch)
 #ds2 = ds1[rowSums(is.na(ds1)) == 0, colSums(is.na(ds1)) == 0, drop = F]
 hc = hclust(as.dist(ds1))

 tr1 = cutree(hc, k =5)
 length(tr1)
 tr1 = data.frame(tr1)
 pcind_ch = as.data.frame(pc_ch$ind$coord)
 id1 = which(rownames(pc_ch$ind$coord) %in% rownames(tr1))

 dim(pcind_ch)
 pcind_ch = pcind_ch[id1,]
 pcind_ch$cl = tr1[,1]
 pcind_ch$cl = as.factor(pcind_ch$cl)
 ggplot(data = pcind_ch, aes(x = Dim.1, Dim.2)) +
   geom_point(colour = pcind_ch$cl)  + theme_bw()
tr1$tr1 = as.factor(tr1$tr1)
 fviz_pca_biplot(pc_ch, geom.ind = "point", col.ind = tr1[,1],)
  h11 = pheatmap(ds1, show_rownames = T, show_colnames = F)
 h11$tree_row$labels

 # kmean clustering
#  ?kmeans
# id =  rownames(Mean_GC1) %in% rownames(tr1)
# Mean_GCK = Mean_GC1[id,]
#  kk = kmeans(Mean_GCK, centers = 4)
#  kk$cluster # cluster number
#  kk$centers # the mean value of the clusters

 ## dendorogram
 # Mean192 = Mean_GC[id,]
 # ?scale
 # Mean1
 # Mean192S = scale(x = Mean192[,-1], scale = T)
 # disk = dist(Mean192S)
 #
 # hc192 = hclust(disk, method = "ward.D")
 # hc192$labels
 # plot(hc192, hang = -1, labels = F)
 #

####################################

 #Predict the BLUEs for chips
 library(emmeans)
 traitnames<- (c( "DMC", "STC_ch", "AMY_ch", "CF_ch", "SC_ch") )
 BLUEMean_ch = tibble()
 for(traits in traitnames){
   eval(parse(text = paste("model_ch= lmer(",traits,"~ (1|Year)  + (1|Location) + (1|Year:Location) + Genotype +  (1|Genotype:Year) + (1|Genotype: Location),
  data=Data_new)")))
   em1_ch = as.data.frame(emmeans(object = model_ch, specs = "Genotype"))

   em2_ch = cbind(Trait = traits,em1_ch[,1:2])
   BLUEMean_ch = rbind(BLUEMean_ch, em2_ch)
 }
 head(BLUEMean_ch)
 library(reshape2)
 BLUEMeanWide_ch = dcast(data = BLUEMean_ch, formula = Genotype ~ Trait, fun.aggregate = mean, value.var = "emmean")
 head(BLUEMeanWide_ch)
 BLUEMeanWide_ch = as.data.frame(BLUEMeanWide_ch)
 rownames(BLUEMeanWide_ch) <- BLUEMeanWide_ch$Genotype
 #gsub(pattern = )
 gc_ch = c()
 for(i in 1:193){
 idx = paste(c("G",i), collapse = "")
 gc_ch = c(gc_ch,idx)
 }
 BLUEMeanWide_ch$Gcode = gc_ch
 rownames(BLUEMeanWide_ch) <- BLUEMeanWide_ch$Gcode
 pc_Ch = PCA(X = BLUEMeanWide_ch[,-c(1,7)], scale.unit = T  )
 factoextra::fviz_pca_ind( pc_Ch,cex=0.8,  labelsize = 1)


 fviz_pca_biplot(X = pc_Ch,labelsize= 2, arrowsize = 0.3,
                 pointsize = 0.1, col.var = "red") +
   theme(text = element_text(size = 10))

 write.csv(x = BLUEMeanWide_ch[,c(1,7)], file = "/Users/ca384/Documents/ChinedoziRepo/GarriQuality/output/Chips_genotype_code.csv",row.names = F)

 pc_Ch$svd # pcplot

 pc_Ch$var$coord # looking for the cordinates
 pc_Ch$ind$coord
 pcsum = summary(pc_Ch)
 pcsum

 print(citation("ggplot2"), style = "text")
 print(citation("tidyverse"), style = "text")
 print(citation("factoextra"), style = "text")
