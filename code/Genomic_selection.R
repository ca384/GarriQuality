packages_used <- c("workflowr", "tidyverse", "ggplot2", "here", "stringr" ,
                   "sommer","lme4","reshape", "rrBLUP","","emmeans")
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



require(here)
here::i_am("analysis/Linear_model_garri_quality.Rmd")


set.seed(123)
Acc_pred <- tibble()
traits<- (c("SI","SP", "WAC", "BD", "AMY","CF", "STA", "SUG"))

library(emmeans)
library(lme4)
mean_Data1<- tibble()
for(Trait in traits){
  eval(parse(text = paste("ansm <- lmer(",Trait,"~ (Genotype) + Location +
    (1|Year) ,
      data= Data1)")))
  em = emmeans(object = ansm, specs = ~ Genotype)
  em = as.data.frame(em)
  m2 = cbind(Trait = paste(Trait), em[,1:2])
  mean_Data1 = rbind(mean_Data1,m2)

}

library(reshape)
Meanwide <-dcast(data = mean_Data1, formula = Genotype ~ Trait, fun.aggregate = mean,
                 value.var = "emmean")
#mean_Data1 <- Data1%>% group_by(Genotype) %>%
#summarise_if(is.numeric, mean, na.rm = TRUE)

for(Trait in traits) {
  # Creating a folder that contain 5 subset with 5times with a total of 5*5*8(traits)= 200
  fold5 = caret::createMultiFolds(y = unique(Meanwide$Genotype), k = 5, times = 5)


  for(i in 1:length(fold5)){
    index = fold5[[i]] # the index of the sample for training set
    #subset the phenotypic data
    train_geno <- droplevels(unique(Meanwide$Genotype)[index])
    train_geno_ind <- which(Meanwide$Genotype %in% train_geno)
    train.data <- droplevels(Meanwide %>%
                               filter(row_number() %in% train_geno_ind)) # subset the training set
    dim(train.data)
    test.data <- droplevels(Meanwide %>%
                              filter(!row_number() %in% train_geno_ind)) # subset the testing set
    dim(test.data)
    test.data1 <- test.data
    test.data[,traits] = NA # change the grain yield of the training set to NA value

    mod_dat <- rbind(train.data, test.data) # combine the the data set for analysis

    #####################


    eval(parse(text = paste("ans4 <- mmer(",Trait,"~1,
  random=~vsr(Genotype, Gu=Gmat1),
      rcov=~units,
      data= mod_dat)")))

    pblub = as.data.frame(ans4$U$`u:Genotype`)
    pblub = cbind(rownames(pblub),pblub)
    colnames(pblub)[1] = "Genotype"
    testind = which(pblub$Genotype %in% unique(test.data$Genotype))
    testpred = pblub[testind,]
    testpred$Genotype %in% unique(test.data$Genotype)
    testpred = testpred[order(testpred$Genotype),]

    obs_test =  test.data1 %>% group_by(Genotype) %>%
      summarise_at(.vars = Trait, .funs = mean)
    obs_test= obs_test[order(obs_test$Genotype),]
    r = cbind( Trait, predictability = round(cor(testpred[,Trait],obs_test[,Trait], use = "pairwise.complete.obs"),3))
    colnames(r)[2] <- c("predictability")
    Acc_pred = rbind(Acc_pred,r)


  }
}
#}
Acc_pred$predictability= as.numeric(Acc_pred$predictability)
ggplot(Acc_pred, mapping = aes(x = Trait, y = predictability, fill = Trait)) +
  geom_boxplot() + theme_bw()



ggplot()
