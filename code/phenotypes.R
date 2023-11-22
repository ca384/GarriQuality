Oto_22_func <- read.csv(here::here("data", "OTOBI_2022PYT_FUNCTIONAL.csv"), head=T)

# Oto_22_func<-read.csv("/Users/ca384/Documents/ChinedoziRepo/GarriQuality/data/OTOBI_2022PYT_FUNCTIONAL.csv", head=T) # data for Otobi wit swelling power, swelling indx, bulk denSIty and water absorbing capacity for 2022
Umu_22_func <- read.csv(here::here("data/UMU_2022_PYT_FUNCTIONAL.csv"), head=T)
# Umu_22_func<-read.csv("/Users/ca384/Documents/ChinedoziRepo/GarriQuality/data/UMU_2022_PYT_FUNCTIONAL.csv",head=T)# data for Umudike wit swelling power, swelling index, bulk denSIty and water absorbing capacity for 2022

str(Oto_22_func)
str(Umu_22_func)
Oto_22_func$Genotype<-as.factor(Oto_22_func$Genotype)
Oto_22_func$wac.<-as.numeric(Oto_22_func$wac.)
Umu_22_func$Genotype<-as.factor(Umu_22_func$Genotype)
head(Umu_22_func)
head(Oto_22_func)
str(Oto_22_func)
str(Umu_22_func)
summary(Oto_22_func)
summary(Umu_22_func)

#Using rename() in dplyr

colnames(Oto_22_func)
colnames(Umu_22_func)
Umu_22_func <- Umu_22_func %>%
  rename(Bulk.Density = AverageBD) # rename bulk density

Umu_22_func <- Umu_22_func %>%
  rename(WAC= wac.) # rename wac. to WAC
Umu_22_func <- Umu_22_func %>%
  rename(swelling.index= Swelling.index) # rename wac. to WAC

Oto_22_func <- Oto_22_func %>%
  rename(WAC= wac.) # rename wac. to WAC for otobi

#remove column, sn and add column for year and location
Umu_22_func <-Umu_22_func %>%
  add_column(location = "Umudike", year="2022") # add column with location and year for Umudike

#Umu_22_func <- subset(Umu_22_func, select=-c(sn)) # Remove the sn

Umu_22_func <- Umu_22_func %>% relocate(location, year) # change the position of the location and year

Oto_22_func <- Oto_22_func %>%
  add_column(location = "Otobi", year="2022") %>% relocate(location, year)
# change the position of location and year from last column to the first column

Oto_22_func <- Oto_22_func %>% dplyr::select("location", "year",  "Genotype","swelling.index", "swelling.power", "WAC", "Bulk.Density") # select the needed columns

Umu_22_func <- Umu_22_func %>% dplyr:: select("location", "year",  "Genotype", "swelling.index", "swelling.power", "WAC", "Bulk.Density") # select the needed columns

#joining the two data frames
Umu_oto_22 <- rbind(Oto_22_func, Umu_22_func)

Umu_oto_22$location <- as.factor(Umu_oto_22$location)
Umu_oto_22$year <- as.factor(Umu_oto_22$year)

colnames(Umu_oto_22)

#importing 2019 functional properties data and renaming column
# change Oto and Umu_19 to Umu_20
Oto_19 <- read.csv(here::here("data", "functional-otobi2019_2020.csv"), head=T)

Umu_19 <- read.csv(here::here("data", "Umu_func_pyt2019_2020.csv"), head=T)

colnames(Oto_19)
colnames(Umu_19)

Oto_19 <- Oto_19 %>%
  rename(Genotype= sample_id)

Umu_19 <- Umu_19 %>%
  rename(Genotype= Genotypes)

Umu_19 <- Umu_19 %>%
  rename(Moisture = X.moisture)

Umu_19 <- Umu_19 %>%
  rename(Dry.Matter =X..Dry.Matter)

Oto_19 <- Oto_19 %>%
  rename(swelling.index=Swelling.Index..final.initial.)

Oto_19 <- Oto_19 %>%
  rename(swelling.power= SP)

Oto_19 <- Oto_19 %>%
  rename(Bulk.Density= Bulk.Density..g.ml.)

colnames(Oto_19)
colnames(Umu_19)
dim(Oto_19)

dim(Umu_19)

#select the need column from Oto_19 and Umu_19, Add year 2019, change the stucture of the data
#Rename Oto_Selected_aggre


Oto_Selected <- Oto_19 %>% dplyr::select("location", "Genotype" ,"swelling.index", "swelling.power", "WAC"          , "Bulk.Density" )


#taking the average of the genotypes duplicates for otobi
Oto_Selected_aggre <- Oto_Selected %>% group_by(Genotype) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

str(Oto_Selected_aggre)


Oto_Selected_aggre$Genotype <- as.factor(Oto_Selected_aggre$Genotype)

#Add columns for years and location
Oto_Selected_aggre  <-  Oto_Selected_aggre  %>%
  add_column(location = "Otobi",year = "2020") # change the year from 2019 to 2020

Umu_Selected <- Umu_19 %>% dplyr:: select('Genotype', 'swelling.index', 'swelling.power', 'WAC', 'Bulk.Density') # select the needed columns

Umu_Selected  <-  Umu_Selected  %>%
  add_column(location = "Umudike",year = "2020") # add location and year to the data frame

#add otobi and umudike 2019 data
Umu_oto_19 <- rbind(Oto_Selected_aggre,Umu_Selected)
str(Umu_oto_19)
Umu_oto_19$location <-  as.factor(Umu_oto_19$location)
Umu_oto_19$year <-  as.factor(Umu_oto_19$year)
Umu_oto_19 <- data.frame(Umu_oto_19 )
Umu_oto_19 <-  Umu_oto_19  %>% relocate(location, year) # change the positions of year and location

#read in 2020_2021 data for Otobi and Umudike
Umu_21 <- read.csv(here::here("data", "Umu_2020_2021_PYT_Functional.csv"), head=T) # input genotypes planted in 2020 and harvested in 2021

Oto_21 <- read.csv(here::here("data", "OTOBI_2020_2021_PYT_Functional.csv"), head=T)

dim(Umu_21)
Umu_21$Genotypes <- as.factor(Umu_21$Genotypes)
Umu_21 <- Umu_21 %>%
  rename(Genotype= Genotypes) # rename Genotypes to genotype
Umu_21 <- Umu_21 %>%
  rename(swelling.index= Swelling.index)

Umu_21$WAC <- as.numeric(Umu_21$WAC)

Oto_21$Genotype <- as.factor(Oto_21$Genotype)
dim(Oto_21)

#combining 2021 Umudike and otobi in one data frame and adding year and location
Umu_21  <-  Umu_21  %>%
  add_column(location = "Umudike",year = "2021") # add year and location to Umu_21
Oto_21  <-  Oto_21  %>%
  add_column(location = "Otobi",year = "2021") # add year and location to oto_21

Umu_oto_21 <- rbind(Oto_21, Umu_21)# combining the two data frame

Umu_oto_21 <-  Umu_oto_21 %>% relocate(location, year) # change the position of year and location

#combining all three years
pyt_functional <- rbind(Umu_oto_19, Umu_oto_21, Umu_oto_22) # combine all three years
dim(pyt_functional)
pyt_functional <- pyt_functional %>% arrange(location)
pyt_functional = as.data.frame(pyt_functional)

#editing genotype names in pyt_functional_1$Genotype to match with snp$X and using metan package to change multiple column as factor

head(pyt_functional_1)

levels(pyt_functional_1$Genotype)[levels(pyt_functional_1$Genotype) =="TME419"] = "TMEB419" # changing genotype name pattern

levels(pyt_functional_1$Genotype)[levels(pyt_functional_1$Genotype) =="IITA-TMS-1BA30572" ] = "IITA-TMS-IBA30572" # changing genotype name pattern
pyt_functional_1$Genotype = gsub(pattern = " ", replacement = "", x = pyt_functional_1$Genotype) # remove space from the genotype name


library(metan)
pyt_functional_1 <-  as.data.frame(pyt_functional_1)
pyt_functional_1 <- as_factor(.data = pyt_functional_1,
                              c("environment","location","year","Genotype"))
str(pyt_functional_1)

length(levels(pyt_functional_1$Genotype))
#exm = as_numeric(.data = pyt_functional_1, c("swelling.index","swelling.power","WAC","Bulk.DensIty"))
#str(exm)

write.csv(x = pyt_functional_1$Genotype, file = "pyt_functional_all_years.csv")

#extract each year to get the year information
oto20  <- droplevels(pyt_functional_1[pyt_functional_1$environment == "Otobi_2020",])
oto21  <- droplevels(pyt_functional_1[pyt_functional_1$environment == "Otobi_2021",])
oto22  <- droplevels(pyt_functional_1[pyt_functional_1$environment == "Otobi_2022",])

summary(oto20)
summary(oto21)
summary(oto22)

Umu20  <- droplevels(pyt_functional_1[pyt_functional_1$environment == "Umudike_2020",])
Umu21  <- droplevels(pyt_functional_1[pyt_functional_1$environment == "Umudike_2021",])
Umu22  <- droplevels(pyt_functional_1[pyt_functional_1$environment == "Umudike_2022",])

summary(Umu20)
summary(Umu21)
summary(Umu22)


#INSERT STARCH AND FREE SUGAR
umu_20_Starch_Sugar <-


# r making a function to run the mixed model for all the traits without removing outliers
traitname<- (c("swelling.index","swelling.power", "WAC", "Bulk.Density")) # matrix of column names

Result.h2org<- matrix(nrow = 1, ncol = length(traitname))
colnames(Result.h2org) <- traitname # naming the colnames of the heritability output so we can identify which trait has that heritability
rownames(Result.h2org) <- "H2Comb" # heritability for the combine analysis
Result.vcomporg = matrix(nrow = length(traitname), ncol = 5) # Result matrix for the variance components
colnames(Result.vcomporg) <- c("vg", "vgl","vgy", "vp","ve") # column names for the variance components
rownames(Result.vcomporg) <- traitname # row names of the variance components
for(trait in traitname){
  tr = paste0(trait) # the loop will be able to select  and analyse each trait
  Na_g = which(is.na(pyt_functional_1[,tr]))# identify traits with
  if(length(Na_g) == 0){                 # to see if there is any null within the length within the length
    DatWithout_na = pyt_functional_1    # if there is no na, it will use the pyt_functional_1 data frame
  }else{
    DatWithout_na = pyt_functional_1[-Na_g,]  # if there is missing data in the data frame, the use a new dataframe pyt_functional_1[-Na_g,] that has all the na removed
  }
  eval(parse(text=paste("model=lmer(formula = ",trait, " ~  Year + Location + Year:Location + (1|Genotype) + (1|Genotype:Year) +(1|Genotype:Location), data = pyt_functional_1)"))) # run the model

  ss <- summary(model)
  sv <- as.data.frame(ss$varcor) # get the variance component
  gv <- sv[sv$grp=="Genotype", "vcov"] # genetic variance
  gl <- sv[sv$grp == "Genotype:Location" , "vcov"] # genotype by location variance
  gy <- sv[sv$grp == "Genotype:Year" , "vcov"] # genotype by year variance
  ve <- sv[sv$grp == "Residual" , "vcov"] # residual variance
  vph <- gv + gl/2 + gy/3 + ve/6 # variance due to phenotype
  # make a matrix of the variance output
  H2 <- gv/vph
  Result.vcomporg[trait,1] = gv
  Result.vcomporg[trait,2] = gl
  Result.vcomporg[trait,3] = gy
  Result.vcomporg[trait,4] = ve
  Result.vcomporg[trait,5] = vph
  Result.h2org[,trait] = H2
}
Result.vcomporg
Result.h2org

# r  making a function to run the mixed model for all the traits after removing outliers

traitname<- (c("swelling.index","swelling.power", "WAC", "Bulk.Density")) # matrix of column names
Result.h2 = matrix(nrow = 1, ncol = length(traitname)) # this is how we want the result of the heritability to appear after analysis
colnames(Result.h2) <- traitname # naming the colnames of the heritability output so we can identify which trait has that heritability
rownames(Result.h2) <- "H2Comb" # heritability for the combine analysis
Result.vcomp = matrix(nrow = length(traitname), ncol = 5) # Result matrix for the variance components
colnames(Result.vcomp) <- c("vg", "vgl","vgy", "vp","ve") # column names for the variance components
rownames(Result.vcomp) <- traitname # row names of the variance components
for(trait in traitname){
  tr = paste0(trait)  # the loop will be able to select  and analyse each trait
  Na_g = which(is.na(pyt_functional_1[,tr]))# identify traits with missing data
  if(length(Na_g) == 0){                 # to see if there is any null within the length within the length
    DatWithout_na = pyt_functional_1    # if there is no na, it will use the pyt_functional_1 data frame
  }else{
    DatWithout_na = pyt_functional_1[-Na_g,]  # if there is missing data in the data frame, the use a new dataframe pyt_functional_1[-Na_g,] that has all the na removed
  }

  eval(parse(text=paste("model2=lmer(formula = ",trait, " ~  Year + Location + Year:Location + (1|Genotype) + (1|Genotype:Year) +(1|Genotype:Location), data = DatWithout_na )")))
  RES_trait<- rstudent(model2) # i am using the model that used data without missing data
  outliers<- which(abs(scale(RES_trait))>3) #Remove outlier
  Dat_no_outlier = DatWithout_na[-outliers,] # data frame without outliers
  eval(parse(text=paste("model=lmer(formula = ",trait, " ~  Year + Location + Year:Location + (1|Genotype) + (1|Genotype:Year) +(1|Genotype:Location), data = Dat_no_outlier )"))) # run the model
  ss <- summary(model)
  sv <- as.data.frame(ss$varcor) # get the variance component
  sv
  gv <- sv[sv$grp=="Genotype", "vcov"] # genetic variance
  gl <- sv[sv$grp == "Genotype:Location" , "vcov"] # genotype by location variance
  gy <- sv[sv$grp == "Genotype:Year" , "vcov"] # genotype by year variance
  ve <- sv[sv$grp == "Residual" , "vcov"] # residual variance
  vph <- gv + gl/2 + gy/3 + ve/6 # variance due to phenotype
  # make a matrix of the variance output
  H2 <- gv/vph
  Result.vcomp[trait,1] = gv
  Result.vcomp[trait,2] = gl
  Result.vcomp[trait,3] = gy
  Result.vcomp[trait,4] = ve
  Result.vcomp[trait,5] = vph
  Result.h2[,trait] = H2
}

Result.vcomp

Result.h2

Result.h2 <- data.frame(Result.h2) # making a data frame for the heritability output
Result.vcomp <- data.frame(Result.vcomp) #making a data frame for the vcomp output
T_h2 <- t(Result.h2) # transpose the heritability output

df_lmer <- cbind(Result.vcomp, T_h2 ) # bind the two data frame
df_lmer
