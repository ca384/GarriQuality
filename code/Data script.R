# r read in the dosage file and modify the row names to match with that in the phenotypic file
snps_dosage <- read.csv(here::here("data/snps_199_genotype.csv"), head=T)

dim(snps_dosage)
head(snps_dosage[,1:15])
snps_dosage1 <- (snps_dosage[,-1])
head(snps_dosage1[1:5,1:5])

# phenotype data
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

# change Oto and Umu_19 to Umu_20
Oto_19 <- read.csv(here::here("data", "functional-otobi2019_2020.csv"), head=T)

Umu_19 <- read.csv(here::here("data", "Umu_func_pyt2019_2020.csv"), head=T)

#read in 2020_2021 data for Otobi and Umudike
Umu_21 <- read.csv(here::here("data", "Umu_2020_2021_PYT_Functional.csv"), head=T) # input genotypes planted in 2020 and harvested in 2021

Oto_21 <- read.csv(here::here("data", "OTOBI_2020_2021_PYT_Functional.csv"), head=T)

All_chem <- read.csv(here::here("data/PYT_chemical_data/All_years_chem.csv"))
dim(All_chem)
pyt_functional <- pyt_functional %>% dplyr::rename(Location=location, Year=year)
pyt_functional <- pyt_functional %>% dplyr::mutate(combo=paste0(Location, Year, Genotype))
All_chem <- All_chem %>% dplyr::mutate(combo=paste0(Location, Year, Genotype))
str(All_chem)
All_chem <- All_chem[-1,] #remove row 1
All_chem <- All_chem[-1]
str(All_chem)

#All_chem <- All_chem %>% dplyr::mutate(Year_n=as.factor(Year)) %>% dplyr::select(-Year) %>% rename(Year=Year_n)
str(All_chem)

#All_chem <- All_chem %>% dplyr::rename(location=Location, year=Year)
All_chem <- All_chem%>% relocate(Location, Year) # change the position of the location and year
All_chem$Location <- as.factor(All_chem$Location)
All_chem$Year <- as.factor(All_chem$Year)
All_chem$Genotype <- as.factor(All_chem$Genotype)
chem_func <- dplyr::full_join(pyt_functional, All_chem, by=c("Location", "Year", "Genotype"))
# chem_func <- merge(pyt_functional, All_chem, all.x=FALSE)

dim(pyt_functional)
dim(All_chem)
str(chem_func)

levels(chem_func$Genotype)[levels(chem_func$Genotype) =="TME419"] = "TMEB419" # changing genotype name pattern

levels(chem_func$Genotype)[levels(chem_func$Genotype) =="IITA-TMS-1BA30572" ] = "IITA-TMS-IBA30572" # changing genotype name pattern
chem_func$Genotype = gsub(pattern = " ", replacement = "", x = chem_func$Genotype) # remove space from the genotype name
colnames(chem_func)

str(chem_func)
chem_func1 <- chem_func %>% dplyr::select(c("swelling.index","swelling.power","WAC","Bulk.Density","Amylose","Crude.Fiber")) # select only the numeric data
