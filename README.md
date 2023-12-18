# GarriQuality

A [workflowr][] project.

[workflowr]: https://github.com/workflowr/workflowr
The r project contain the scrips for the data analysis in the following order.
There is a chunck that contains all the packages used.
All the file containing the functional traits for garri to be evaluated for 2022 and the two locations
We bind the two locations. (chunk 34 to chunk 114)
Data wrangling was done in from chunk 114 to 222
All three years (pyt_functional) data were combined in chunk 225 
Data wrangling from chunk 234 to 272
boxplot for the functional trait of garri 272 to 308
Remove outliers using boxplot method chunk 309 to 322
Chunks 325 to 347 get summary and histogram for each year
Chunks 351 to 401 was to make different plots for the garri data
Chunks 403 to 441 - trying out linear models using lmer function of lme4 for the garri data and calculating the variance component and heritability of each trait
Chunks 444 to 496-  made a function to run the mixed model for all the traits after removing outlier
Chunks 498 to 541- read in the dosage file and modify the row names to match with that in the phenotypic file, run PCA, A.mat
Chunks 543- 589 to removing missing data from the snps file, looking for snps that are not in the phenotype data, and rerun the PCA and A-mat
Chunks 590 to 624 - run mixed model with the new genotype data 
Chunks 628 to 688 -make a loop all mixed models using the A matrix for all traits and calculating the heritability
Chunks 671 and 696 crating a table for the narrow sense and broad sense heritabilities and their variance components.
Chunks 697 to 837 - merging the functional data with the chemical data of garri for both locations
Chunks 840 to 854 - clean the new data frame chem_pheno
Chunks 857 to 886- make plots (correlation and box plot of functional with amylose and crude fiber)
Chunk 890 to 943- use the data with amylose and crude fiber for the lmer and broad sense heritability removing outliers
Chunks 946 to 997 - using data with no outlier to make boxplot and correlation plot
Chunks 1002 to 1077 -  work on the genotype data again and selecting all phenotype with snps information and clean up the snps file using the ASRgenomics (This information was the same in Genotype.data.R)
Chunks 1081 to 1134- loop for lmer  with heritabilities and variance components for the physicochemical properties.
Chunks 1139 to 1234 - function to loop all mixed models using the GA matrix for all traits using data without outliers, also made the boxplot of the physicochemical in a loop


This is the main data used for the final analysis used in the paper
Chunks 1235 to 1337 - function to loop all mixed models using the GA matrix for all traits using data without outliers for only 2021 and 2022 data, also made the boxplot and correlation of the physicochemical in a loop

Genomic prediction
Chunk 1349 to 1642 - genomic prediction for the chemical_functional properties of garri
Chunk 1643 to 1690 -  Predict the BLUEs for garri to to make a biplot for the physicochemical traits of garri

The Chips_data_linear_model.R file contains similar analysis for the chips data 
