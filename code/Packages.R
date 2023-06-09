packages_used <- c("workflowr", "tidyverse", "ggplot2", "here", "stringr" , "sommer","lme4","metan", "rrBLUP","FactoMineR","corrplot")
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
