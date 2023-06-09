boxplot(Oto_22_func$wac. , main="Water absorbing capacity of garri OTO",ylab = "WAC",xlab="Genotypes")

boxplot(Umu_22_func$wac. , main="Water absorbing capacity of garri UMU",ylab="WAC", xlab="Genotypes") #Removing outliers using the boxplot method

boxplot(Oto_22_func$swelling.index, xlab="Genotypes", ylab="SI")
boxplot(Oto_22_func$swelling.power, ylab="SP",xlab="Genotypes")
boxplot(Oto_22_func$wac.,ylab="WAC", xlab="Genotypes")
boxplot(Oto_22_func$Bulk.Density, xlab="Genotypes", ylab="Bulk density")

plot(Umu_oto_22)

#removing outliers for swelling index using the quatile and interquatile ranges
quartiles_SI <- quantile(pyt_functional_1$swelling.index, probs=c(.25, .75), na.rm = TRUE) # qualtile for swelling index
IQR_SI <- IQR(pyt_functional_1$swelling.index, na.rm = TRUE) # interquatile range

Lower_SI <- quartiles_SI[1] - 1.5*IQR_SI
Upper_SI <- quartiles_SI[2] + 1.5*IQR_SI
data_no_outlier_SI <- subset(pyt_functional_1, pyt_functional_1$swelling.index > Lower_SI & pyt_functional_1$swelling.index < Upper_SI)

dim(data_no_outlier_SI)
ggplot(data = data_no_outlier_SI,aes(x=year,y= swelling.index,fill=location)) +
  geom_boxplot() + facet_grid(.~location)


# Boxplot
narrow_pyt_functional_1 <- pyt_functional_1 %>%
  dplyr::select(Location,Year, Genotype, swelling.index, swelling.power, WAC, Bulk.Density) %>%
  gather(key = trait, value = y, -c(Location,Year, Genotype))


boxplot_env_trait <- ggplot(data = narrow_pyt_functional_1 , aes(x=Location,y=y, fill=trait)) +
  geom_boxplot( stat = "boxplot",  outlier.colour = "red", outlier.shape = 16, outlier.size = 0.5, na.rm=TRUE) +
  stat_summary(fun=mean, colour="blue", geom="point", na.rm = TRUE) +
  labs(x= "Location", y= "Response values") + theme_bw() +
  theme(axis.title = element_text(colour="black",face="bold", size=12),
        plot.title = element_text(hjust = 0.5,lineheight=.5,colour="black",face="bold", size=12),
        axis.text = element_text(face="bold", size=7), axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")+
  facet_grid(trait~Year, scales = "free")

# save the plot to a file
ggsave("output/boxplot_env_trait.jpeg",height=5, width=8, units="in", dpi=300)
print(boxplot_env_trait)




#trait correlation
corr_functional <-pyt_functional_1[,-c(1:4)]
corr_functional <- corr_functional[,-c(5)]
corr_func <- cor(corr_functional,  use="complete.obs") # correlation
# correlation
corr_func
pairs(corr_func)

covar_functional  <- cov(corr_functional, use="complete.obs")
covar_functional
write.csv(corr_functional,file=here("./output/correlation_functional_all_year.csv"))

# plot using the chem_func1 that contins all the chemical data
library(corrplot)
corr <- cor(chem_func1,use = "pairwise.complete.obs")
cor_chem_func <- corrplot(corr, method = "number",type = "upper",order = "alphabet", is.corr = TRUE)



#Box plot for all traits
narrow_chem_func2 <- chem_func2 %>%
  dplyr::select(Location,Year, Genotype, swelling.index, swelling.power, WAC, Bulk.Density) %>%
  gather(key = trait, value = y, -c(Location,Year, Genotype))

#box plot of functional with amylose and crude fiber

boxplot_env_trait2 <- ggplot(data = narrow_chem_func2 , aes(x=Location,y=y, fill=trait)) +
  geom_boxplot( stat = "boxplot",  outlier.colour = "red", outlier.shape = 16, outlier.size = 0.5, na.rm=TRUE) +
  stat_summary(fun=mean, colour="blue", geom="point", na.rm = TRUE) +
  labs(x= "Location", y= "Response values") + theme_bw() +
  theme(axis.title = element_text(colour="black",face="bold", size=9),
        plot.title = element_text(hjust = 0.3,lineheight=.5,colour="black",face="bold", size=9),
        axis.text = element_text(face="bold", size=3), axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")+
  facet_grid(trait~Year, scales = "free")

# save the plot to a file
ggsave("output/boxplot_env_trait.jpeg",height=5, width=8, units="in", dpi=300)
print(boxplot_env_trait2)

# Boxplot for Amylose and Crude.Fiber
narrow_chem_func3 <- chem_func2 %>%
  dplyr::select(Location,Year, Genotype, Amylose,Crude.Fiber) %>%
  gather(key = trait, value = y, -c(Location,Year,Genotype))

boxplot_env_trait3 <- ggplot(data = narrow_chem_func3 , aes(x=Location,y=y, fill=trait)) +
  geom_boxplot( stat = "boxplot",  outlier.colour = "red", outlier.shape = 16, outlier.size = 0.5, na.rm=TRUE) +
  stat_summary(fun=mean, colour="blue", geom="point", na.rm = TRUE) +
  labs(x= "Location", y= "Response values") + theme_bw() +
  theme(axis.title = element_text(colour="black",face="bold", size=9),
        plot.title = element_text(hjust = 0.3,lineheight=.5,colour="black",face="bold", size=9),
        axis.text = element_text(face="bold", size=3), axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")+
  facet_grid(trait~Year, scales = "free")

ggsave("output/boxplot_env_trait.jpeg",height=5, width=8, units="in", dpi=300)
print(boxplot_env_trait3)


chem_func2 <-chem_func %>% dplyr::select(c( "Location","Year","Genotype", "swelling.index","swelling.power","WAC", "Bulk.Density","Amylose","Crude.Fiber"))
chem_func2$Location <- as.factor(chem_func2$Location)
str(chem_func2)
#chem_func2 <- chem_func2[-1,]# remove row 1, it is empty and was created when I sorted my data

ggplot(data = chem_func2,aes(x=Year,y= swelling.index,fill=Location)) +
  geom_boxplot()
ggplot(data = chem_func2,aes(x=Year,y= swelling.power,fill=Location)) +
  geom_boxplot()

ggplot(data = chem_func2,aes(x=Year,y= WAC,fill=Location))+
  geom_boxplot()

ggplot(data = chem_func2,aes(x=Year,y= Bulk.Density,fill=Location))+
  geom_boxplot()

ggplot(data = chem_func2,aes(x=Year,y=Crude.Fiber,fill=Location))+
  geom_boxplot()

ggplot(data = chem_func2,aes(x=Year,y= Amylose,fill=Location))+
  geom_boxplot()
