
#========================================================================#
#  Tittle:"Multivariate Analysis of Risk Factors of the COVID-19 Pandemic 
#          in the Community of Madrid (Spain)"
#  Author: Víctor Pérez-Segura
#  Date: 10/04/2021
#  Institutional  affiliation: University Institute of Studies on Migrations. 
#                              Chair of Catastrophes. Comillas Pontifical 
#                              University
#========================================================================#


setwd("/Users/JohnDoe/Desktop/paper_covid")
library(readxl)
library(extrafont)
library(tidyverse)
library(nortest)
library(psych) 
library(cluster)    
library(factoextra) 
library(patchwork)
library(MASS)
library(lmtest) 
library(missForest)
library(sandwich)
library(agricolae)
library(reshape2)
library(sf)
library(ggspatial)


#---------------------------- CONTENTS -----------------------------#
#                                                                 
# 1. Editing and Data Processing
# 1.1. Loading the raw database
# 1.2. Imputation processing
# 1.3. Box-cox transoformation
# 2. Principal Components Analysis
# 2.1. PCA`s analysis matrix
# 2.2. Matrix adecuacy PCA tests
# 2.3. Parrallel analisys
# 2.4. PCA
# 3. Cluster
# 3.1. Silhouette method
# 3.2. K-means cluster analisys
# 4. GLM & ANOVA database
# 4.1. GLM database
# 5. GLM C.Madrid
# 5.1. CA Madrid fit model
# 5.2. Linear Regression Model Diagnostics
# 6. GLM by cluster
# 6.1.1 Madrid-City fit model
# 6.1.2 Linear Regression Model Diagnostics
# 6.2.1 Madrid-Surroundings fit model
# 6.2.2 Linear Regression Model Diagnostics
# 6.3.1 North-East fit model
# 6.3.2 Linear Regression Model Diagnostics
# 7. ANOVA & Scheffe´s method
# 8. Figures
# 8.1 Figure 1: Cluster Map
# 8.2 Figure 2: Components by cluster barplot
#  8.3 Figura 3: Boxplot of COVID-19 cases by cluster


###############################################################
#  1. Editing and Data Processing  
###############################################################

###############################################################
# 1.1. Loading the raw database

# Independent variables

tfgCAM<-read_excel("mad_cov.xlsx")

tfgcam<-tfgCAM[c(-200:-204),]

raw_database <- dplyr::select(tfgcam, "Municipio / Distrito", "%Pob 65+", 
                              "SO2_Nivel_Medio_Abril", "CO_Nivel_Medio_Abril", 
                              "NO_Nivel_Medio", "Ozono_Nivel_Medio","% M Sist Resp",
                              "Densidad", "PM2.5","NO2_Nivel_Medio", "TempMedia",
                              "HumRel", "Renta per capita", "%Total afiliados","Población")

names(raw_database) <- c("municipio_distrito",  "Pob65", "SO2", "CO","NO", "Ozono", "MSistResp",
                         "Densidad", "PM25", "NO2", "TempMedia", "HumRel", "Rentapercapita", 
                         "Totalafiliados", "Pop")

population <- dplyr::select(raw_database,"municipio_distrito","Pop")


# Dependent variable
total <- read.csv("covid19_tia_muni_y_distritos.csv", sep=";")
total$fecha_informe <- as.Date(gsub(" 09:00:00", "", gsub("07:00:00", "", total$fecha_informe)))
total <- filter(total, fecha_informe == "2020-05-13")
total$casos_confirmados_totales <- as.numeric(total$casos_confirmados_totales)
total <- dplyr::select(total,"municipio_distrito", "casos_confirmados_totales")

names(total)<-c("municipio_distrito", "total")

# Merge: Dependent and independent variables
df <- merge(total, raw_database, by = "municipio_distrito")
remove(total)
remove(tfgCAM)
remove(tfgcam)
remove(raw_database)

###############################################################
# 1.2 Imputation processing

df$total_imputed <- df$total
df$total_imputed[is.na(df$total_imputed)] <- 5

###############################################################
# 1.3 Box-cox transoformation

# Logaritmic transformation
df$total_imputed_log <- log(df$total_imputed)


###############################################################
#  2. Principal Components Analysis
###############################################################

###############################################################
# 2.1. Construction of the analysis matrix

# Selection of risk factor variables
pca_df <- dplyr::select(df, -"total", -"municipio_distrito", -"total_imputed_log", 
                        -"total_imputed", -"Pop")

# The missing value in NO2 is imputed using Random Forest
imp <- missForest(pca_df)
pca_df<-imp[["ximp"]]

pca_df<-data.frame(lapply(pca_df, scale))

###############################################################
# 2.2. Matrix adecuacy PCA tests

# KMO(pca_df)
# bartlett.test(pca_df)

###############################################################
# 2.3. Parrallel analisys

 #fa.parallel(pca_df, fm="ml",fa = 'pc', main = 'Figure A.1. Parrallel Analysis Scree Plot')


###############################################################
# 2.4. PCA

pca_3 <- principal(pca_df, nfactors=3, rotate = "varimax")
#print(pca_3)

# Predict PCA values for al territories
pred.3pca<-predict(pca_3,pca_df)


###############################################################
#  3. Cluster                                                
###############################################################

rownames(pred.3pca) <- df$municipio_distrito

km3 <- kmeans(pred.3pca, centers = 3, nstart=25)


# 3.1. Elbow method

fviz_nbclust(pred.3pca, kmeans, method = "wss") +### TIMES NEW ROMAN
  labs(title= "Figure B.2. Elbow Method with the components") +
  theme(plot.title = element_text(face = "bold")) + 
  theme(plot.title = element_text(hjust = 0.5))


###############################################################
# 3.2. K-means cluster analisys

km3 <- kmeans(pred.3pca, centers = 3, nstart=25)
################################################################
#  4. GLM & ANOVA database
###############################################################

###############################################################
#  4.1. GLM database

# Components database
km3 <- as.data.frame(km3$cluster)
km3 <- rownames_to_column(km3, var = "municipio_distrito")
pred.3pca <- as.data.frame(pred.3pca)
pred.3pca <- rownames_to_column(pred.3pca, var = "municipio_distrito")
df_cluster<-merge(km3, pred.3pca, by = "municipio_distrito")
names(df_cluster)<-c("municipio_distrito","cluster","Pollution_Density","PM_Temp", "Socioeconomic")


# Target database
df_total <- dplyr::select(df, "municipio_distrito", "total_imputed_log")


# Merge: Target + components

df_models <- merge(df_total, df_cluster,by = "municipio_distrito")

# Relabing cluster variable
df_models$cluster[df_models$cluster == 1 & sum(df_models$cluster == 1)==21] <- "Madrid-City"
df_models$cluster[df_models$cluster == 1 & sum(df_models$cluster == 1)==96] <- "Madrid-Surroundings"
df_models$cluster[df_models$cluster == 1 & sum(df_models$cluster == 1)==82] <-  "North-East"

df_models$cluster[df_models$cluster == 2 & sum(df_models$cluster == 2)==21] <- "Madrid-City"
df_models$cluster[df_models$cluster == 2 & sum(df_models$cluster == 2)==96] <- "Madrid-Surroundings"
df_models$cluster[df_models$cluster == 2 & sum(df_models$cluster == 2)==82] <-  "North-East"

df_models$cluster[df_models$cluster == 3 & sum(df_models$cluster == 3)==21] <- "Madrid-City"
df_models$cluster[df_models$cluster == 3 & sum(df_models$cluster == 3)==96] <- "Madrid-Surroundings"
df_models$cluster[df_models$cluster == 3 & sum(df_models$cluster == 3)==82] <-  "North-East"



###############################################################
# 4.2. ANOVA database

clus_vdep <- dplyr::select(df_models, "municipio_distrito", "total_imputed_log", "cluster")
vindep <- dplyr::select(df, -"total", -"total_imputed")
df.var <- merge(clus_vdep, vindep, by = "municipio_distrito")

# Box-cox inverse transformation
df.var$casos<-exp(df.var$total_imputed_log.x)

# deletion of unneeded databases
remove(clus_vdep)
remove (vindep)

###############################################################
#  5. GLM C.Madrid
###############################################################
# CA Madrid model database preparation
df_C_Madrid <- df_models %>%
  merge(population, by = "municipio_distrito") %>%
  dplyr::select(-"municipio_distrito", -"cluster")

###############################################################
#  5.1. CA Madrid fit model

fit_C_Madrid<-lm(total_imputed_log~.,data=df_C_Madrid) 
summary(fit_C_Madrid)

AIC(fit_C_Madrid)
BIC(fit_C_Madrid)
logLik(fit_C_Madrid)

# HC3 coefcients
coeftest(fit_C_Madrid, vcov = vcovHC(fit_C_Madrid, type = "HC3"))

###############################################################
# 5.2. Residual Model Diagnostics 

# Independency: Durbin-Watson test
dwtest(total_imputed_log ~ ., data = df_C_Madrid)# 2.044 con el estadístico de Durbin-Watson. Si éste

# Homocedasticity: White test 
bptest(fit_C_Madrid, varformula = ~ I(Pollution_Density^2) + I(PM_Temp^2) + I(Socioeconomic^2) +
         I(Pop^2) + Pollution_Density * PM_Temp * Socioeconomic * Pop, data = df_C_Madrid)

# Normality: Shapiro-Wilk normality test
residuals <- rstandard(fit_C_Madrid)
shapiro.test(residuals) # Normality is acepted
cvm.test(residuals)
lillie.test(residuals)
qqnorm(residuals, main="Figure.B.2 Q-Q Plot")
mean(residuals)

#   Multicolinearity test: VIF
car::vif(fit_C_Madrid)

###############################################################
#  6. GLM by cluster
###############################################################

###############################################################
# 6.1. Cluster Madrid-City

# Madrid-City model database preparation
df_Mad_City<-df_models %>% 
  merge(population, by="municipio_distrito")%>%
  filter(cluster=="Madrid-City")%>% 
  dplyr::select("total_imputed_log","Pollution_Density","PM_Temp","Socioeconomic", "Pop")

###############################################################
# 6.1.1 Madrid-City fit model
fit_Mad_City<-lm(total_imputed_log~.,data=df_Mad_City)

summary(fit_Mad_City) 

# HC3 coefcients ????
coeftest(fit_Mad_City, vcov = vcovHC(fit_Mad_City, type = "HC3"))

# Model information
AIC(fit_Mad_City)
BIC(fit_Mad_City)
logLik(fit_Mad_City)


###############################################################
# 6.1.2. Residual Diagnostics 

# Independency: Durbin-Watson test
dwtest(total_imputed_log ~ ., data = df_Mad_City)# DW = 1.7683

# Homocedasticity: White test 
bptest(fit_Mad_City, ~I(Pollution_Density^2)+I(PM_Temp^2)+I(Socioeconomic^2)+I(Pop^2)
         +Pop*Pollution_Density*PM_Temp*Socioeconomic, data=df_Mad_City)

# Normality: Shapiro-Wilk normality test
residuals <- rstandard(fit_Mad_City)
shapiro.test(residuals) # No son normales. Comentarselo a Raquel
cvm.test(residuals)
lillie.test(residuals)
mean(residuals)


#   Multicolinearity test: VIF
car::vif(fit_Mad_City) # CUAL ES EL MÁXIMO QUE PUEDE TNER DE VIF. 13.01


###############################################################
# 6.2. Cluster Madrid-Surroundings

# Madrid-Surroundings model database preparation
df_Mad_Surroundings<-df_models %>%
  merge(population, by="municipio_distrito")%>%
  filter(cluster=="Madrid-Surroundings")%>%
  dplyr::select("total_imputed_log","Pollution_Density","PM_Temp","Socioeconomic","Pop")#

###############################################################
# 6.2.1 Madrid-Surroundings fit model  North-East
fit_Mad_Surroundings<-lm(total_imputed_log~.,data=df_Mad_Surroundings)
summary(fit_Mad_Surroundings)
coeftest(fit_Mad_City, vcov = vcovHC(fit_Mad_City, type = "HC3"))

AIC(fit_Mad_Surroundings)
BIC(fit_Mad_Surroundings)
logLik(fit_Mad_Surroundings)

###############################################################
# 6.2.2. Reasidual Diagnostics 

# Independency: Durbin-Watson test
dwtest(total_imputed_log ~ ., data = df_Mad_Surroundings)# DW = 1.7683

# Homocedasticity: White test 
bptest(fit_Mad_Surroundings, ~ I(Pollution_Density^2) + I(PM_Temp^2) + I(Socioeconomic^2) +
        I(Pop^2) + Pop * Pollution_Density * PM_Temp * Socioeconomic, data = df_Mad_Surroundings)

# Normality: Shapiro-Wilk normality test
residuals <- rstandard(fit_Mad_Surroundings)
shapiro.test(residuals) # No son normales. Comentarselo a Raquel
cvm.test(residuals)
lillie.test(residuals)
mean(residuals)

#   Multicolinearity test: VIF
car::vif(fit_Mad_Surroundings) # CUAL ES EL MÁXIMO QUE PUEDE TNER DE VIF. 13.01


###############################################################
# 6.3. Cluster North-East

# North-East model database preparation
df_North_East<-df_models %>%
  merge(population, by="municipio_distrito")%>%
  filter(cluster=="North-East")%>%#"North-East"
  dplyr::select("total_imputed_log","Pollution_Density","PM_Temp","Socioeconomic","Pop")#"Pop"

###############################################################
# 6.3.1 Madrid-Surroundings fit model

fit_North_East<-lm(total_imputed_log~.,data=df_North_East)
summary(fit_North_East)

# HC3 coefcients
coeftest(fit_North_East, vcov = vcovHC(fit_North_East, type = "HC3"))

#
AIC(fit_North_East)
BIC(fit_North_East)
logLik(fit_North_East)


###############################################################
# 6.3.2. Reasidul Diagnostics 

# Independency: Durbin-Watson test
dwtest(total_imputed_log ~ ., data = df_North_East)# DW = 1.7683

# Homocedasticity: White test 
bptest(fit_North_East, ~ I(Pollution_Density^2) + I(PM_Temp^2) + I(Socioeconomic^2) + I(Pop^2) +
         Pop * Pollution_Density * PM_Temp * Socioeconomic, data = df_North_East)

# Normality: Shapiro-Wilk normality test
residuals <- rstandard(fit_North_East)
shapiro.test(residuals) # No son normales. Comentarselo a Raquel
cvm.test(residuals)
lillie.test(residuals)
mean(residuals)

#   Multicolinearity test: VIF
car::vif(fit_North_East) # CUAL ES EL MÁXIMO QUE PUEDE TNER DE VIF. 13.01




###############################################################
#  7. ANOVA & Scheffe´s method
###############################################################

# Average Cases by cluster
aov_Cases = aov(casos ~ cluster, data = df.var)
summary(aov_Cases) 
scheffe.test(aov_Cases, "cluster",console=TRUE) 


# Average PM25 by cluster
aov_PM25 = aov(PM25 ~ cluster, data = df.var)
summary(aov_PM25)
scheffe.test(aov_PM25, "cluster",console=TRUE) 

# Average PM25 by cluster
aov_SO2 = aov(SO2 ~ cluster, data = df.var)
summary(aov_SO2)
scheffe.test(aov_SO2, "cluster",console=TRUE) 

# Average CO by cluster
aov_CO = aov(CO ~ cluster,data = df.var)
summary(aov_CO)
scheffe.test(aov_CO, "cluster",console=TRUE) 

# Average NO by cluster
aovNO = aov(NO ~ cluster, data = df.var)
summary(aovNO)
scheffe.test(aovNO, "cluster",console=TRUE) 

# Average NO2 by cluster
aov_NO2 = aov(NO2 ~ cluster, data = df.var)
summary(aov_NO2)
scheffe.test(aov_NO2, "cluster",console=TRUE) 

# Average Ozono by cluster
aov_Ozono = aov(Ozono ~ cluster, data = df.var)
summary(aov_Ozono)
scheffe.test(aov_Ozono, "cluster",console=TRUE) 

# Average Income by cluster
aov_Income = aov(Rentapercapita ~ cluster, data = df.var)
summary(aov_Income)
scheffe.test(aov_Income, "cluster",console=TRUE) 

# Average Workers by cluster
aov_Workers = aov(Totalafiliados ~ cluster, data = df.var)
summary(aov_Workers)
scheffe.test(aov_Workers, "cluster",console=TRUE) 

# Average Relative Humidity by cluster
aov_RH = aov(HumRel ~ cluster, data = df.var)
summary(aov_RH)
scheffe.test(aov_RH, "cluster",console=TRUE) 

# Average temperature by cluster
aov_Temp = aov(TempMedia ~ cluster, data = df.var)
summary(aov_Temp)
scheffe.test(aov_Temp, "cluster",console=TRUE) 

# Average density by cluster
aov_Density = aov(Densidad ~ cluster, data = df.var)
summary(aov_Density)
scheffe.test(aov_Density, "cluster",console=TRUE) 




###############################################################
# 8. Figures 
###############################################################

# Load fontd
loadfonts(device = "win")      

###############################################################
#  8.1 Figure 1: Cluster Map         
###############################################################

mad <- read_sf("municipios_y_distritos_madrid.shp")

names(mad)[2] <- "municipio_distrito"

map_df_models <- merge(mad, df_models, by="municipio_distrito")
map_df_models$cluster <- as.factor(map_df_models$cluster)
names(map_df_models)[5] <- "Cluster"

# Figure 1. Cluster Map
tiff('figure1.tiff', units= "mm", width= 140, height = 120, res=300, compression = 'lzw')
ggplot(map_df_models) +
  geom_sf(aes(fill = Cluster)) +
  theme_bw() + labs(title = "Figure 1. Cluster with components",
                    x = "", y = "")+
  scale_fill_brewer(palette = "Set3")+
  theme(axis.text=element_text(size=12), legend.position="bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.25, "in"), pad_y = unit(2.5, "in"),
                         style = north_arrow_fancy_orienteering) +     
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  theme(plot.title = element_text(face = "bold")) + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()



###############################################################
#  8.2 Figure 2: Components by cluster barplot
###############################################################

df.scale.var <- dplyr::select(df_models, "cluster", "Pollution_Density", "PM_Temp",      
                            "Socioeconomic")

vmed <- df.scale.var %>% group_by(cluster) %>%
  summarise(Pollution_Density = mean(Pollution_Density),
            PM_Temp = mean(PM_Temp),
            Socioeconomic = mean(Socioeconomic))

names(vmed)[2:4] <- c("Pollution & Density", "Particular Matter & Temperature", 
                      "Socio-Economic")

# Matrix transformation
miel<-melt(vmed, id.vars = "cluster")
miel$cluster<-as.factor((miel$cluster))

tiff('figure2.tiff', units= "mm", width= 190, height = 130, res=300, compression = 'lzw')
ggplot(miel, aes(fill = variable, y = value, x = cluster)) + 
  geom_bar(position = "dodge", stat = "identity") +
  labs(title="Figure 2. Components by cluster", x =" ", y = " ") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw()+
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 12),
        legend.title = element_blank())+     
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  theme(plot.title = element_text(face = "bold")) + 
  theme(plot.title = element_text(hjust = 0.5))

dev.off()

###############################################################
#  8.3 Figure 3: Boxplot of COVID-19 cases by cluster
###############################################################

vmed <- df_models %>% group_by(cluster) %>%
  summarise(total_imputed_log = mean(exp(total_imputed_log)))

names(vmed) <- c("Cluster", "Cases")
vmed[1]

map_vmed<-merge(map_df_models,vmed, by.x="cluster", by.y="Cluster")
#map_vmed$Cases<-exp(map_vmed$cases_mean)


tiff('figure3.tiff', units= "mm", width= 140, height = 80, res=300, compression = 'lzw')

ggplot(df_models,aes(y=exp(total_imputed_log), x=as.factor(cluster)))+
  geom_boxplot()+
  theme_bw()+labs(title="Figure 3. Boxplot of COVID-19 cases by cluster",
                  x =" ", y = "Cases")+
  theme(axis.text=element_text(size=12))+
  theme(plot.title = element_text(face = "bold")) + 
  theme(plot.title = element_text(hjust = 0.5))


dev.off()



