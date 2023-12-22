# Prevalence of untreated prediabetes and alterations of glucose metabolism in Mexico
# An analysis of nationally representative surveys spanning 2016-2021 
# Data Analysis: Carlos Alberto Fermin-Martinez, Cesar Daniel Paz-Cabrera & Omar Yaxmehen Bello-Chavolla
# Latest version of Analysis September, 2023
# For any question regarding analysis contact Omar Yaxmehen Bello-Chavolla at oyaxbell@yahoo.com.mx

#### Package load ####
#setwd("~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/Prediabetes ENSANUT")
#setwd("C:/Users/Investigacion/OneDrive - UNIVERSIDAD NACIONAL AUT?NOMA DE M?XICO/Prediabetes ENSANUT")
setwd("/Users/carlosfermin/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/Prediabetes ENSANUT")

pacman::p_load(readr, haven, tidyverse, nortest, ggstance, sf, spdep, biscale, shadowtext, jtools, spdep,
               na.tools, survey, jtools, ggpubr, ggsci, rmapshaper, gapminder, performance, gtsummary, nephro,
               remotes, ggplot2, reshape2, cowplot, ggthemes, lme4, lmerTest, glmmTMB, flextable, officer)

R.getIRR <- function(x){
  y <- (jtools::summ(x, confint = TRUE, digits = 3) %>% with(coeftable))[-1,c(1:3)]; z <- y %>% exp
  a <- paste0(sprintf("%#.3f",z[,1])," (",sprintf("%#.3f",z[,2])," – ",sprintf("%#.3f",z[,3]),")"); c("Ref",a)}
R.getOR <- function(x){
  y<-summary(x); a<-y$coefficients[2,1]; a<-exp(a)
  z<-confint(x); b<-z[2,1:2]; b<-exp(b)
  c<-c(a,b); names(c) <- c("OR", "L95", "U95"); c}
R.getOR.DM <- function(x){
  y<-summary(x); a<-y$coefficients[3,1]; a<-exp(a)
  z<-confint(x); b<-z[3,1:2]; b<-exp(b); c<-c(a,b); c}
R.getOR2 <- function(x){
  y <- (jtools::summ(x, confint = TRUE, digits = 3) %>% with(coeftable))[-1,c(1:3)]
  z <- y %>% exp; z[1,] %>% `names<-`(c("OR", "L95", "U95"))}
R.getOR.DM2 <- function(x){
  y <- (jtools::summ(x, confint = TRUE, digits = 3) %>% with(coeftable))[-1,c(1:3)]
  z <- y %>% exp; z[2,] %>% `names<-`(c("OR", "L95", "U95"))}
boxplot_indigenous <- function(x){
  x+geom_boxplot(color="black", linetype=5, outlier.alpha = 0.1) + theme_pubclean() + labs(x=NULL, fill=NULL) +
    scale_fill_manual(values=c("#CA4E46", "#353134")) +
    stat_compare_means(label.x.npc=0.25, size=3.5, label = after_stat("p.format")) +
    theme(legend.position = "bottom", axis.title.y = element_text(size=14), axis.text.x = element_blank())}

#### ENSANUT 2016 ####
ensanut_2016 <- read_csv("Bases/ensanut2016.csv"); ensanut_2016.unfiltered <- ensanut_2016
##HOMA-IR
#ensanut_2016 %>% transmute("FOLIO_INT"=row_number(), "GLU_MMOL"=valor.GLU_SUERO/18,
#                           "Ins_PMOL"=valor.INSULINA*6) %>% write_csv("Bases/ensanut.2016.HOMA.csv")
ensanut_2016 <- ensanut_2016 %>% mutate("FOLIO_INT"=row_number()); homa_2016<-read.csv("Bases/ensanut.2016.HOMA.csv")
homa_2016[homa_2016=="Unable to calculate"]<-NA; homa_2016[homa_2016=="#VALUE!"]<-NA
homa_2016$HOMA2IR<-as.numeric(homa_2016$HOMA2IR); homa_2016$HOMA2B<-as.numeric(homa_2016$HOMA2B)
homa_2016$HOMA2IR[homa_2016$HOMA2IR==0]<-NA; ensanut_2016<-ensanut_2016 %>% left_join(homa_2016, by="FOLIO_INT")


###Variable definitions
##Sex
ensanut_2016$sexo_2016[ensanut_2016$sexo.x==1]<-0 #Men
ensanut_2016$sexo_2016[ensanut_2016$sexo.x==2]<-1 #Women
table(ensanut_2016$sexo_2016, useNA = "always")  #0: Men // 1: Women

##Age
ensanut_2016$edad.x<-as.numeric(ensanut_2016$edad.x)
ensanut_2016$edad_cat_2016[ensanut_2016$edad.x>=20 & ensanut_2016$edad.x<=39]<-0
ensanut_2016$edad_cat_2016[ensanut_2016$edad.x>=40 & ensanut_2016$edad.x<=59]<-1 
ensanut_2016$edad_cat_2016[ensanut_2016$edad.x>=60]<-2
table(ensanut_2016$edad_cat_2016, useNA = "always") #0: 20-39 // 1: 40-59 // 2: >60

##BMI
ensanut_2016$imc_cat_2016[ensanut_2016$imc>=18.5 & ensanut_2016$imc<25]<-0 #Normal weight
ensanut_2016$imc_cat_2016[ensanut_2016$imc>=25 & ensanut_2016$imc<30]<-1 #Overweight
ensanut_2016$imc_cat_2016[ensanut_2016$imc>=30]<-2 #Obese
table(ensanut_2016$imc_cat_2016, useNA = "always") #0: Normal // 1: Over // 2: Obese

##Regions (4)
ensanut_2016$entidad_2016<-as.numeric(ensanut_2016$entidad.x)
ensanut_2016$cve_entidad_2016<-ensanut_2016$entidad_2016
ensanut_2016$regiones_2016[ensanut_2016$entidad_2016 %in% c(2,3,5,8,19,25,26,28)]<-1 #North
ensanut_2016$regiones_2016[ensanut_2016$entidad_2016 %in% c(1,6,10,11,14,16,18,24,32)]<-2 #Central-West
ensanut_2016$regiones_2016[ensanut_2016$entidad_2016 %in% c(9,13,15,17,21,22,29)]<-3 #Central
ensanut_2016$regiones_2016[ensanut_2016$entidad_2016 %in% c(4,7,12,20,23,27,30,31)]<-4 #South-Southeast
table(ensanut_2016$regiones_2016)

##Regions (9)
ensanut_2016$regiones2_2016[ensanut_2016$entidad_2016 %in% c(2,3,18,25,26)]<-1 #Pacific-North
ensanut_2016$regiones2_2016[ensanut_2016$entidad_2016 %in% c(5,8,19,28)]<-2 #Frontier
ensanut_2016$regiones2_2016[ensanut_2016$entidad_2016 %in% c(6,14,16)]<-3 #Pacific-Center
ensanut_2016$regiones2_2016[ensanut_2016$entidad_2016 %in% c(1,10,11,22,24,32)]<-4 #Center-North
ensanut_2016$regiones2_2016[ensanut_2016$entidad_2016 %in% c(13,29,30)]<-5 #Center
ensanut_2016$regiones2_2016[ensanut_2016$entidad_2016 %in% c(9)]<-6 #Mexico City
ensanut_2016$regiones2_2016[ensanut_2016$entidad_2016 %in% c(15)]<-7 #State of Mexico
ensanut_2016$regiones2_2016[ensanut_2016$entidad_2016 %in% c(12,17,20,21)]<-8 #Pacific-South
ensanut_2016$regiones2_2016[ensanut_2016$entidad_2016 %in% c(4,7,23,27,31)]<-9 #Peninsula
table(ensanut_2016$regiones2_2016)

##Area
ensanut_2016$area_2016[ensanut_2016$estrato_nvo.x==1]<-0 #Urban
ensanut_2016$area_2016[ensanut_2016$estrato_nvo.x==2]<-1 #Rural
ensanut_2016$area_2016[ensanut_2016$estrato_nvo.x==3]<-2 #Rural
table(ensanut_2016$area_2016)
ensanut_2016$area2_2016 <- ifelse(ensanut_2016$rural==1, 0, 1)
ensanut_2016$area2_2016 %>% table(useNA = "always") #0 = Rural // 1 = Urban

##Currently smoking
ensanut_2016$smoking <- with(ensanut_2016, ifelse(a1301a==1, 2, 1)); ensanut_2016$smoking[ensanut_2016$a1301==3]<-0
ensanut_2016$smoking[with(ensanut_2016, a1301a==3|a1301==9)]<-NA
ensanut_2016$smoking %>% table(useNA = "always") #0:Has never smoked // 1:No // 2:Yes

##Smoking intensity
ensanut_2016$smoking2 <- with(ensanut_2016, ifelse(a1301==1, 2, 1)); ensanut_2016$smoking2[ensanut_2016$a1301==3]<-0
ensanut_2016$smoking2[with(ensanut_2016, a1301==9)]<-NA
ensanut_2016$smoking2 %>% table(useNA = "always") #0:Has never smoked // 1:<100 cigarettes // 2:>100 cigarettes

#Indigenous language
ensanut_2016$lengua_indigena <- ifelse(ensanut_2016$h212==1, 1, 0)
ensanut_2016$lengua_indigena %>% table(useNA = "always") #0:No // 1:Yes

#Diagnosis of type 2 diabetes
ensanut_2016$previous_diabetes_2016<-ifelse(ensanut_2016$a301.x==1 | (ensanut_2016$a307 %in% c(1,2,3)), 1,0) #Previous DX
ensanut_2016$diabetes_biochem<-ifelse(((ensanut_2016$sanvenh>=8 & ensanut_2016$valor.GLU_SUERO>=126) | #Glucose ≥126 (fasting)
                                         (ensanut_2016$sanvenh<8 & ensanut_2016$valor.GLU_SUERO>=200) | #Glucose ≥200 (no fasting)
                                         (ensanut_2016$valor.HB1AC>=6.5)), 1,0) #HbA1c ≥6.5
ensanut_2016$diabetes_fin <- ensanut_2016$previous_diabetes_2016
ensanut_2016$diabetes_fin[!is.na(ensanut_2016$diabetes_biochem)] <- ifelse(
  (ensanut_2016$previous_diabetes_2016[!is.na(ensanut_2016$diabetes_biochem)]==1 |
     ensanut_2016$diabetes_biochem[!is.na(ensanut_2016$diabetes_biochem)]==1),1,0) 
#Age at diagnosis
ensanut_2016$a3025[ensanut_2016$a3025==99]<-NA
ensanut_2016$EDAD_DIABETES<-ensanut_2016$a3025
ensanut_2016$edad.x<-as.numeric(ensanut_2016$edad.x.x)
ensanut_2016$EDAD_DIABETES<-na.tools::na.replace(ensanut_2016$EDAD_DIABETES, ensanut_2016$edad.x)
ensanut_2016$EDAD_DIABETES[ensanut_2016$diabetes_fin==0|is.na(ensanut_2016$diabetes_fin)]<-NA
#Undiagnosed diabetes
ensanut_2016$undx_diabetes<-ifelse(((ensanut_2016$previous_diabetes_2016==1) | (ensanut_2016$diabetes_fin==0)),0,1)
table(ensanut_2016$undx_diabetes)

#Recoding
ensanut_2016$YEAR<-c(2016)

##Impaired fasting glucose Prediabetes
ensanut_2016$prediabetes_ifg_2016<-ifelse(ensanut_2016$sanvenh>=8 & ensanut_2016$valor.GLU_SUERO>=100 & 
                                            ensanut_2016$valor.GLU_SUERO<126,1,0)
ensanut_2016$prediabetes_ifg_2016[ensanut_2016$diabetes_fin==1]<-0
table(ensanut_2016$prediabetes_ifg_2016)

##Hb1Ac Prediabetes
ensanut_2016$prediabetes_hb1ac_2016<-ifelse(ensanut_2016$valor.HB1AC>=5.7 & 
                                              ensanut_2016$valor.HB1AC<6.5,1,0)
ensanut_2016$prediabetes_hb1ac_2016[ensanut_2016$diabetes_fin==1]<-0
table(ensanut_2016$prediabetes_hb1ac_2016)

#WHO and IEC criteria
ensanut_2016$ifg_who<-ifelse(ensanut_2016$sanvenh>=8 & ensanut_2016$valor.GLU_SUERO>=110 &
                               ensanut_2016$valor.GLU_SUERO<126,1,0)
ensanut_2016$a1c_iec<-ifelse(ensanut_2016$valor.HB1AC>=6 & ensanut_2016$valor.HB1AC<6.5,1,0)
ensanut_2016$ifg_who[ensanut_2016$diabetes_fin==1]<-0; ensanut_2016$a1c_iec[ensanut_2016$diabetes_fin==1]<-0

##IFG + Hb1Ac Prediabetes
ensanut_2016$prediabetes_ifg_hb1ac_2016<-ensanut_2016$prediabetes_ifg_2016 + ensanut_2016$prediabetes_hb1ac_2016*2
ensanut_2016$prediabetes_ifg_hb1ac_2016[ensanut_2016$diabetes_fin==1]<-0
ensanut_2016$prediabetes.prev_ifg_hb1ac_2016<-ifelse(ensanut_2016$prediabetes_ifg_hb1ac_2016==3,1,0); 
ensanut_2016$prediabetes.prev_hb1ac_2016<-ifelse(ensanut_2016$prediabetes_ifg_hb1ac_2016==2,1,0);
ensanut_2016$prediabetes.prev_ifg_2016<-ifelse(ensanut_2016$prediabetes_ifg_hb1ac_2016==1,1,0); 

## Prediabetes all criteria 
ensanut_2016$prediabetes_allcriteria_2016<-ifelse(ensanut_2016$prediabetes_ifg_hb1ac_2016==0,0,1)
table(ensanut_2016$prediabetes_allcriteria_2016)

## Any alteration in glucose metabolism
ensanut_2016$any_glucose<-ensanut_2016$prediabetes_ifg_hb1ac_2016
ensanut_2016$any_glucose[ensanut_2016$diabetes_fin==1]<-4
table(ensanut_2016$any_glucose)

## Insulin resistance
ensanut_2016$homa_cat<-ifelse(ensanut_2016$HOMA2IR>=2.5 & ensanut_2016$diabetes_fin==0,1,0)
ensanut_2016$homa_cat[is.na(ensanut_2016$HOMA2IR)] <- NA
ensanut_2016$homa_cat2<-ifelse(ensanut_2016$HOMA2IR>=2.5,1,0)

## Central obesity
ensanut_2016$ob_central<-ifelse((ensanut_2016$sexo_2016==1 & ensanut_2016$cintura>=80)| (ensanut_2016$sexo_2016==0 & ensanut_2016$cintura>=90),1,0)

## SLI and DISLI
SLI <- readxl::read_xlsx("Bases/irs.xlsx") %>% mutate("entidad"=as.numeric(id)) %>% transmute(entidad, IRS, IRS_cat) %>% 
  mutate(dens=c(253.9,52.8,10.8,16.1,20.8,130.0,75.6,15.1, 6163.3,14.9,201.5,55.7,148.1,106.2,760.2,81.0, 404.1,44.4,90.2,
                44.1,191.9,202.6,41.6,46.2, 52.8,16.4,97.1,44.0,336.0,112.3,58.7,21.5), entidad2=as.factor(entidad)) %>%
  mutate("DISLI"=(lm(formula=IRS~dens) %>% residuals)) %>% mutate("DISLI_cat"=cut(DISLI, breaks = c(
    -Inf, quantile(x = DISLI, probs = c(.25,.5,.75)), Inf), labels=c("Very-Low", "Low", "Moderate", "High"))) %>%
  mutate("DISLI_cat2"=ifelse(DISLI_cat %in% c("Very-Low","Low", "Moderate"),0,1))

SLI$region<-SLI$entidad; SLI$region[SLI$entidad %in% c(2,3,5,8,19,25,26,28)]<-1 #North
SLI$region[SLI$entidad %in% c(1,6,10,11,14,16,18,24,32)]<-2 #Central-West
SLI$region[SLI$entidad %in% c(9,13,15,17,21,22,29)]<-3 #Central
SLI$region[SLI$entidad %in% c(4,7,12,20,23,27,30,31)]<-4 #South-Southeast

SLI$region2[SLI$entidad %in% c(2,3,18,25,26)]<-1 #Pacific-North
SLI$region2[SLI$entidad %in% c(5,8,19,28)]<-2 #Frontier
SLI$region2[SLI$entidad %in% c(6,14,16)]<-3 #Pacific-Center
SLI$region2[SLI$entidad %in% c(1,10,11,22,24,32)]<-4 #Center-North
SLI$region2[SLI$entidad %in% c(13,29,30)]<-5 #Center
SLI$region2[SLI$entidad %in% c(9)]<-6 #Mexico City
SLI$region2[SLI$entidad %in% c(15)]<-7 #State of Mexico
SLI$region2[SLI$entidad %in% c(12,17,20,21)]<-8 #Pacific-South
SLI$region2[SLI$entidad %in% c(4,7,23,27,31)]<-9 #Peninsula

ensanut_2016$entidad <- ensanut_2016$entidad_2016; e16.SLI <- merge(ensanut_2016, SLI, by="entidad", all.x = T)
ensanut_2016$DISLI <- e16.SLI$DISLI; ensanut_2016$DISLI_cat <- e16.SLI$DISLI_cat; ensanut_2016$DISLI_cat2 <- e16.SLI$DISLI_cat2
ensanut_2016$DISLI_cat2 %>% table(useNA = "always")


#### ENSANUT 2018 #####
ensanut_2018 <- read_csv("Bases/ensanut2018.csv"); ensanut_2018.unfiltered <- ensanut_2018
##HOMA-IR
#ensanut_2018 %>% transmute("FOLIO_INT"=id, "GLU_MMOL"=VALOR_GLU_SUERO/18,
#                           "Ins_PMOL"=VALOR_INSULINA*6) %>% write_csv("Bases/ensanut.2018.HOMA.csv")
ensanut_2018$FOLIO_INT <- ensanut_2018$id; homa_2018<-read.csv("Bases/ensanut.2018.HOMA.csv")
homa_2018[homa_2018=="Unable to calculate"]<-NA; homa_2018[homa_2018=="#VALUE!"]<-NA
homa_2018$HOMA2IR<-as.numeric(homa_2018$HOMA2IR); homa_2018$HOMA2B<-as.numeric(homa_2018$HOMA2B)
homa_2018$HOMA2IR[homa_2018$HOMA2IR==0]<-NA; ensanut_2018<-ensanut_2018%>% left_join(homa_2018, by="FOLIO_INT")

###Variable definitions
##Sex
ensanut_2018$sexo_2018[ensanut_2018$SEXO.x==1]<-0 #Hombre
ensanut_2018$sexo_2018[ensanut_2018$SEXO.x==2]<-1 #Mujer
table(ensanut_2018$sexo_2018, useNA = "always")

##Age
ensanut_2018$edad_cat_2018[ensanut_2018$EDAD.x>=20 & ensanut_2018$EDAD.x<=39]<-0
ensanut_2018$edad_cat_2018[ensanut_2018$EDAD.x>=40 & ensanut_2018$EDAD.x<=59]<-1 
ensanut_2018$edad_cat_2018[ensanut_2018$EDAD.x>=60]<-2
table(ensanut_2018$edad_cat_2018, useNA = "always")

##Weight and height have different variables for participants <60 and >60 years old
ensanut_2018 %>% filter(!is.na(TALLA4_1)) %>% select(EDAD.x) %>% summary #Height in <60
ensanut_2018 %>% filter(!is.na(TALLA15_1)) %>% select(EDAD.x) %>% summary #Height in >60
ensanut_2018 %>% filter(!is.na(PESO1_1)) %>% select(EDAD.x) %>% summary #Weight in <60
ensanut_2018 %>% filter(!is.na(PESO12_1)) %>% select(EDAD.x) %>% summary #Weight in <60
##Unified Height
ensanut_2018$talla_2018 <- ensanut_2018$TALLA4_1
ensanut_2018$talla_2018[ensanut_2018$EDAD.x>=60&!is.na(ensanut_2018$TALLA15_1)] <- with(ensanut_2018, TALLA15_1[!is.na(TALLA15_1)])
ensanut_2018$talla_metros_2018<-(ensanut_2018$talla_2018/100) #En metros
##Unified Weight
ensanut_2018$peso_2018 <- ensanut_2018$PESO1_1
ensanut_2018$peso_2018[ensanut_2018$EDAD.x>=60&!is.na(ensanut_2018$PESO12_1)] <- with(ensanut_2018, PESO12_1[!is.na(PESO12_1)])

##IMC
ensanut_2018$imc_calculo_2018<-(ensanut_2018$peso_2018/ensanut_2018$talla_metros_2018^2)
ensanut_2018$imc_cat_2018[ensanut_2018$imc_calculo_2018>=18.5 & ensanut_2018$imc_calculo_2018<25]<-0 #Normal
ensanut_2018$imc_cat_2018[ensanut_2018$imc_calculo_2018>=25 & ensanut_2018$imc_calculo_2018<30]<-1 #Sobrepeso
ensanut_2018$imc_cat_2018[ensanut_2018$imc_calculo_2018>=30]<-2 #Obesidad
table(ensanut_2018$imc_cat_2018, useNA = "always")

##Regions(4)
ensanut_2018$entidad_2018<-as.numeric(ensanut_2018$ENT.x)
ensanut_2018$cve_entidad_2018<-ensanut_2018$entidad_2018
ensanut_2018$regiones_2018[ensanut_2018$entidad_2018 %in% c(2,3,5,8,19,25,26,28)]<-1 #North
ensanut_2018$regiones_2018[ensanut_2018$entidad_2018 %in% c(1,6,10,11,14,16,18,24,32)]<-2 #Central-West
ensanut_2018$regiones_2018[ensanut_2018$entidad_2018 %in% c(9,13,15,17,21,22,29)]<-3 #Central
ensanut_2018$regiones_2018[ensanut_2018$entidad_2018 %in% c(4,7,12,20,23,27,30,31)]<-4 #South-Southeast
table(ensanut_2018$regiones_2018, useNA = "always")

##Regions (9)
ensanut_2018$regiones2_2018[ensanut_2018$entidad_2018 %in% c(2,3,18,25,26)]<-1 #Pacific-North
ensanut_2018$regiones2_2018[ensanut_2018$entidad_2018 %in% c(5,8,19,28)]<-2 #Frontier
ensanut_2018$regiones2_2018[ensanut_2018$entidad_2018 %in% c(6,14,16)]<-3 #Pacific-Center
ensanut_2018$regiones2_2018[ensanut_2018$entidad_2018 %in% c(1,10,11,22,24,32)]<-4 #Center-North
ensanut_2018$regiones2_2018[ensanut_2018$entidad_2018 %in% c(13,29,30)]<-5 #Center
ensanut_2018$regiones2_2018[ensanut_2018$entidad_2018 %in% c(9)]<-6 #Mexico City
ensanut_2018$regiones2_2018[ensanut_2018$entidad_2018 %in% c(15)]<-7 #State of Mexico
ensanut_2018$regiones2_2018[ensanut_2018$entidad_2018 %in% c(12,17,20,21)]<-8 #Pacific-South
ensanut_2018$regiones2_2018[ensanut_2018$entidad_2018 %in% c(4,7,23,27,31)]<-9 #Peninsula
table(ensanut_2018$regiones2_2018)

##Area
ensanut_2018$area_2018[ensanut_2018$DOMINIO.x==1]<-0 #U
ensanut_2018$area_2018[ensanut_2018$DOMINIO.x==2]<-1 #R
ensanut_2018$area_2018 <- ifelse(ensanut_2018$area_2018==1, 0, 1)
ensanut_2018$area_2018 %>% table(useNA = "always") #0 = Rural // 1 = Urban

##Currently smoking
ensanut_2018$smoking <- with(ensanut_2018, ifelse(P13_2 %in% c(1,2), 2, 1)); ensanut_2018$smoking[ensanut_2018$P13_4==3]<-0
ensanut_2018$smoking[with(ensanut_2018, P13_2==8|is.na(P13_2))]<-NA
ensanut_2018$smoking %>% table(useNA = "always") #0:Has never smoked // 1:No // 2:Yes

##Smoking intensity
ensanut_2018$smoking2 <- ifelse(ensanut_2018$P13_1==1,2,1); ensanut_2018$smoking2[ensanut_2018$P13_4==3]<-0
ensanut_2018$smoking2[ensanut_2018$P13_1==9]<-NA
ensanut_2018$smoking2 %>% table(useNA = "always")

#Diagnosis of type 2 diabetes
ensanut_2018$previous_diabetes_2018<-ifelse(ensanut_2018$P3_1==1 | (ensanut_2018$P3_8 %in% c(1,2,3)), 1,0) #Previous DX
ensanut_2018$diabetes_biochem<-ifelse(((ensanut_2018$P5_1.y>=8 & ensanut_2018$VALOR_GLU_SUERO>=126) | #Glucose ≥126 (fasting)
                                         (ensanut_2018$P5_1.y<8 & ensanut_2018$VALOR_GLU_SUERO>=200) | #Glucose ≥200 (no fasting)
                                         (ensanut_2018$VALOR_HB1AC>=6.5)), 1,0) #HbA1c ≥6.5
ensanut_2018$diabetes_fin <- ensanut_2018$previous_diabetes_2018
ensanut_2018$diabetes_fin[!is.na(ensanut_2018$diabetes_biochem)] <- ifelse(
  (ensanut_2018$previous_diabetes_2018[!is.na(ensanut_2018$diabetes_biochem)]==1 |
     ensanut_2018$diabetes_biochem[!is.na(ensanut_2018$diabetes_biochem)]==1),1,0)
#Age at diagnosis
ensanut_2018$P3_2[ensanut_2018$P3_2==99]<-NA
ensanut_2018$EDAD_DIABETES<-ensanut_2018$P3_2
ensanut_2018$EDAD_DIABETES<-na.tools::na.replace(ensanut_2018$EDAD_DIABETES, ensanut_2018$EDAD.x)
ensanut_2018$EDAD_DIABETES[ensanut_2018$diabetes_fin==0|is.na(ensanut_2018$diabetes_fin)]<-NA
#Undiagnosed diabetes
ensanut_2018$undx_diabetes<-ifelse(((ensanut_2018$previous_diabetes_2018==1) | (ensanut_2018$diabetes_fin==0)),0,1)
table(ensanut_2018$undx_diabetes)

###Recoding
ensanut_2018$YEAR<-c(2018)

##Impaired fasting glucose Prediabetes
ensanut_2018$prediabetes_ifg_2018<-ifelse(
  ensanut_2018$P5_1.y>=8 & ensanut_2018$VALOR_GLU_SUERO>=100 & ensanut_2018$VALOR_GLU_SUERO<126,1,0)
ensanut_2018$prediabetes_ifg_2018[ensanut_2018$diabetes_fin==1]<-0
table(ensanut_2018$prediabetes_ifg_2018)

##Hb1Ac Prediabetes
ensanut_2018$prediabetes_hb1ac_2018<-ifelse(ensanut_2018$VALOR_HB1AC>=5.7 & ensanut_2018$VALOR_HB1AC<6.5,1,0)
ensanut_2018$prediabetes_hb1ac_2018[ensanut_2018$diabetes_fin==1]<-0
table(ensanut_2018$prediabetes_hb1ac_2018)

#WHO and IEC criteria
ensanut_2018$ifg_who<-ifelse(
  ensanut_2018$P5_1.y>=8 & ensanut_2018$VALOR_GLU_SUERO>=110 & ensanut_2018$VALOR_GLU_SUERO<126,1,0)
ensanut_2018$a1c_iec<-ifelse(ensanut_2018$VALOR_HB1AC>=6 & ensanut_2018$VALOR_HB1AC<6.5,1,0)
ensanut_2018$ifg_who[ensanut_2018$diabetes_fin==1]<-0; ensanut_2018$a1c_iec[ensanut_2018$diabetes_fin==1]<-0

##IFG + Hb1Ac Prediabetes
ensanut_2018$prediabetes_ifg_hb1ac_2018<-ensanut_2018$prediabetes_ifg_2018 +
  ensanut_2018$prediabetes_hb1ac_2018*2
ensanut_2018$prediabetes_ifg_hb1ac_2018[ensanut_2018$diabetes_fin==1]<-0
ensanut_2018$prediabetes.prev_ifg_hb1ac_2018<-ifelse(ensanut_2018$prediabetes_ifg_hb1ac_2018==3,1,0)
ensanut_2018$prediabetes.prev_hb1ac_2018<-ifelse(ensanut_2018$prediabetes_ifg_hb1ac_2018==2,1,0)
ensanut_2018$prediabetes.prev_ifg_2018<-ifelse(ensanut_2018$prediabetes_ifg_hb1ac_2018==1,1,0)
table(ensanut_2018$prediabetes_ifg_hb1ac_2018)

##Prediabetes all criteria 
ensanut_2018$prediabetes_allcriteria_2018<-ifelse(ensanut_2018$prediabetes_ifg_hb1ac_2018==0,0,1)
table(ensanut_2018$prediabetes_allcriteria_2018)

## Any alteration in glucose metabolism
ensanut_2018$any_glucose<-ensanut_2018$prediabetes_ifg_hb1ac_2018
ensanut_2018$any_glucose[ensanut_2018$diabetes_fin==1]<-4
table(ensanut_2018$any_glucose)

## Insulin resistance
ensanut_2018$homa_cat<-ifelse(ensanut_2018$HOMA2IR>=2.5 & ensanut_2018$diabetes_fin==0,1,0)
ensanut_2018$homa_cat[is.na(ensanut_2018$HOMA2IR)] <- NA
ensanut_2018$homa_cat2<-ifelse(ensanut_2018$HOMA2IR>=2.5,1,0)

##Waist circumference has different variables for participants <60 and >60 years old
ensanut_2018$waist.A <- ensanut_2018 %>% select(CIRCUNFERENCIA8_1, CIRCUNFERENCIA8_2) %>% apply(1, mean, na.rm=T)
ensanut_2018$waist.B <- ensanut_2018 %>% select(CINTURA21_1, CINTURA21_2) %>% apply(1, mean, na.rm=T)
ensanut_2018 %>% filter(!is.na(waist.A)) %>% select(EDAD.x) %>% summary #Height in <60
ensanut_2018 %>% filter(!is.na(waist.B)) %>% select(EDAD.x) %>% summary #Height in >60
##Unified Waist
ensanut_2018$waist_2018 <- ensanut_2018$waist.A; ensanut_2018$waist_2018[ensanut_2018$EDAD.x>=60&!is.na(
  ensanut_2018$waist.B)] <- with(ensanut_2018, waist.B[!is.na(waist.B)]); ensanut_2018$waist_2018[ensanut_2018$waist_2018<50]<-NA
## Central obesity
ensanut_2018$ob_central<-ifelse((ensanut_2018$sexo_2018==1 & ensanut_2018$waist_2018>=80)|(
  ensanut_2018$sexo_2018==0 & ensanut_2018$waist_2018>=90),1,0); ensanut_2018$ob_central %>% table(useNA = "always")

#Indigenous language
lengua_2018<-read_csv("Bases/LENGUA_2018.csv") %>%
  mutate("ID_lengua"=paste0(UPM,"_",VIV_SEL,"_",HOGAR,"_",NUMREN)) %>% select(ID_lengua,P3_11)
ensanut_2018 <- ensanut_2018 %>%
  mutate("ID_lengua"=paste0(UPM.x,"_",VIV_SEL.x,"_",HOGAR.x,"_",NUMREN.x))
ensanut_2018<-merge(ensanut_2018, lengua_2018, by="ID_lengua", all.x = T)
ensanut_2018$lengua_indigena <- ifelse(ensanut_2018$P3_11.y==1, 1, 0)
ensanut_2018$lengua_indigena %>% table(useNA = "always") #0:No // 1:Yes

## Social Lag Index
ensanut_2018$entidad <- ensanut_2018$entidad_2018; e18.SLI <- merge(ensanut_2018, SLI, by="entidad", all.x = T)
ensanut_2018$DISLI <- e18.SLI$DISLI; ensanut_2018$DISLI_cat <- e18.SLI$DISLI_cat; ensanut_2018$DISLI_cat2 <- e18.SLI$DISLI_cat2
ensanut_2018$DISLI_cat2 %>% table(useNA = "always")


#### ENSANUT 2020 ####
ensanut_2020 <- read_sav("Bases/ensanut2020.sav"); ensanut_2020.unfiltered <- ensanut_2020
##HOMA-IR
#ensanut_2020 %>% transmute("FOLIO_INT"=FOLIO_INT, "GLU_MMOL"=valor.GLU_SUERO/18,
#                           "Ins_PMOL"=valor.INSULINA*6) %>% write_csv("Bases/ensanut.2020.HOMA.csv")
ensanut_2020$FOLIO_INT <- ensanut_2020$FOLIO_INT; homa_2020<-read.csv("Bases/ensanut.2020.HOMA.csv")
homa_2020[homa_2020=="Unable to calculate"]<-NA; homa_2020[homa_2020=="#VALUE!"]<-NA
homa_2020$HOMA2IR<-as.numeric(homa_2020$HOMA2IR); homa_2020$HOMA2B<-as.numeric(homa_2020$HOMA2B)
homa_2020$HOMA2IR[homa_2020$HOMA2IR==0]<-NA; ensanut_2020<-ensanut_2020%>% left_join(homa_2020, by="FOLIO_INT")

###Variable definitions
##Sex
ensanut_2020$sexo_2020[ensanut_2020$H0302.x==1]<-0 #Men
ensanut_2020$sexo_2020[ensanut_2020$H0302.x==2]<-1 #Women
table(ensanut_2020$sexo_2020, useNA = "always") #0: Men // 1: Women

##Age
ensanut_2020$H0303.x[ensanut_2020$H0303.x==999]<-NA
ensanut_2020$edad_num<-as.numeric(ensanut_2020$H0303.x)
ensanut_2020$edad_cat_2020[ensanut_2020$edad_num>=20 & ensanut_2020$edad_num<40]<-0
ensanut_2020$edad_cat_2020[ensanut_2020$edad_num>=40 & ensanut_2020$edad_num<60]<-1
ensanut_2020$edad_cat_2020[ensanut_2020$edad_num>=60]<-2
table(ensanut_2020$edad_cat_2020, useNA = "always") #0: 20-39 // 1: 40-59 // 2: >60

##Smoking
ensanut_2020 <- ensanut_2020 %>% mutate("smoking"=case_when(ADUL1A04A==3~0, ADUL1A01==3~1, ADUL1A01%in%1:2~2))
table(ensanut_2020$smoking, useNA = "always") #0:Never // 1:Former // 2:Current (NA in laboratory subset)

##BMI
ensanut_2020$imc_calculo_2020<-(ensanut_2020$an01_1/((ensanut_2020$an04_01/100)^2))
ensanut_2020$imc_cat_2020[ensanut_2020$imc_calculo_2020>=18.5 & ensanut_2020$imc_calculo_2020<25]<-0 #Normal
ensanut_2020$imc_cat_2020[ensanut_2020$imc_calculo_2020>=25 & ensanut_2020$imc_calculo_2020<30]<-1 #Overweight
ensanut_2020$imc_cat_2020[ensanut_2020$imc_calculo_2020>=30]<-2 #Obesity
table(ensanut_2020$imc_cat_2020, useNA = "always") #0: Normal // 1: Over // 2: Obese

##Regions (4)
ensanut_2020$entidad_2020<-as.numeric(ensanut_2020$ENTIDAD.x)
ensanut_2020$cve_entidad_2020<-ensanut_2020$entidad_2020
ensanut_2020$regiones_2020[ensanut_2020$entidad_2020 %in% c(2,3,5,8,19,25,26,28)]<-1 #North
ensanut_2020$regiones_2020[ensanut_2020$entidad_2020 %in% c(1,6,10,11,14,16,18,24,32)]<-2 #Central-West
ensanut_2020$regiones_2020[ensanut_2020$entidad_2020 %in% c(9,13,15,17,21,22,29)]<-3 #Central
ensanut_2020$regiones_2020[ensanut_2020$entidad_2020 %in% c(4,7,12,20,23,27,30,31)]<-4 #South-Southeast
table(ensanut_2020$regiones_2020)

##Regions (9)
ensanut_2020$regiones2_2020[ensanut_2020$entidad_2020 %in% c(2,3,18,25,26)]<-1 #Pacific-North
ensanut_2020$regiones2_2020[ensanut_2020$entidad_2020 %in% c(5,8,19,28)]<-2 #Frontier
ensanut_2020$regiones2_2020[ensanut_2020$entidad_2020 %in% c(6,14,16)]<-3 #Pacific-Center
ensanut_2020$regiones2_2020[ensanut_2020$entidad_2020 %in% c(1,10,11,22,24,32)]<-4 #Center-North
ensanut_2020$regiones2_2020[ensanut_2020$entidad_2020 %in% c(13,29,30)]<-5 #Center
ensanut_2020$regiones2_2020[ensanut_2020$entidad_2020 %in% c(9)]<-6 #Mexico City
ensanut_2020$regiones2_2020[ensanut_2020$entidad_2020 %in% c(15)]<-7 #State of Mexico
ensanut_2020$regiones2_2020[ensanut_2020$entidad_2020 %in% c(12,17,20,21)]<-8 #Pacific-South
ensanut_2020$regiones2_2020[ensanut_2020$entidad_2020 %in% c(4,7,23,27,31)]<-9 #Peninsula
table(ensanut_2020$regiones2_2020)
runif(1,9,10)
##Area
ensanut_2020$area_2020[ensanut_2020$area_20.x==1]<-0 #Rural
ensanut_2020$area_2020[ensanut_2020$area_20.x==2]<-1 #Urban
ensanut_2020$area_2020 %>% table(useNA = "always") #0: Rural // 1: Urban

#Indigenous language
ensanut_2020$lengua_indigena <- ifelse(ensanut_2020$H0311==1, 1, 0)
ensanut_2020$lengua_indigena %>% table(useNA = "always") #0:No // 1:Yes

#Diagnosis of type 2 diabetes
ensanut_2020$previous_diabetes_2020<-ifelse((ensanut_2020$H0902A.x=="1" | ensanut_2020$H0902B=="1" | #Previous diagnosis
                                               ensanut_2020$H0902C=="1" | ensanut_2020$H0902D=="1" | 
                                               ensanut_2020$H0902E=="1" | ensanut_2020$H0902F=="1" |
                                               ensanut_2020$H0902G=="1" | ensanut_2020$H0902H=="1")==T,1,0)
ensanut_2020$diabetes_biochem<-ifelse(((ensanut_2020$san04>=8 & ensanut_2020$valor.GLU_SUERO>=126) | #Glucose ≥126 (fasting)
                                         (ensanut_2020$san04<8 & ensanut_2020$valor.GLU_SUERO>=200) | #Glucose ≥200 (no fasting)
                                         (ensanut_2020$HB1AC.Valor>=6.5)), 1,0) #HbA1c ≥6.5
ensanut_2020$diabetes_fin <- ensanut_2020$previous_diabetes_2020
ensanut_2020$diabetes_fin[!is.na(ensanut_2020$diabetes_biochem)] <- ifelse(
  (ensanut_2020$previous_diabetes_2020[!is.na(ensanut_2020$diabetes_biochem)]==1 |
     ensanut_2020$diabetes_biochem[!is.na(ensanut_2020$diabetes_biochem)]==1),1,0)
#Undiagnosed diabetes
ensanut_2020$undx_diabetes<-ifelse(((ensanut_2020$previous_diabetes_2020==1) | (ensanut_2020$diabetes_fin==0)),0,1)
table(ensanut_2020$undx_diabetes)

##Impaired fasting glucose prediabetes
ensanut_2020$prediabetes_ifg_2020<-ifelse(ensanut_2020$san04>=8 & ensanut_2020$valor.GLU_SUERO>=100 & 
                                            ensanut_2020$valor.GLU_SUERO<126,1,0)
ensanut_2020$prediabetes_ifg_2020[ensanut_2020$diabetes_fin==1]<-0
table(ensanut_2020$prediabetes_ifg_2020)

## HBA1C prediabetes
ensanut_2020$prediabetes_hb1ac_2020<-ifelse(ensanut_2020$HB1AC.Valor>=5.7 & ensanut_2020$HB1AC.Valor<6.5,1,0)
ensanut_2020$prediabetes_hb1ac_2020[ensanut_2020$diabetes_fin==1]<-0
table(ensanut_2020$prediabetes_hb1ac_2020)

#WHO and IEC criteria
ensanut_2020$ifg_who<-ifelse(ensanut_2020$san04>=8 & ensanut_2020$valor.GLU_SUERO>=110 &
                               ensanut_2020$valor.GLU_SUERO<126,1,0)
ensanut_2020$a1c_iec<-ifelse(ensanut_2020$HB1AC.Valor>=6 & ensanut_2020$HB1AC.Valor<6.5,1,0)
ensanut_2020$ifg_who[ensanut_2020$diabetes_fin==1]<-0; ensanut_2020$a1c_iec[ensanut_2020$diabetes_fin==1]<-0

###Recoding
ensanut_2020$YEAR<-c(2020)

##IFG + HBA1C Prediabetes
ensanut_2020$prediabetes_ifg_hb1ac_2020<-ensanut_2020$prediabetes_ifg_2020 +
  ensanut_2020$prediabetes_hb1ac_2020*2
ensanut_2020$prediabetes_ifg_hb1ac_2020[ensanut_2020$diabetes_fin==1]<-0
ensanut_2020$prediabetes.prev_ifg_hb1ac_2020<-ifelse(ensanut_2020$prediabetes_ifg_hb1ac_2020==3,1,0)
ensanut_2020$prediabetes.prev_hb1ac_2020<-ifelse(ensanut_2020$prediabetes_ifg_hb1ac_2020==2,1,0)
ensanut_2020$prediabetes.prev_ifg_2020<-ifelse(ensanut_2020$prediabetes_ifg_hb1ac_2020==1,1,0)

##Prediabetes all criteria 
ensanut_2020$prediabetes_allcriteria_2020<-ifelse(ensanut_2020$prediabetes_ifg_hb1ac_2020==0,0,1)
table(ensanut_2020$prediabetes_allcriteria_2020)

## Any alteration in glucose metabolism
ensanut_2020$any_glucose<-ensanut_2020$prediabetes_ifg_hb1ac_2020
ensanut_2020$any_glucose[ensanut_2020$diabetes_fin==1]<-4
table(ensanut_2020$any_glucose, useNA = "always")

## Insulin resistance
ensanut_2020$homa_cat<-ifelse(ensanut_2020$HOMA2IR>=2.5 & ensanut_2020$diabetes_fin==0,1,0)
ensanut_2020$homa_cat[is.na(ensanut_2020$HOMA2IR)] <- NA
ensanut_2020$homa_cat2<-ifelse(ensanut_2020$HOMA2IR>=2.5,1,0)

## Social Lag Index
ensanut_2020$entidad <- ensanut_2020$entidad_2020; e20.SLI <- merge(ensanut_2020, SLI, by="entidad", all.x = T)
ensanut_2020$DISLI <- e20.SLI$DISLI; ensanut_2020$DISLI_cat <- e20.SLI$DISLI_cat; ensanut_2020$DISLI_cat2 <- e20.SLI$DISLI_cat2
ensanut_2020$DISLI_cat2 %>% table(useNA = "always")


#### ENSANUT 2021 #### 
ensanut_2021 <- read_csv("Bases/ensanut2021.csv"); ensanut_2021.unfiltered <- ensanut_2021

##HOMA-IR
#ensanut_2021 %>% transmute("FOLIO_INT"=FOLIO_INT, "GLU_MMOL"=valor_GLU_SUERO/18,
#                           "Ins_PMOL"=valor_INSULINA*6) %>% write_csv("Bases/ensanut.2021.HOMA.csv")
ensanut_2021$FOLIO_INT <- ensanut_2021$FOLIO_INT; homa_2021<-read.csv("Bases/ensanut.2021.HOMA.csv")
homa_2021[homa_2021=="Unable to calculate"]<-NA; homa_2021[homa_2021=="#VALUE!"]<-NA
homa_2021$HOMA2IR<-as.numeric(homa_2021$HOMA2IR); homa_2021$HOMA2B<-as.numeric(homa_2021$HOMA2B)
homa_2021$HOMA2IR[homa_2021$HOMA2IR==0]<-NA; ensanut_2021<-ensanut_2021%>% left_join(homa_2021, by="FOLIO_INT")
homa_2016
###Variable definitions
##Sex
ensanut_2021$sexo_2021[ensanut_2021$sexo==1]<-0 #Men
ensanut_2021$sexo_2021[ensanut_2021$sexo==2]<-1 #Women
table(ensanut_2021$sexo_2021, useNA = "always") #0: Men // 1: Women

##Age
ensanut_2021$edad_num<-as.numeric(ensanut_2021$edad)
ensanut_2021$edad_cat_2021[ensanut_2021$edad_num>=20 & ensanut_2021$edad_num<=39]<-0
ensanut_2021$edad_cat_2021[ensanut_2021$edad_num>=40 & ensanut_2021$edad_num<=59]<-1
ensanut_2021$edad_cat_2021[ensanut_2021$edad_num>=60]<-2
table(ensanut_2021$edad_cat_2021, useNA = "always") #0: 20-39 // 1: 40-59 // 2: >60 

##Weight and height have different variables for participants <60 and >60 years old
ensanut_2021 <- ensanut_2021 %>% filter(edad_num>=20)
ensanut_2021 %>% filter(!is.na(an04_1)) %>% select(edad_num) %>% summary #Height in <60
ensanut_2021 %>% filter(!is.na(an15_1)) %>% select(edad_num) %>% summary #Height in >60
ensanut_2021 %>% filter(!is.na(an01_1)) %>% select(edad_num) %>% summary #Weight in <60
ensanut_2021 %>% filter(!is.na(an12_1)) %>% select(edad_num) %>% summary #Weight in >60
##Unified Height
ensanut_2021$talla_2021 <- ensanut_2021$an04_1
ensanut_2021$talla_2021[ensanut_2021$edad_num>=60&!is.na(ensanut_2021$an15_1)] <- with(ensanut_2021, an15_1[!is.na(an15_1)])
ensanut_2021$talla_metros_2021<-(ensanut_2021$talla_2021/100) #En metros
##Unified Weight
ensanut_2021$peso_2021 <- ensanut_2021$an01_1
ensanut_2021$peso_2021[ensanut_2021$edad_num>=60&!is.na(ensanut_2021$an12_1)] <- with(ensanut_2021, an12_1[!is.na(an12_1)])
##BMI
ensanut_2021$imc_calculo_2021<-(ensanut_2021$peso_2021/ensanut_2021$talla_metros_2021^2)
ensanut_2021$imc_cat_2021[ensanut_2021$imc_calculo_2021>=18.5 & ensanut_2021$imc_calculo_2021<25]<-0 #Normal weight
ensanut_2021$imc_cat_2021[ensanut_2021$imc_calculo_2021>=25 & ensanut_2021$imc_calculo_2021<30]<-1 #Overweight
ensanut_2021$imc_cat_2021[ensanut_2021$imc_calculo_2021>=30]<-2 #Obesity
table(ensanut_2021$imc_cat_2021, useNA = "always")

##Area
ensanut_2021$area_2021[ensanut_2021$estrato.x==1]<-0 #Rural
ensanut_2021$area_2021[ensanut_2021$estrato.x==2]<-1 #Urban
ensanut_2021$area_2021[ensanut_2021$estrato.x==3]<-2 #Metropolitan
ensanut_2021$area2_2021 <- ifelse(ensanut_2021$area_2021==0,0,1)
ensanut_2021$area2_2021 %>% table(useNA = "always") #0 = Rural // 1 = Urban

##Currently smoking
ensanut_2021$smoking <- with(ensanut_2021, ifelse(a1301%in%c(1,2), 2, 1)); ensanut_2021$smoking[ensanut_2021$a1305==3]<-0
ensanut_2021$smoking[with(ensanut_2021, a1301==9|a1305==9|is.na(a1301))]<-NA
ensanut_2021$smoking %>% table(useNA = "always") #0:Has never smoked // 1:No // 2:Yes

## Indigenous language
ensanut_2021$lengua_indigena <- ifelse(ensanut_2021$h0311==1, 1, 0)
ensanut_2021$lengua_indigena %>% table(useNA = "always") #0:No // 1:Yes

##Regions (4)
ensanut_2021$entidad_2021<-as.numeric(ensanut_2021$entidad.x)
ensanut_2021$cve_entidad_2021<-ensanut_2021$entidad_2021
ensanut_2021$regiones_2021[ensanut_2021$entidad_2021 %in% c(2,3,5,8,19,25,26,28)]<-1 #North
ensanut_2021$regiones_2021[ensanut_2021$entidad_2021 %in% c(1,6,10,11,14,16,18,24,32)]<-2 #Central-West
ensanut_2021$regiones_2021[ensanut_2021$entidad_2021 %in% c(9,13,15,17,21,22,29)]<-3 #Central
ensanut_2021$regiones_2021[ensanut_2021$entidad_2021 %in% c(4,7,12,20,23,27,30,31)]<-4 #South-Southeast
table(ensanut_2021$regiones_2021)

##Regions (9)
ensanut_2021$regiones2_2021[ensanut_2021$entidad_2021 %in% c(2,3,18,25,26)]<-1 #Pacific-North
ensanut_2021$regiones2_2021[ensanut_2021$entidad_2021 %in% c(5,8,19,28)]<-2 #Frontier
ensanut_2021$regiones2_2021[ensanut_2021$entidad_2021 %in% c(6,14,16)]<-3 #Pacific-Center
ensanut_2021$regiones2_2021[ensanut_2021$entidad_2021 %in% c(1,10,11,22,24,32)]<-4 #Center-North
ensanut_2021$regiones2_2021[ensanut_2021$entidad_2021 %in% c(13,29,30)]<-5 #Center
ensanut_2021$regiones2_2021[ensanut_2021$entidad_2021 %in% c(9)]<-6 #Mexico City
ensanut_2021$regiones2_2021[ensanut_2021$entidad_2021 %in% c(15)]<-7 #State of Mexico
ensanut_2021$regiones2_2021[ensanut_2021$entidad_2021 %in% c(12,17,20,21)]<-8 #Pacific-South
ensanut_2021$regiones2_2021[ensanut_2021$entidad_2021 %in% c(4,7,23,27,31)]<-9 #Peninsula
table(ensanut_2021$regiones2_2021)

#Diagnosis of type 2 diabetes
ensanut_2021$previous_diabetes_2021<-ifelse((ensanut_2021$a0301.x==1 | #Previous diagnosis
                                               ensanut_2021$a0307 %in% c(1,2,3)),1,0) #Diabetes medications
ensanut_2021$diabetes_biochem<-ifelse((ensanut_2021$san04>=8 & ensanut_2021$valor_GLU_SUERO>=126 | #≥126 (fasting)
                                         ensanut_2021$san04<8 & ensanut_2021$valor_GLU_SUERO>=200 | #≥200 (no fasting)
                                         ensanut_2021$valor_HB1AC>=6.5),1,0) #HbA1c
ensanut_2021$diabetes_fin <- ensanut_2021$previous_diabetes_2021
ensanut_2021$diabetes_fin[!is.na(ensanut_2021$diabetes_biochem)] <- ifelse(
  (ensanut_2021$previous_diabetes_2021[!is.na(ensanut_2021$diabetes_biochem)]==1 |
     ensanut_2021$diabetes_biochem[!is.na(ensanut_2021$diabetes_biochem)]==1),1,0) 
#Age at diagnosis
ensanut_2021$a0302[ensanut_2021$a0302==99]<-NA
ensanut_2021$EDAD_DIABETES<-ensanut_2021$a0302
ensanut_2021$EDAD_DIABETES<-na.tools::na.replace(ensanut_2021$EDAD_DIABETES, ensanut_2021$edad_num)
ensanut_2021$EDAD_DIABETES[ensanut_2021$diabetes_fin==0|is.na(ensanut_2021$diabetes_fin)]<-NA
#Undiagnosed diabetes
ensanut_2021$undx_diabetes<-ifelse(((ensanut_2021$previous_diabetes_2021==1) | (ensanut_2021$diabetes_fin==0)),0,1)
table(ensanut_2021$undx_diabetes)

###Recoding
ensanut_2021$YEAR<-c(2021)
##Impaired fasting glucose Prediabetes
ensanut_2021$prediabetes_ifg_2021<-ifelse(
  ensanut_2021$san04>=8 & ensanut_2021$valor_GLU_SUERO>=100 & ensanut_2021$valor_GLU_SUERO<126,1,0)
ensanut_2021$prediabetes_ifg_2021[ensanut_2021$diabetes_fin==1]<-0
table(ensanut_2021$prediabetes_ifg_2021)

##Hb1Ac Prediabetes
ensanut_2021$prediabetes_hb1ac_2021<-ifelse(
  ensanut_2021$valor_HB1AC>=5.7 & ensanut_2021$valor_HB1AC<6.5,1,0)
ensanut_2021$prediabetes_hb1ac_2021[ensanut_2021$diabetes_fin==1]<-0
table(ensanut_2021$ifg_who)

#WHO and IEC criteria
ensanut_2021$ifg_who<-ifelse(
  ensanut_2021$san04>=8 & ensanut_2021$valor_GLU_SUERO>=110 & ensanut_2021$valor_GLU_SUERO<126,1,0)
ensanut_2021$a1c_iec<-ifelse(ensanut_2021$valor_HB1AC>=6 & ensanut_2021$valor_HB1AC<6.5,1,0)
ensanut_2021$ifg_who[ensanut_2021$diabetes_fin==1]<-0; ensanut_2021$a1c_iec[ensanut_2021$diabetes_fin==1]<-0

##IFG + Hb1Ac Prediabetes
ensanut_2021$prediabetes_ifg_hb1ac_2021<-ensanut_2021$prediabetes_ifg_2021 + ensanut_2021$prediabetes_hb1ac_2021*2
ensanut_2021$prediabetes_ifg_hb1ac_2021[ensanut_2021$diabetes_fin==1]<-0
ensanut_2021$prediabetes.prev_ifg_hb1ac_2021<-ifelse(ensanut_2021$prediabetes_ifg_hb1ac_2021==3,1,0)
ensanut_2021$prediabetes.prev_hb1ac_2021<-ifelse(ensanut_2021$prediabetes_ifg_hb1ac_2021==2,1,0)
ensanut_2021$prediabetes.prev_ifg_2021<-ifelse(ensanut_2021$prediabetes_ifg_hb1ac_2021==1,1,0)
table(ensanut_2021$prediabetes_ifg_hb1ac_2021)

##Prediabetes all criteria 
ensanut_2021$prediabetes_allcriteria_2021<-ifelse(ensanut_2021$prediabetes_ifg_hb1ac_2021==0,0,1)
table(ensanut_2021$prediabetes_allcriteria_2021)

## Any alteration in glucose metabolism
ensanut_2021$any_glucose<-ensanut_2021$prediabetes_ifg_hb1ac_2021
ensanut_2021$any_glucose[ensanut_2021$diabetes_fin==1]<-4
table(ensanut_2021$any_glucose)

## Insulin resistance
ensanut_2021$homa_cat<-ifelse(ensanut_2021$HOMA2IR>=2.5 & ensanut_2021$diabetes_fin==0,1,0)
ensanut_2021$homa_cat[is.na(ensanut_2021$HOMA2IR)] <- NA
ensanut_2021$homa_cat2<-ifelse(ensanut_2021$HOMA2IR>=2.5,1,0)

##Waist circumference has different variables for participants <60 and >60 years old
ensanut_2021$waist.A <- ensanut_2021 %>% select(an08_1, an08_2) %>% apply(1, mean, na.rm=T)
ensanut_2021$waist.B <- ensanut_2021 %>% select(an21_1, an21_2) %>% apply(1, mean, na.rm=T)
ensanut_2021 %>% filter(!is.na(waist.A)) %>% select(edad_num) %>% summary #Height in <60
ensanut_2021 %>% filter(!is.na(waist.B)) %>% select(edad_num) %>% summary #Height in >60
##Unified Waist
ensanut_2021$waist_2021 <- ensanut_2021$waist.A; ensanut_2021$waist_2021[ensanut_2021$edad_num>=60&!is.na(
  ensanut_2021$waist.B)] <- with(ensanut_2021, waist.B[!is.na(waist.B)])
ensanut_2021$waist_2021[ensanut_2021$waist_2021<50]<-NA
## Central obesity
ensanut_2021$ob_central<-ifelse((ensanut_2021$sexo_2021==1 & ensanut_2021$waist_2021>=80)|(
  ensanut_2021$sexo_2021==0 & ensanut_2021$waist_2021>=90),1,0); ensanut_2021$ob_central %>% table(useNA = "always")

## Social Lag Index
ensanut_2021$entidad <- ensanut_2021$entidad_2021; e21.SLI <- merge(ensanut_2021, SLI, by="entidad", all.x = T)
ensanut_2021$DISLI <- e21.SLI$DISLI; ensanut_2021$DISLI_cat <- e21.SLI$DISLI_cat; ensanut_2021$DISLI_cat2 <- e21.SLI$DISLI_cat2
ensanut_2021$DISLI_cat2 %>% table(useNA = "always")


#### ENSANUT 2022 #### 
ensanut_2022 <- read_csv("Bases/ensanut2022.csv"); ensanut_2022.unfiltered <- ensanut_2022
##HOMA-IR
#ensanut_2022 %>% transmute("FOLIO_INT"=FOLIO_INT, "GLU_MMOL"=valor_GLU_SUERO/18,
#                           "Ins_PMOL"=valor_INSULINA*6) %>% write_csv("Bases/ensanut.2022.HOMA.csv")
ensanut_2022$FOLIO_INT <- ensanut_2022$FOLIO_INT; homa_2022<-read.csv("Bases/ensanut.2022.HOMA.csv")
homa_2022[homa_2022=="Unable to calculate"]<-NA; homa_2022[homa_2022=="#VALUE!"]<-NA
homa_2022$HOMA2IR<-as.numeric(homa_2022$HOMA2IR); homa_2022$HOMA2B<-as.numeric(homa_2022$HOMA2B)
homa_2022$HOMA2IR[homa_2022$HOMA2IR==0]<-NA; ensanut_2022<-ensanut_2022%>% left_join(homa_2022, by="FOLIO_INT")

###Variable definitions
##Sex
ensanut_2022$sexo_2022[ensanut_2022$sexo==1]<-0 #Men
ensanut_2022$sexo_2022[ensanut_2022$sexo==2]<-1 #Women
table(ensanut_2022$sexo_2022, useNA = "always") #0: Men // 1: Women

##Age
ensanut_2022$edad_num<-as.numeric(ensanut_2022$edad)
ensanut_2022$edad_cat_2022[ensanut_2022$edad_num>=20 & ensanut_2022$edad_num<=39]<-0
ensanut_2022$edad_cat_2022[ensanut_2022$edad_num>=40 & ensanut_2022$edad_num<=59]<-1
ensanut_2022$edad_cat_2022[ensanut_2022$edad_num>=60]<-2
table(ensanut_2022$edad_cat_2022, useNA = "always") #0: 20-39 // 1: 40-59 // 2: >60 

##Weight, height and waist
#Remove missing values
ensanut_2022 <- ensanut_2022 %>% mutate(
  "an01_1"=case_when(an01_1==222.22~NA, an01_1!=222.22~an01_1),
  "an12_1"=case_when(an12_1==222.22~NA, an12_1!=222.22~an12_1),
  "an04_1"=case_when(an04_1==222.2~NA, an04_1!=222.2~an04_1),
  "an15_1"=case_when(an15_1==222.2~NA, an15_1!=222.2~an15_1),
  "an08_1"=case_when(an08_1==222.2~NA, an08_1!=222.2~an08_1),
  "an21_1"=case_when(an21_1==222.2~NA, an21_1!=222.2~an21_1),
  "an08_2"=case_when(an08_2==222.2~NA, an08_2!=222.2~an08_2),
  "an21_2"=case_when(an21_2==222.2~NA, an21_2!=222.2~an21_2))

##Average of two waist measurements
ensanut_2022$waist.A <- ensanut_2022 %>% select(an08_1 , an08_2 ) %>% apply(1, mean, na.rm=T)
ensanut_2022$waist.B <- ensanut_2022 %>% select(an21_1 , an21_2 ) %>% apply(1, mean, na.rm=T)
##Weight, height and waist have different variables for participants <60 and >60 years old
##Unified Height
ensanut_2022$talla_2022 <- ensanut_2022$an04_1 
ensanut_2022$talla_2022[ensanut_2022$edad_num>=60&!is.na(ensanut_2022$an15_1 )] <- with(
  ensanut_2022, an15_1 [ensanut_2022$edad_num>=60&!is.na(ensanut_2022$an15_1 )])
ensanut_2022$talla_metros_2022<-(ensanut_2022$talla_2022/100) #En metros
##Unified Weight
ensanut_2022$peso_2022 <- ensanut_2022$an01_1 
ensanut_2022$peso_2022[ensanut_2022$edad_num>=60&!is.na(ensanut_2022$an12_1 )] <- with(
  ensanut_2022, an12_1 [ensanut_2022$edad_num>=60&!is.na(ensanut_2022$an12_1 )])
##Unified Waist
ensanut_2022$waist_2022 <- ensanut_2022$waist.A
ensanut_2022$waist_2022[ensanut_2022$edad_num>=60&!is.na(ensanut_2022$waist.B)] <- with(
  ensanut_2022, waist.B[ensanut_2022$edad_num>=60&!is.na(ensanut_2022$waist.B)])
## Central obesity
ensanut_2022$ob_central<-ifelse((ensanut_2022$sexo_2022==1 & ensanut_2022$waist_2022>=80)|(
  ensanut_2022$sexo_2022==0 & ensanut_2022$waist_2022>=90),1,0); ensanut_2022$ob_central %>% table(useNA = "always")
##BMI
ensanut_2022$imc_calculo_2022<-(ensanut_2022$peso_2022/ensanut_2022$talla_metros_2022^2)
ensanut_2022$imc_cat_2022[ensanut_2022$imc_calculo_2022>=18.5 & ensanut_2022$imc_calculo_2022<25]<-0 #Normal weight
ensanut_2022$imc_cat_2022[ensanut_2022$imc_calculo_2022>=25 & ensanut_2022$imc_calculo_2022<30]<-1 #Overweight
ensanut_2022$imc_cat_2022[ensanut_2022$imc_calculo_2022>=30]<-2 #Obesity
table(ensanut_2022$imc_cat_2022, useNA = "always")
ensanut_2022$imc_cat2 <- ensanut_2022$imc_calculo_2022 %>% cut(c(-Inf, 18.5, 25, 30, Inf))
table(ensanut_2022$imc_cat2, useNA = "always")

##Area
ensanut_2022$area_2022[ensanut_2022$estrato.y==1]<-0 #Rural
ensanut_2022$area_2022[ensanut_2022$estrato.y==2]<-1 #Urban
ensanut_2022$area_2022[ensanut_2022$estrato.y==3]<-2 #Metropolitan
ensanut_2022$area2_2022 <- ifelse(ensanut_2022$area_2022==0,0,1)
ensanut_2022$area2_2022 %>% table(useNA = "always") #0 = Rural // 1 = Urban

##Region
ensanut_2022$entidad_2022<-as.numeric(ensanut_2022$entidad1)
ensanut_2022$cve_entidad_2022<-ensanut_2022$entidad_2022
ensanut_2022$regiones_2022[ensanut_2022$entidad_2022 %in% c(2,3,5,8,19,25,26,28)]<-1 #North
ensanut_2022$regiones_2022[ensanut_2022$entidad_2022 %in% c(1,6,10,11,14,16,18,24,32)]<-2 #Central-West
ensanut_2022$regiones_2022[ensanut_2022$entidad_2022 %in% c(9,13,15,17,21,22,29)]<-3 #Central
ensanut_2022$regiones_2022[ensanut_2022$entidad_2022 %in% c(4,7,12,20,23,27,30,31)]<-4 #South-Southeast
table(ensanut_2022$regiones_2022)

##Currently smoking
ensanut_2022$smoking <- with(ensanut_2022, ifelse(a1301%in%c(1,2), 2, 1)); ensanut_2022$smoking[ensanut_2022$a1305==3]<-0
ensanut_2022$smoking[with(ensanut_2022, a1301==9|a1305==9|is.na(a1301))]<-NA
ensanut_2022$smoking %>% table(useNA = "always") #0:Has never smoked // 1:No // 2:Yes

## Indigenous language
ensanut_2022$lengua_indigena <- ifelse(ensanut_2022$h0311==1, 1, 0)
ensanut_2022$lengua_indigena %>% table(useNA = "always") #0:No // 1:Yes

#Diagnosis of type 2 diabetes
ensanut_2022$previous_diabetes_2022<-ifelse((ensanut_2022$a0301==1 | #Previous diagnosis
                                               ensanut_2022$a0307 %in% c(1,2,3)),1,0) #Diabetes medications
ensanut_2022$diabetes_biochem<-ifelse((ensanut_2022$san04>=8 & ensanut_2022$valor_GLU_SUERO>=126 | #≥126 (fasting)
                                         ensanut_2022$san04<8 & ensanut_2022$valor_GLU_SUERO>=200 | #≥200 (no fasting)
                                         ensanut_2022$valor_HB1AC>=6.5),1,0) #HbA1c
ensanut_2022$diabetes_fin <- ensanut_2022$previous_diabetes_2022
ensanut_2022$diabetes_fin[!is.na(ensanut_2022$diabetes_biochem)] <- ifelse(
  (ensanut_2022$previous_diabetes_2022[!is.na(ensanut_2022$diabetes_biochem)]==1 |
     ensanut_2022$diabetes_biochem[!is.na(ensanut_2022$diabetes_biochem)]==1),1,0) 
#Age at diagnosis
ensanut_2022$a0302[ensanut_2022$a0302==99]<-NA
ensanut_2022$EDAD_DIABETES<-ensanut_2022$a0302
ensanut_2022$EDAD_DIABETES<-na.tools::na.replace(ensanut_2022$EDAD_DIABETES, ensanut_2022$edad_num)
ensanut_2022$EDAD_DIABETES[ensanut_2022$diabetes_fin==0|is.na(ensanut_2022$diabetes_fin)]<-NA
#Undiagnosed diabetes
ensanut_2022$undx_diabetes<-ifelse(((ensanut_2022$previous_diabetes_2022==1) | (ensanut_2022$diabetes_fin==0)),0,1)
table(ensanut_2022$previous_diabetes_2022, useNA = "always")
table(ensanut_2022$diabetes_biochem, useNA = "always")
table(ensanut_2022$diabetes_fin, useNA = "always")

###Recoding
ensanut_2022$YEAR<-c(2022)

#PREDIABETES
##Impaired fasting glucose Prediabetes
ensanut_2022$prediabetes_ifg_2022<-ifelse(
  ensanut_2022$san04>=8 & between(ensanut_2022$valor_GLU_SUERO, 100, 125), 1,0)
ensanut_2022$prediabetes_ifg_2022[ensanut_2022$diabetes_fin==1]<-0
table(ensanut_2022$prediabetes_ifg_2022)
##Hb1Ac Prediabetes
ensanut_2022$prediabetes_hb1ac_2022<-ifelse(
  between(ensanut_2022$valor_HB1AC, 5.7, 6.4),1,0)
ensanut_2022$prediabetes_hb1ac_2022[ensanut_2022$diabetes_fin==1]<-0
table(ensanut_2022$prediabetes_hb1ac_2022)
#WHO and IEC criteria
ensanut_2022$ifg_who<-ifelse(
  ensanut_2022$san04>=8 & ensanut_2022$valor_GLU_SUERO>=110 & ensanut_2022$valor_GLU_SUERO<126,1,0)
ensanut_2022$a1c_iec<-ifelse(ensanut_2022$valor_HB1AC>=6 & ensanut_2022$valor_HB1AC<6.5,1,0)
ensanut_2022$ifg_who[ensanut_2022$diabetes_fin==1]<-0; ensanut_2022$a1c_iec[ensanut_2022$diabetes_fin==1]<-0
table(ensanut_2022$ifg_who); table(ensanut_2022$a1c_iec)
##IFG + Hb1Ac Prediabetes
ensanut_2022$prediabetes_ifg_hb1ac_2022<-ensanut_2022$prediabetes_ifg_2022 + ensanut_2022$prediabetes_hb1ac_2022*2
ensanut_2022$prediabetes_ifg_hb1ac_2022[ensanut_2022$diabetes_fin==1]<-0
ensanut_2022$prediabetes.prev_ifg_hb1ac_2022<-ifelse(ensanut_2022$prediabetes_ifg_hb1ac_2022==3,1,0) #Both
ensanut_2022$prediabetes.prev_hb1ac_2022<-ifelse(ensanut_2022$prediabetes_ifg_hb1ac_2022==2,1,0) #A1c only
ensanut_2022$prediabetes.prev_ifg_2022<-ifelse(ensanut_2022$prediabetes_ifg_hb1ac_2022==1,1,0) #IFG only
table(ensanut_2022$prediabetes_ifg_hb1ac_2022)
##Prediabetes all criteria 
ensanut_2022$prediabetes_allcriteria_2022<-ifelse(ensanut_2022$prediabetes_ifg_hb1ac_2022==0,0,1)
table(ensanut_2022$prediabetes_allcriteria_2022)

## Any alteration in glucose metabolism
ensanut_2022$any_glucose<-ensanut_2022$prediabetes_ifg_hb1ac_2022
ensanut_2022$any_glucose[ensanut_2022$diabetes_fin==1]<-4
table(ensanut_2022$any_glucose)

## Insulin resistance
ensanut_2022$homa_cat<-ifelse(ensanut_2022$HOMA2IR>=2.5 & ensanut_2022$diabetes_fin==0,1,0)
ensanut_2022$homa_cat[is.na(ensanut_2022$HOMA2IR)] <- NA
ensanut_2022$homa_cat2<-ifelse(ensanut_2022$HOMA2IR>=2.5,1,0)

## Social Lag Index
ensanut_2022$entidad <- ensanut_2022$entidad_2022; e21.SLI <- merge(ensanut_2022, SLI, by="entidad", all.x = T)
ensanut_2022$DISLI <- e21.SLI$DISLI; ensanut_2022$DISLI_cat <- e21.SLI$DISLI_cat; ensanut_2022$DISLI_cat2 <- e21.SLI$DISLI_cat2
ensanut_2022$DISLI_cat2 %>% table(useNA = "always")


#### SVY DESIGNS  #### 

#Family history of diabetes 2016
ensanut_2016$T2D.AHF <- (ensanut_2016 %>% transmute(
  "T2D.PAP"=case_when(a701a==1~1,a701a==2~0), "T2D.MAM"=case_when(a701b==1~1,a701b==2~0),
  "T2D.HER"=NA) %>% apply(1, sum, na.rm=T)>0) %>% ifelse(1,0)
table(ensanut_2016$T2D.AHF, useNA = "always") %>% prop.table( )
#Family history of diabetes 2018
ensanut_2018$T2D.AHF <- (ensanut_2018 %>% transmute(
  "T2D.PAP"=case_when(P7_1_1.x==1~1,P7_1_1.x==2~0), "T2D.MAM"=case_when(P7_1_2.x==1~1,P7_1_2.x==2~0),
  "T2D.HER"=case_when(P7_1_3.x==1~1,P7_1_3.x==2~0)) %>% apply(1, sum, na.rm=T)>0) %>% ifelse(1,0)
table(ensanut_2018$T2D.AHF, useNA = "always") %>% prop.table( )
#Family history of diabetes 2021
ensanut_2021$T2D.AHF <- (ensanut_2021 %>% transmute(
  "T2D.PAP"=case_when(a0701p==1~1,a0701p==2~0), "T2D.MAM"=case_when(a0701m==1~1,a0701m==2~0),
  "T2D.HER"=case_when(a0701h==1~1,a0701h==2~0)) %>% apply(1, sum, na.rm=T)>0) %>% ifelse(1,0)
table(ensanut_2021$T2D.AHF, useNA = "always") %>% prop.table( )
#Family history of diabetes 2022
ensanut_2022$T2D.AHF <- (ensanut_2022 %>% transmute(
  "T2D.PAP"=case_when(a0701p==1~1,a0701p==2~0), "T2D.MAM"=case_when(a0701m==1~1,a0701m==2~0),
  "T2D.HER"=case_when(a0701h==1~1,a0701h==2~0)) %>% apply(1, sum, na.rm=T)>0) %>% ifelse(1,0)
table(ensanut_2022$T2D.AHF, useNA = "always") %>% prop.table( )

### 2016 ###
ensanut_2016<-ensanut_2016 %>% mutate("ID" = row_number()) %>% filter(edad.x>=20) #F1
ensanut_2016_2<-ensanut_2016%>%filter(!is.na(ponde_f_vv), !is.na(est_var)) #F2
ensanut_2016_survey2 <- svydesign(data=ensanut_2016_2, id=~ID, strata=~est_var, weights=~ponde_f_vv, nest=TRUE) #SD1
ensanut_2016_survey<- subset(ensanut_2016_survey2, (sanvenh>=8)&(a301.x!=2)&!is.na(valor.HB1AC)&!is.na(valor.GLU_SUERO)) #SD2
ensanut_2016_survey %>% nrow
#Diabetes prevalence
svymean(~previous_diabetes_2016, ensanut_2016_survey,na.rm = T)*100 #Diagnosed
svymean(~diabetes_fin, ensanut_2016_survey,na.rm = T)*100 #Total

### 2018 ###
ensanut_2018 <- ensanut_2018 %>% mutate("ID" = row_number()) %>% dplyr::filter(!is.na(F_20MAS)) %>% filter(EDAD.x>=20) #F1
ensanut_2018_2 <- ensanut_2018 %>% filter(!is.na(ponderador_glucosa), !is.na(ESTRATO.y)) #F2
ensanut_2018_survey2 <- svydesign(data=ensanut_2018_2, id=~ID, strata=~ESTRATO.y, weights=~ponderador_glucosa, nest=TRUE) #SD1
ensanut_2018_survey<-subset(ensanut_2018_survey2, (P5_1.y>=8)&(P3_1!=2)&!is.na(VALOR_HB1AC)&!is.na(VALOR_GLU_SUERO)) #SD2
ensanut_2018_survey %>% nrow
#Diabetes prevalence
svymean(~previous_diabetes_2018, ensanut_2018_survey,na.rm = T)*100 #Diagnosed
svymean(~diabetes_fin, ensanut_2018_survey,na.rm = T)*100 #Total

### 2020 ###
ensanut_2020 <- ensanut_2020 %>% mutate("ID" = row_number()) %>% filter(H0303.x>=20) #F1
ensanut_2020_2 <- ensanut_2020 %>% filter(!is.na(ponde_g20.y), !is.na(est_sel.x)) %>% #F2
  filter(!est_sel.x%in%(which((est_sel.x%>%table)<2)%>%names%>%as.numeric))
ensanut_2020_survey2 <- svydesign(data=ensanut_2020_2, id=~ID, strata=~est_sel.x, weights=~ponde_g20.y, nest=TRUE) #SD2
ensanut_2020_survey<- subset(ensanut_2020_survey2, (san04>=8)&!is.na(HB1AC.Valor)&!is.na(valor.GLU_SUERO)) #SD2
ensanut_2020_survey %>% nrow
#Diabetes prevalence
svymean(~previous_diabetes_2020, ensanut_2020_survey,na.rm = T)*100 #Diagnosed
svymean(~diabetes_fin, ensanut_2020_survey,na.rm = T)*100 #Total

### 2021 ###
ensanut_2021<-ensanut_2021 %>% mutate("ID" = row_number()) %>% filter(edad_num>=20) #F1
ensanut_2021_2 <- ensanut_2021 %>% filter(!is.na(ponde_vv), !is.na(est_sel.lab)) %>% #F2
  filter(!est_sel.lab%in%(which((est_sel.lab%>%table)<2)%>%names%>%as.numeric)) 
ensanut_2021_survey2<-svydesign(data=ensanut_2021_2, id=~ID, strata=~est_sel.lab, weights=~ponde_vv, nest=TRUE) #SD1
ensanut_2021_survey<-subset(ensanut_2021_survey2, (san04>=8)&(a0301.x!=2)&!is.na(valor_HB1AC)&!is.na(valor_GLU_SUERO)) #SD2
nrow(ensanut_2021_survey)
#Diabetes prevalence
svymean(~previous_diabetes_2021, ensanut_2021_survey,na.rm = T)*100 #Diagnosed
svymean(~diabetes_fin, ensanut_2021_survey,na.rm = T)*100 #Total

### 2022 ###
ensanut_2022 <- ensanut_2022 %>% mutate("ID" = row_number()) %>% filter(edad_num>=20) #F1
ensanut_2022_2 <- ensanut_2022 %>% filter(!is.na(ponde_v), !is.na(est_sel.x.x.x)) #F2
ensanut_2022_survey2<-svydesign(data=ensanut_2022_2, id=~ID, strata=~est_sel.x.x.x, weights=~ponde_v, nest=TRUE) #SD1
ensanut_2022_survey<-subset(ensanut_2022_survey2, (san04>=8)&(a0301!=2)&!is.na(valor_HB1AC)&!is.na(valor_GLU_SUERO)) #SD2
nrow(ensanut_2022_survey)
#Diabetes prevalence
svymean(~previous_diabetes_2022, ensanut_2022_survey,na.rm = T)*100 #Diagnosed
svymean(~diabetes_fin, ensanut_2022_survey,na.rm = T)*100 #Total



####---------------------------- PREDIABETES --------------------------#### ----####
#### HbA1c (prevalence by year) ####
###ENSANUT 2016
prediabetes_prevalence_2016_year<-svyby(~prediabetes_hb1ac_2016, by=~YEAR, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="High HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2016_year

###ENSANUT 2018
prediabetes_prevalence_2018_year<-svyby(~prediabetes_hb1ac_2018, by=~YEAR, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="High HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2018_year

###ENSANUT 2020
prediabetes_prevalence_2020_year<-svyby(~prediabetes_hb1ac_2020, by=~YEAR, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="High HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2020_year

###ENSANUT 2021
prediabetes_prevalence_2021_year<-svyby(~prediabetes_hb1ac_2021, by=~YEAR, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="High HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2021_year

###ENSANUT 2022
prediabetes_prevalence_2022_year<-svyby(~prediabetes_hb1ac_2022, by=~YEAR, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="High HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2022_year

prev_hba1c<-rbind(prediabetes_prevalence_2016_year,prediabetes_prevalence_2018_year,prediabetes_prevalence_2020_year,
                  prediabetes_prevalence_2021_year,prediabetes_prevalence_2022_year)
prev_hba1c<-as.data.frame(prev_hba1c)
prev_hba1c$YEAR<-rownames(prev_hba1c)
prev_hba1c

#### IFG (prevalence by year) ####
###ENSANUT 2016
prediabetes_prevalence_2016_year<-svyby(~prediabetes_ifg_2016, by=~YEAR, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2016_year

###ENSANUT 2018
prediabetes_prevalence_2018_year<-svyby(~prediabetes_ifg_2018, by=~YEAR, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2018_year

###ENSANUT 2020
prediabetes_prevalence_2020_year<-svyby(~prediabetes_ifg_2020, by=~YEAR, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2020_year

###ENSANUT 2021
prediabetes_prevalence_2021_year<-svyby(~prediabetes_ifg_2021, by=~YEAR, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2021_year

###ENSANUT 2022
prediabetes_prevalence_2022_year<-svyby(~prediabetes_ifg_2022, by=~YEAR, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2022_year

prev_ifg<-rbind(prediabetes_prevalence_2016_year,prediabetes_prevalence_2018_year,prediabetes_prevalence_2020_year,
                prediabetes_prevalence_2021_year,prediabetes_prevalence_2022_year)
prev_ifg<-as.data.frame(prev_ifg)
prev_ifg$YEAR<-rownames(prev_ifg)
prev_ifg

#### IFG AND High HbA1c (prevalence by year) ####
###ENSANUT 2016
prediabetes_prevalence_2016_year<-svyby(~prediabetes.prev_ifg_hb1ac_2016, by=~YEAR, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_hb1ac_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="High HbA1c-IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2016_year

###ENSANUT 2018
prediabetes_prevalence_2018_year<-svyby(~prediabetes.prev_ifg_hb1ac_2018, by=~YEAR, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_hb1ac_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="High HbA1c-IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2018_year

###ENSANUT 2020
prediabetes_prevalence_2020_year<-svyby(~prediabetes.prev_ifg_hb1ac_2020, by=~YEAR, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_hb1ac_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="High HbA1c-IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2020_year

###ENSANUT 2021
prediabetes_prevalence_2021_year<-svyby(~prediabetes.prev_ifg_hb1ac_2021, by=~YEAR, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_hb1ac_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="High HbA1c-IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2021_year

###ENSANUT 2022
prediabetes_prevalence_2022_year<-svyby(~prediabetes.prev_ifg_hb1ac_2022, by=~YEAR, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_hb1ac_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="High HbA1c-IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2022_year

prev_ifg_hba1c<-rbind(prediabetes_prevalence_2016_year,prediabetes_prevalence_2018_year,prediabetes_prevalence_2020_year,
                      prediabetes_prevalence_2021_year,prediabetes_prevalence_2022_year)
prev_ifg_hba1c<-as.data.frame(prev_ifg_hba1c)
prev_ifg_hba1c$YEAR<-rownames(prev_ifg_hba1c)
prev_ifg_hba1c

#### Discordant HBA1c (prevalence by year) ####
###ENSANUT 2016
prediabetes_prevalence_2016_year<-svyby(~prediabetes.prev_hb1ac_2016, by=~YEAR, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_hb1ac_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="High HbA1c-NFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2016_year

###ENSANUT 2018
prediabetes_prevalence_2018_year<-svyby(~prediabetes.prev_hb1ac_2018, by=~YEAR, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_hb1ac_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="High HbA1c-NFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2018_year

###ENSANUT 2020
prediabetes_prevalence_2020_year<-svyby(~prediabetes.prev_hb1ac_2020, by=~YEAR, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_hb1ac_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="High HbA1c-NFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2020_year

###ENSANUT 2021
prediabetes_prevalence_2021_year<-svyby(~prediabetes.prev_hb1ac_2021, by=~YEAR, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_hb1ac_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="High HbA1c-NFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2021_year

###ENSANUT 2022
prediabetes_prevalence_2022_year<-svyby(~prediabetes.prev_hb1ac_2022, by=~YEAR, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_hb1ac_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="High HbA1c-NFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2022_year

prev_disc_hba1c<-rbind(prediabetes_prevalence_2016_year,prediabetes_prevalence_2018_year,prediabetes_prevalence_2020_year,
                       prediabetes_prevalence_2021_year,prediabetes_prevalence_2022_year)
prev_disc_hba1c<-as.data.frame(prev_disc_hba1c)
prev_disc_hba1c$YEAR<-rownames(prev_disc_hba1c)
prev_disc_hba1c

#### Discordant IFG (prevalence by year) ####
###ENSANUT 2016
prediabetes_prevalence_2016_year<-svyby(~prediabetes.prev_ifg_2016, by=~YEAR, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG-Normal HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2016_year

###ENSANUT 2018
prediabetes_prevalence_2018_year<-svyby(~prediabetes.prev_ifg_2018, by=~YEAR, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG-Normal HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2018_year

###ENSANUT 2020
prediabetes_prevalence_2020_year<-svyby(~prediabetes.prev_ifg_2020, by=~YEAR, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG-Normal HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2020_year

###ENSANUT 2021
prediabetes_prevalence_2021_year<-svyby(~prediabetes.prev_ifg_2021, by=~YEAR, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG-Normal HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2021_year

###ENSANUT 2022
prediabetes_prevalence_2022_year<-svyby(~prediabetes.prev_ifg_2022, by=~YEAR, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG-Normal HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2022_year

prev_disc_ifg<-rbind(prediabetes_prevalence_2016_year,prediabetes_prevalence_2018_year,prediabetes_prevalence_2020_year,
                     prediabetes_prevalence_2021_year,prediabetes_prevalence_2022_year)
prev_disc_ifg<-as.data.frame(prev_disc_ifg)
prev_disc_ifg$YEAR<-rownames(prev_disc_ifg)
prev_disc_ifg

#### All criteria (prevalence by year) ####
###ENSANUT 2016
prediabetes_prevalence_2016_year<-svyby(~prediabetes_allcriteria_2016, by=~YEAR, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="All") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2016_year

###ENSANUT 2018
prediabetes_prevalence_2018_year<-svyby(~prediabetes_allcriteria_2018, by=~YEAR, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="All") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2018_year

###ENSANUT 2020
prediabetes_prevalence_2020_year<-svyby(~prediabetes_allcriteria_2020, by=~YEAR, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="All") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2020_year

###ENSANUT 2021
prediabetes_prevalence_2021_year<-svyby(~prediabetes_allcriteria_2021, by=~YEAR, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="All") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2021_year

###ENSANUT 2022
prediabetes_prevalence_2022_year<-svyby(~prediabetes_allcriteria_2022, by=~YEAR, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="All") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2022_year

prediabetes<-rbind(prediabetes_prevalence_2016_year,prediabetes_prevalence_2018_year,prediabetes_prevalence_2020_year,
                   prediabetes_prevalence_2021_year,prediabetes_prevalence_2022_year)
prediabetes<-as.data.frame(prediabetes)
prediabetes$YEAR<-rownames(prediabetes)
prediabetes

####------------------------------ OTHER ------------------------------#### ----####
#### Previous diabetes (prevalence by year) ####
###ENSANUT 2016
prediabetes_prevalence_2016_year<-svyby(~previous_diabetes_2016, by=~YEAR, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(previous_diabetes_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Diabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2016_year

###ENSANUT 2018
prediabetes_prevalence_2018_year<-svyby(~previous_diabetes_2018, by=~YEAR, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(previous_diabetes_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Diabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2018_year

###ENSANUT 2020
prediabetes_prevalence_2020_year<-svyby(~previous_diabetes_2020, by=~YEAR, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(previous_diabetes_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Diabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2020_year

###ENSANUT 2021
prediabetes_prevalence_2021_year<-svyby(~previous_diabetes_2021, by=~YEAR, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(previous_diabetes_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Diabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2021_year

###ENSANUT 2022
prediabetes_prevalence_2022_year<-svyby(~previous_diabetes_2022, by=~YEAR, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(previous_diabetes_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Diabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2022_year

previous_diabetes<-rbind(prediabetes_prevalence_2016_year,prediabetes_prevalence_2018_year,prediabetes_prevalence_2020_year,
                         prediabetes_prevalence_2021_year,prediabetes_prevalence_2022_year)
previous_diabetes<-as.data.frame(previous_diabetes)
previous_diabetes$YEAR<-rownames(previous_diabetes)
previous_diabetes

#### Total diabetes (prevalence by year) ####
###ENSANUT 2016
prediabetes_prevalence_2016_year<-svyby(~diabetes_fin, by=~YEAR, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(diabetes_fin*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Diabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2016_year

###ENSANUT 2018
prediabetes_prevalence_2018_year<-svyby(~diabetes_fin, by=~YEAR, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(diabetes_fin*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Diabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2018_year

###ENSANUT 2020
prediabetes_prevalence_2020_year<-svyby(~diabetes_fin, by=~YEAR, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(diabetes_fin*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Diabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2020_year

###ENSANUT 2021
prediabetes_prevalence_2021_year<-svyby(~diabetes_fin, by=~YEAR, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(diabetes_fin*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Diabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2021_year

###ENSANUT 2022
prediabetes_prevalence_2022_year<-svyby(~diabetes_fin, by=~YEAR, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(diabetes_fin*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Diabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2022_year

diabetes<-rbind(prediabetes_prevalence_2016_year,prediabetes_prevalence_2018_year,prediabetes_prevalence_2020_year,
                prediabetes_prevalence_2021_year,prediabetes_prevalence_2022_year)
diabetes<-as.data.frame(diabetes)
diabetes$YEAR<-rownames(diabetes)
diabetes

#### IR HOMA2-IR (prevalence by year) ####
###ENSANUT 2016
prediabetes_prevalence_2016_year<-svyby(~homa_cat, by=~YEAR, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(homa_cat*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Insulin Resistance") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2016_year

###ENSANUT 2018
prediabetes_prevalence_2018_year<-svyby(~homa_cat, by=~YEAR, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(homa_cat*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Insulin Resistance") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2018_year

###ENSANUT 2020
prediabetes_prevalence_2020_year<-svyby(~homa_cat, by=~YEAR, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(homa_cat*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Insulin Resistance") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2020_year

###ENSANUT 2021
prediabetes_prevalence_2021_year<-svyby(~homa_cat, by=~YEAR, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(homa_cat*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Insulin Resistance") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2021_year

###ENSANUT 2022
prediabetes_prevalence_2022_year<-svyby(~homa_cat, by=~YEAR, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(homa_cat*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Insulin Resistance") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster)
prediabetes_prevalence_2022_year

homa_cat<-rbind(prediabetes_prevalence_2016_year,prediabetes_prevalence_2018_year,prediabetes_prevalence_2020_year,
                prediabetes_prevalence_2021_year,prediabetes_prevalence_2022_year)
homa_cat<-as.data.frame(homa_cat)
homa_cat$YEAR<-rownames(homa_cat)
homa_cat

#### Healthy (prevalence by year) ####
###ENSANUT 2016
healthy_prevalence_2016_year <- svyby(~any_glucose==0, by=~YEAR, design=ensanut_2016_survey, svymean, na.rm=T)[c(1,3,5)] %>%
  `names<-`(c("YEAR", "ANY", "SE")) %>% mutate(prop=round(ANY*100, digits=1), IC95=round((SE*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="None") %>% dplyr::select(prop,IC95,lIC95,uIC95, cluster)
healthy_prevalence_2016_year

###ENSANUT 2018
healthy_prevalence_2018_year <- svyby(~any_glucose==0, by=~YEAR, design=ensanut_2018_survey, svymean, na.rm=T)[c(1,3,5)] %>%
  `names<-`(c("YEAR", "ANY", "SE")) %>% mutate(prop=round(ANY*100, digits=1), IC95=round((SE*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="None") %>% dplyr::select(prop,IC95,lIC95,uIC95, cluster)
healthy_prevalence_2018_year

###ENSANUT 2020
healthy_prevalence_2020_year <- svyby(~any_glucose==0, by=~YEAR, design=ensanut_2020_survey, svymean, na.rm=T)[c(1,3,5)] %>%
  `names<-`(c("YEAR", "ANY", "SE")) %>% mutate(prop=round(ANY*100, digits=1), IC95=round((SE*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="None") %>% dplyr::select(prop,IC95,lIC95,uIC95, cluster)
healthy_prevalence_2020_year

###ENSANUT 2021
healthy_prevalence_2021_year <- svyby(~any_glucose==0, by=~YEAR, design=ensanut_2021_survey, svymean, na.rm=T)[c(1,3,5)] %>%
  `names<-`(c("YEAR", "ANY", "SE")) %>% mutate(prop=round(ANY*100, digits=1), IC95=round((SE*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="None") %>% dplyr::select(prop,IC95,lIC95,uIC95, cluster)
healthy_prevalence_2021_year

###ENSANUT 2022
healthy_prevalence_2022_year <- svyby(~any_glucose==0, by=~YEAR, design=ensanut_2022_survey, svymean, na.rm=T)[c(1,3,5)] %>%
  `names<-`(c("YEAR", "ANY", "SE")) %>% mutate(prop=round(ANY*100, digits=1), IC95=round((SE*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="None") %>% dplyr::select(prop,IC95,lIC95,uIC95, cluster)
healthy_prevalence_2022_year

healthy<-rbind(healthy_prevalence_2016_year,healthy_prevalence_2018_year,healthy_prevalence_2020_year,
               healthy_prevalence_2021_year,healthy_prevalence_2022_year)
healthy<-as.data.frame(healthy)
healthy$YEAR<-rownames(healthy)
healthy

####------------------------------ STRATA -----------------------------#### ----####
#### Stratification by Age (3) ####
##--- PREDIABETES all criteria ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_allcriteria_2016, by=~YEAR+edad_cat_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_allcriteria_2018, by=~YEAR+edad_cat_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_allcriteria_2020, by=~YEAR+edad_cat_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_allcriteria_2021, by=~YEAR+edad_cat_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_allcriteria_2022, by=~YEAR+edad_cat_2022, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

all_age<-rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
all_age <- all_age %>% as.data.frame %>% mutate("age_cat"= c(rep(c(1,2,3),5)))

##--- PREDIABETES HbA1c ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_hb1ac_2016, by=~YEAR+edad_cat_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_hb1ac_2018, by=~YEAR+edad_cat_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_hb1ac_2020, by=~YEAR+edad_cat_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_hb1ac_2021, by=~YEAR+edad_cat_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_hb1ac_2022, by=~YEAR+edad_cat_2022, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

hba1c_age <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
hba1c_age <- hba1c_age %>% as.data.frame %>% mutate("age_cat"= c(rep(c(1,2,3),5)))

##--- PREDIABETES IFG ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_ifg_2016, by=~YEAR+edad_cat_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_ifg_2018, by=~YEAR+edad_cat_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_ifg_2020, by=~YEAR+edad_cat_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_ifg_2021, by=~YEAR+edad_cat_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_ifg_2022, by=~YEAR+edad_cat_2022, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

ifg_age <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
ifg_age <- ifg_age %>% as.data.frame %>% mutate("age_cat"= c(rep(c(1,2,3),5)))


#### Stratification by Sex (2) ####
##--- PREDIABETES all criteria ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_allcriteria_2016, by=~YEAR+sexo_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_allcriteria_2018, by=~YEAR+sexo_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_allcriteria_2020, by=~YEAR+sexo_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_allcriteria_2021, by=~YEAR+sexo_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_allcriteria_2022, by=~YEAR+sexo_2022, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

all_sex <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
all_sex <- all_sex %>% as.data.frame %>% mutate("sex"= c(rep(c(1,2),5)))

##--- PREDIABETES HbA1c ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_hb1ac_2016, by=~YEAR+sexo_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_hb1ac_2018, by=~YEAR+sexo_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_hb1ac_2020, by=~YEAR+sexo_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_hb1ac_2021, by=~YEAR+sexo_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_hb1ac_2022, by=~YEAR+sexo_2022, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

hba1c_sex <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
hba1c_sex <- hba1c_sex %>% as.data.frame %>% mutate("sex"= c(rep(c(1,2),5)))

##--- PREDIABETES IFG ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_ifg_2016, by=~YEAR+sexo_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_ifg_2018, by=~YEAR+sexo_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_ifg_2020, by=~YEAR+sexo_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_ifg_2021, by=~YEAR+sexo_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_ifg_2022, by=~YEAR+sexo_2022, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

ifg_sex <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
ifg_sex <- ifg_sex %>% as.data.frame %>% mutate("sex"= c(rep(c(1,2),5)))


#### Stratification by BMI (3) ####
##--- PREDIABETES all criteria ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_allcriteria_2016, by=~YEAR+imc_cat_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_allcriteria_2018, by=~YEAR+imc_cat_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_allcriteria_2020, by=~YEAR+imc_cat_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_allcriteria_2021, by=~YEAR+imc_cat_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_allcriteria_2022, by=~YEAR+imc_cat_2022, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

all_imc<-rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
all_imc <- all_imc %>% as.data.frame %>% mutate("imc_cat"= c(rep(c(1,2,3),5)))

##--- PREDIABETES HbA1c ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_hb1ac_2016, by=~YEAR+imc_cat_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_hb1ac_2018, by=~YEAR+imc_cat_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_hb1ac_2020, by=~YEAR+imc_cat_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_hb1ac_2021, by=~YEAR+imc_cat_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_hb1ac_2022, by=~YEAR+imc_cat_2022, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

hba1c_imc<-rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
hba1c_imc <- hba1c_imc %>% as.data.frame %>% mutate("imc_cat"= c(rep(c(1,2,3),5)))

##--- PREDIABETES IFG ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_ifg_2016, by=~YEAR+imc_cat_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_ifg_2018, by=~YEAR+imc_cat_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_ifg_2020, by=~YEAR+imc_cat_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_ifg_2021, by=~YEAR+imc_cat_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_ifg_2022, by=~YEAR+imc_cat_2022, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

ifg_imc<-rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
ifg_imc <- ifg_imc %>% as.data.frame %>% mutate("imc_cat"= c(rep(c(1,2,3),5)))


#### Stratification by WC  (2) ####
##--- PREDIABETES all criteria ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_allcriteria_2016, by=~YEAR+ob_central, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_allcriteria_2018, by=~YEAR+ob_central, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_allcriteria_2021, by=~YEAR+ob_central, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_allcriteria_2022, by=~YEAR+ob_central, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

all_obc <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2021,prediabetes_prev_2022)
all_obc <- all_obc %>% as.data.frame %>% mutate("ob_central"= c(rep(c(1,2),4)))

##--- PREDIABETES HbA1c ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_hb1ac_2016, by=~YEAR+ob_central, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_hb1ac_2018, by=~YEAR+ob_central, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_hb1ac_2021, by=~YEAR+ob_central, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_hb1ac_2022, by=~YEAR+ob_central, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

hba1c_obc <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2021,prediabetes_prev_2022)
hba1c_obc <- hba1c_obc %>% as.data.frame %>% mutate("ob_central"= c(rep(c(1,2),4)))

##--- PREDIABETES IFG ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_ifg_2016, by=~YEAR+ob_central, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_ifg_2018, by=~YEAR+ob_central, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_ifg_2021, by=~YEAR+ob_central, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_ifg_2022, by=~YEAR+ob_central, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

ifg_obc <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2021,prediabetes_prev_2022)
ifg_obc <- ifg_obc %>% as.data.frame %>% mutate("ob_central"= c(rep(c(1,2),4)))


#### Stratification by T2D Family HX    (2) ####
##--- PREDIABETES all criteria ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_allcriteria_2016, by=~YEAR+T2D.AHF, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_allcriteria_2018, by=~YEAR+T2D.AHF, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_allcriteria_2021, by=~YEAR+T2D.AHF, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_allcriteria_2022, by=~YEAR+T2D.AHF, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

all_ahf <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2021,prediabetes_prev_2022)
all_ahf <- all_ahf %>% as.data.frame %>% mutate("T2D.AHF"= c(rep(c(1,2),4)))

##--- PREDIABETES HbA1c ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_hb1ac_2016, by=~YEAR+T2D.AHF, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_hb1ac_2018, by=~YEAR+T2D.AHF, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_hb1ac_2021, by=~YEAR+T2D.AHF, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_hb1ac_2022, by=~YEAR+T2D.AHF, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

hba1c_ahf <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2021,prediabetes_prev_2022)
hba1c_ahf <- hba1c_ahf %>% as.data.frame %>% mutate("T2D.AHF"= c(rep(c(1,2),4)))

##--- PREDIABETES IFG ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_ifg_2016, by=~YEAR+T2D.AHF, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_ifg_2018, by=~YEAR+T2D.AHF, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_ifg_2021, by=~YEAR+T2D.AHF, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_ifg_2022, by=~YEAR+T2D.AHF, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

ifg_ahf <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2021,prediabetes_prev_2022)
ifg_ahf <- ifg_ahf %>% as.data.frame %>% mutate("T2D.AHF"= c(rep(c(1,2),4)))


#### Stratification by Smoking Status   (3) ####
##--- PREDIABETES all criteria ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_allcriteria_2016, by=~YEAR+smoking, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_allcriteria_2018, by=~YEAR+smoking, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_allcriteria_2021, by=~YEAR+smoking, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_allcriteria_2022, by=~YEAR+smoking, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

all_smo <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2021,prediabetes_prev_2022)
all_smo <- all_smo %>% as.data.frame %>% mutate("smoking"= c(rep(c(1,2,3),4)))

##--- PREDIABETES HbA1c ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_hb1ac_2016, by=~YEAR+smoking, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_hb1ac_2018, by=~YEAR+smoking, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_hb1ac_2021, by=~YEAR+smoking, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_hb1ac_2022, by=~YEAR+smoking, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

hba1c_smo <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2021,prediabetes_prev_2022)
hba1c_smo <- hba1c_smo %>% as.data.frame %>% mutate("smoking"= c(rep(c(1,2,3),4)))

##--- PREDIABETES IFG ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_ifg_2016, by=~YEAR+smoking, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_ifg_2018, by=~YEAR+smoking, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_ifg_2021, by=~YEAR+smoking, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_ifg_2022, by=~YEAR+smoking, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

ifg_smo <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2021,prediabetes_prev_2022)
ifg_smo <- ifg_smo %>% as.data.frame %>% mutate("smoking"= c(rep(c(1,2,3),4)))


#### Stratification by Indigenous Lang. (2) ####
##--- PREDIABETES all criteria ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_allcriteria_2016, by=~YEAR+lengua_indigena, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_allcriteria_2018, by=~YEAR+lengua_indigena, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_allcriteria_2020, by=~YEAR+lengua_indigena, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_allcriteria_2021, by=~YEAR+lengua_indigena, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_allcriteria_2022, by=~YEAR+lengua_indigena, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

all_inl <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
all_inl <- all_inl %>% as.data.frame %>% mutate("inl"= c(rep(c(1,2),5)))

##--- PREDIABETES HbA1c ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_hb1ac_2016, by=~YEAR+lengua_indigena, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_hb1ac_2018, by=~YEAR+lengua_indigena, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_hb1ac_2020, by=~YEAR+lengua_indigena, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_hb1ac_2021, by=~YEAR+lengua_indigena, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_hb1ac_2022, by=~YEAR+lengua_indigena, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

hba1c_inl <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
hba1c_inl <- hba1c_inl %>% as.data.frame %>% mutate("inl"= c(rep(c(1,2),5)))

##--- PREDIABETES IFG ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_ifg_2016, by=~YEAR+lengua_indigena, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_ifg_2018, by=~YEAR+lengua_indigena, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_ifg_2020, by=~YEAR+lengua_indigena, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_ifg_2021, by=~YEAR+lengua_indigena, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_ifg_2022, by=~YEAR+lengua_indigena, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

ifg_inl <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
ifg_inl <- ifg_sex %>% as.data.frame %>% mutate("inl"= c(rep(c(1,2),5)))


#### Stratification by Density Ind. SLI (2) ####
##--- PREDIABETES all criteria ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_allcriteria_2016, by=~YEAR+DISLI_cat2, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_allcriteria_2018, by=~YEAR+DISLI_cat2, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_allcriteria_2020, by=~YEAR+DISLI_cat2, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_allcriteria_2021, by=~YEAR+DISLI_cat2, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_allcriteria_2022, by=~YEAR+DISLI_cat2, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

all_sli <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
all_sli <- all_sli %>% as.data.frame %>% mutate("SLI"= c(rep(c(1,2),5)))

##--- PREDIABETES HbA1c ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_hb1ac_2016, by=~YEAR+DISLI_cat2, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_hb1ac_2018, by=~YEAR+DISLI_cat2, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_hb1ac_2020, by=~YEAR+DISLI_cat2, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_hb1ac_2021, by=~YEAR+DISLI_cat2, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_hb1ac_2022, by=~YEAR+DISLI_cat2, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

hba1c_sli <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
hba1c_sli <- hba1c_sli %>% as.data.frame %>% mutate("SLI"= c(rep(c(1,2),5)))

##--- PREDIABETES IFG ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_ifg_2016, by=~YEAR+DISLI_cat2, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_ifg_2018, by=~YEAR+DISLI_cat2, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_ifg_2020, by=~YEAR+DISLI_cat2, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_ifg_2021, by=~YEAR+DISLI_cat2, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_ifg_2022, by=~YEAR+DISLI_cat2, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

ifg_sli <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
ifg_sli <- ifg_sli %>% as.data.frame %>% mutate("SLI"= c(rep(c(1,2),5)))


#### Stratification by Area Rural/Urban (2) ####
##--- PREDIABETES all criteria ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_allcriteria_2016, by=~YEAR+area2_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_allcriteria_2018, by=~YEAR+area_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_allcriteria_2020, by=~YEAR+area_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_allcriteria_2021, by=~YEAR+area2_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_allcriteria_2022, by=~YEAR+area2_2022, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

all_aru <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
all_aru <- all_aru %>% as.data.frame %>% mutate("area"= c(rep(c(1,2),5)))

##--- PREDIABETES HbA1c ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_hb1ac_2016, by=~YEAR+area2_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_hb1ac_2018, by=~YEAR+area_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_hb1ac_2020, by=~YEAR+area_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_hb1ac_2021, by=~YEAR+area2_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_hb1ac_2022, by=~YEAR+area2_2022, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

hba1c_aru <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
hba1c_aru <- hba1c_aru %>% as.data.frame %>% mutate("area"= c(rep(c(1,2),5)))

##--- PREDIABETES IFG ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_ifg_2016, by=~YEAR+area2_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_ifg_2018, by=~YEAR+area_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_ifg_2020, by=~YEAR+area_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_ifg_2021, by=~YEAR+area2_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~prediabetes_ifg_2022, by=~YEAR+area2_2022, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2022*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)

ifg_aru <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021,prediabetes_prev_2022)
ifg_aru <- ifg_aru %>% as.data.frame %>% mutate("area"= c(rep(c(1,2),5)))


#### Stratification by Region (4) ####
##--- PREDIABETES all criteria ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_allcriteria_2016, by=~YEAR+regiones_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_allcriteria_2018, by=~YEAR+regiones_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_allcriteria_2020, by=~YEAR+regiones_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_allcriteria_2021, by=~YEAR+regiones_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
all_reg <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021)
all_reg <- all_reg %>% as.data.frame %>% mutate("region"= c(rep(c(1,2,3,4),4)))

##--- PREDIABETES HbA1c ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_hb1ac_2016, by=~YEAR+regiones_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_hb1ac_2018, by=~YEAR+regiones_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_hb1ac_2020, by=~YEAR+regiones_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_hb1ac_2021, by=~YEAR+regiones_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
a1c_reg <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021)
a1c_reg <- a1c_reg %>% as.data.frame %>% mutate("region"= c(rep(c(1,2,3,4),4)))

##--- PREDIABETES IFG ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_ifg_2016, by=~YEAR+regiones_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_ifg_2018, by=~YEAR+regiones_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_ifg_2020, by=~YEAR+regiones_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_ifg_2021, by=~YEAR+regiones_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
ifg_reg <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021)
ifg_reg <- ifg_reg %>% as.data.frame %>% mutate("region"= c(rep(c(1,2,3,4),4)))

##--- PREDIABETES BOTH ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes.prev_ifg_hb1ac_2016, by=~YEAR+regiones_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_hb1ac_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Both") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes.prev_ifg_hb1ac_2018, by=~YEAR+regiones_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_hb1ac_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Both") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes.prev_ifg_hb1ac_2020, by=~YEAR+regiones_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_hb1ac_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Both") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes.prev_ifg_hb1ac_2021, by=~YEAR+regiones_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_hb1ac_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Both") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
bth_reg <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021)
bth_reg <- bth_reg %>% as.data.frame %>% mutate("region"= c(rep(c(1,2,3,4),4)))

#### Stratification by Region (9) ####
##--- PREDIABETES all criteria ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_allcriteria_2016, by=~YEAR+regiones2_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_allcriteria_2018, by=~YEAR+regiones2_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_allcriteria_2020, by=~YEAR+regiones2_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_allcriteria_2021, by=~YEAR+regiones2_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_allcriteria_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Prediabetes") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
all_reg9 <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021)
all_reg9 <- all_reg9 %>% as.data.frame %>% mutate("region2"= c(rep(c(1:9),4)))

##--- PREDIABETES HbA1c ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_hb1ac_2016, by=~YEAR+regiones2_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_hb1ac_2018, by=~YEAR+regiones2_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_hb1ac_2020, by=~YEAR+regiones2_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_hb1ac_2021, by=~YEAR+regiones2_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_hb1ac_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
a1c_reg9 <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021)
a1c_reg9 <- a1c_reg9 %>% as.data.frame %>% mutate("region2"= c(rep(c(1:9),4)))

##--- PREDIABETES IFG ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes_ifg_2016, by=~YEAR+regiones2_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes_ifg_2018, by=~YEAR+regiones2_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes_ifg_2020, by=~YEAR+regiones2_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes_ifg_2021, by=~YEAR+regiones2_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes_ifg_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
ifg_reg9 <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021)
ifg_reg9 <- ifg_reg9 %>% as.data.frame %>% mutate("region2"= c(rep(c(1:9),4)))

##--- PREDIABETES BOTH ---##
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~prediabetes.prev_ifg_hb1ac_2016, by=~YEAR+regiones2_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_hb1ac_2016*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Both") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~prediabetes.prev_ifg_hb1ac_2018, by=~YEAR+regiones2_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_hb1ac_2018*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Both") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~prediabetes.prev_ifg_hb1ac_2020, by=~YEAR+regiones2_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_hb1ac_2020*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Both") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~prediabetes.prev_ifg_hb1ac_2021, by=~YEAR+regiones2_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(prediabetes.prev_ifg_hb1ac_2021*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Both") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
bth_reg9 <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,prediabetes_prev_2021)
bth_reg9 <- bth_reg9 %>% as.data.frame %>% mutate("region2"= c(rep(c(1:9),4)))

####---------------------------- REGRESSIONS ---------------------------#### ----####
#### Regression Outcomes ####
## --- # --- # METABOLIC SYNDROME COMPONENTS # --- # --- ##
##-- Arterial hypertension (DX | TX | ≥140/90) --##
#2016
ensanut_2016$PAS_fin <- ensanut_2016 %>% select(sistol3.y, sistol4.y) %>% apply(1, mean, na.rm=T)
ensanut_2016$PAD_fin <- ensanut_2016 %>% select(diastol3.y, diastol4.y) %>% apply(1, mean, na.rm=T)
ensanut_2016$HTA_fin <- (with(ensanut_2016, (
  a401.x==1 | (a405.y==1&!is.na(a405.y)) | PAS_fin>=140&!is.na(PAS_fin) |
    PAD_fin>=90&!is.na(PAD_fin) ))) %>% ifelse(1,0); ensanut_2016$HTA_fin %>% table(useNA = "always")
#2018
ensanut_2018$P27_1_1[ensanut_2018$P27_1_1<80]<-NA; ensanut_2018$P27_2_1[ensanut_2018$P27_2_1<80]<-NA
ensanut_2018$P27_1_2[ensanut_2018$P27_1_2<49]<-NA; ensanut_2018$P27_2_2[ensanut_2018$P27_2_2<49]<-NA
ensanut_2018$P27_1_2[ensanut_2018$P27_1_2==222]<-NA
ensanut_2018$PAS_fin <- ensanut_2018 %>% select(P27_1_1, P27_2_1) %>% apply(1, mean, na.rm=T)
ensanut_2018$PAD_fin <- ensanut_2018 %>% select(P27_1_2, P27_2_2) %>% apply(1, mean, na.rm=T)
ensanut_2018$HTA_fin <- (with(ensanut_2018, (
  P4_1.x==1 | (P4_4==1&!is.na(P4_4)) | (PAS_fin>=140&!is.na(PAS_fin)) |
    (PAD_fin>=90&!is.na(PAD_fin)) ))) %>% ifelse(1,0); ensanut_2018$HTA_fin %>% table(useNA = "always")
#2020 (NO INFO ON TX)
ensanut_2020<-ensanut_2020%>%mutate("HX_HBP"=ifelse((H0902A.x=="3"|H0902B=="3"|H0902C=="3"|H0902D=="3")==T,1,0))
ensanut_2020$an08_01s[ensanut_2020$an08_01s==999]<-NA; ensanut_2020$an08_01d[ensanut_2020$an08_01d==999]<-NA
ensanut_2020$an08_02s[ensanut_2020$an08_02s==999]<-NA; ensanut_2020$an08_02d[ensanut_2020$an08_02d==999]<-NA
ensanut_2020$an08_03s[ensanut_2020$an08_03s==999]<-NA; ensanut_2020$an08_03d[ensanut_2020$an08_03d==999]<-NA
ensanut_2020$PAS_fin <- ensanut_2020 %>% select(an08_01s, an08_02s, an08_03s) %>% apply(1, mean, na.rm=T)
ensanut_2020$PAD_fin <- ensanut_2020 %>% select(an08_01d, an08_02d, an08_03d) %>% apply(1, mean, na.rm=T)
ensanut_2020$HTA_fin <- (with(ensanut_2020, (HX_HBP==1 | PAS_fin>=140&!is.na(PAS_fin) | PAD_fin>=90&!is.na(PAD_fin)))) %>%
  ifelse(1,0); ensanut_2020$HTA_fin %>% table(useNA = "always")
#2021
ensanut_2021 <- ensanut_2021 %>% mutate("a0401.x"=case_when(a0401%in%2:3~0, a0401==1~1))
ensanut_2021$an27_01s[ensanut_2021$an27_01s==999]<-NA; ensanut_2021$an27_01d[ensanut_2021$an27_01d==999]<-NA
ensanut_2021$an27_02s[ensanut_2021$an27_02s==999]<-NA; ensanut_2021$an27_02d[ensanut_2021$an27_02d==999]<-NA
ensanut_2021$an27_03s[ensanut_2021$an27_03s==999]<-NA; ensanut_2021$an27_03d[ensanut_2021$an27_03d==999]<-NA
ensanut_2021$PAS_fin <- ensanut_2021 %>% select(an27_01s, an27_02s, an27_03s) %>% apply(1, mean, na.rm=T)
ensanut_2021$PAD_fin <- ensanut_2021 %>% select(an27_01d, an27_02d, an27_03d) %>% apply(1, mean, na.rm=T)
ensanut_2021$HTA_fin <- (with(ensanut_2021, (
  a0401.x==1 | (a0404==1&!is.na(a0404)) | PAS_fin>=140&!is.na(PAS_fin) |
    PAD_fin>=90&!is.na(PAD_fin) ))) %>% ifelse(1,0); ensanut_2021$HTA_fin %>% table(useNA = "always")
#2022
ensanut_2022 <- ensanut_2022 %>% mutate("a0401.x"=case_when(a0401%in%2:3~0, a0401==1~1))
ensanut_2022$an27_01s[ensanut_2022$an27_01s==999]<-NA; ensanut_2022$an27_01d[ensanut_2022$an27_01d==999]<-NA
ensanut_2022$an27_02s[ensanut_2022$an27_02s==999]<-NA; ensanut_2022$an27_02d[ensanut_2022$an27_02d==999]<-NA
ensanut_2022$an27_03s[ensanut_2022$an27_03s==999]<-NA; ensanut_2022$an27_03d[ensanut_2022$an27_03d==999]<-NA
ensanut_2022$PAS_fin <- ensanut_2022 %>% select(an27_01s, an27_02s, an27_03s) %>% apply(1, mean, na.rm=T)
ensanut_2022$PAD_fin <- ensanut_2022 %>% select(an27_01d, an27_02d, an27_03d) %>% apply(1, mean, na.rm=T)
ensanut_2022$HTA_fin <- (with(ensanut_2022, (
  a0401.x==1 | (a0404==1&!is.na(a0404)) | PAS_fin>=140&!is.na(PAS_fin) |
    PAD_fin>=90&!is.na(PAD_fin) ))) %>% ifelse(1,0); ensanut_2022$HTA_fin %>% table(useNA = "always")

##-- Hypertriglyceridemia (TG ≥150 or treatment) --##
#2016
ensanut_2016$S.TG_fin <- ensanut_2016$valor.TRIG; ensanut_2016$S.TG_fin[ensanut_2016$sanvenh<8]<-NA
ensanut_2016$HTG_biochem <- (with(ensanut_2016, ((S.TG_fin>=150)))) %>% ifelse(1,0)
ensanut_2016$HTG_fin <- ( 
  (ensanut_2016 %>% transmute(a610a==1, S.TG_fin>=150) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2016$HTG_fin[is.na(ensanut_2016$a610a) & is.na(ensanut_2016$S.TG_fin)] <- NA
ensanut_2016$HTG_biochem %>% table(useNA = "always"); ensanut_2016$HTG_fin %>% table(useNA = "always")
#2018
ensanut_2018$S.TG_fin <- ensanut_2018$VALOR_TRIG; ensanut_2018$S.TG_fin[ensanut_2018$P5_1.y<8]<-NA
ensanut_2018$HTG_biochem <- (with(ensanut_2018, ((S.TG_fin>=150)))) %>% ifelse(1,0)
ensanut_2018$HTG_fin <- ( 
  (ensanut_2018 %>% transmute(P6_7_1==1, S.TG_fin>=150) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2018$HTG_fin[is.na(ensanut_2018$P6_7_1) & is.na(ensanut_2018$S.TG_fin)] <- NA
ensanut_2018$HTG_biochem %>% table(useNA = "always"); ensanut_2018$HTG_fin %>% table(useNA = "always")
#2020
ensanut_2020$S.TG_fin <- ensanut_2020$valor.TRIG; ensanut_2020$S.TG_fin[ensanut_2020$san04<8]<-NA
ensanut_2020$HTG_biochem <- (with(ensanut_2020, ((S.TG_fin>=150)))) %>% ifelse(1,0)
ensanut_2020$HTG_fin <- ensanut_2020$HTG_biochem
ensanut_2020$HTG_biochem %>% table(useNA = "always"); ensanut_2020$HTG_fin %>% table(useNA = "always")
#2021
ensanut_2021$S.TG_fin <- ensanut_2021$valor_TRIG; ensanut_2021$S.TG_fin[ensanut_2021$san04<8]<-NA
ensanut_2021$HTG_biochem <- (with(ensanut_2021, ((S.TG_fin>=150)))) %>% ifelse(1,0)
ensanut_2021$HTG_fin <- ( 
  (ensanut_2021 %>% transmute(A0607A==1, S.TG_fin>=150) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2021$HTG_fin[is.na(ensanut_2021$A0607A) & is.na(ensanut_2021$S.TG_fin)] <- NA
ensanut_2021$HTG_biochem %>% table(useNA = "always"); ensanut_2021$HTG_fin %>% table(useNA = "always")
#2022
ensanut_2022$S.TG_fin <- ensanut_2022$valor_TRIG; ensanut_2022$S.TG_fin[ensanut_2022$san04<8]<-NA
ensanut_2022$HTG_biochem <- (with(ensanut_2022, ((S.TG_fin>=150)))) %>% ifelse(1,0)
ensanut_2022$HTG_fin <- ( 
  (ensanut_2022 %>% transmute(A0607A1==1, S.TG_fin>=150) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2022$HTG_fin[is.na(ensanut_2022$A0607A1) & is.na(ensanut_2022$S.TG_fin)] <- NA
ensanut_2022$HTG_biochem %>% table(useNA = "always"); ensanut_2022$HTG_fin %>% table(useNA = "always")

##-- Hypercholesterolemia (TC ≥200 or treatment) --##
#2016 
ensanut_2016$S.TC_fin <- ensanut_2016$valor.COLEST; ensanut_2016$S.TC_fin[ensanut_2016$sanvenh<8]<-NA
ensanut_2016$HCT_biochem <- (with(ensanut_2016, ((S.TC_fin>=200)))) %>% ifelse(1,0)
ensanut_2016$HCT_fin <- ( 
  (ensanut_2016 %>% transmute(a608a==1, S.TC_fin>=200) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2016$HCT_fin[is.na(ensanut_2016$a608a) & is.na(ensanut_2016$S.TG_fin)] <- NA
ensanut_2016$HCT_biochem %>% table(useNA = "always"); ensanut_2016$HCT_fin %>% table(useNA = "always")
#2018
ensanut_2018$S.TC_fin <- ensanut_2018$VALOR_COLEST; ensanut_2018$S.TC_fin[ensanut_2018$P5_1.y<8]<-NA
ensanut_2018$HCT_biochem <- (with(ensanut_2018, ((S.TC_fin>=200)))) %>% ifelse(1,0)
ensanut_2018$HCT_fin <- ( 
  (ensanut_2018 %>% transmute(P6_5_1==1, S.TC_fin>=200) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2018$HCT_fin[is.na(ensanut_2018$P6_5_1) & is.na(ensanut_2018$S.TG_fin)] <- NA
ensanut_2018$HCT_biochem %>% table(useNA = "always"); ensanut_2018$HCT_fin %>% table(useNA = "always")
#2020
ensanut_2020$S.TC_fin <- ensanut_2020$valor.COLEST; ensanut_2020$S.TC_fin[ensanut_2020$san04<8]<-NA
ensanut_2020$HCT_biochem <- (with(ensanut_2020, ((S.TC_fin>=200)))) %>% ifelse(1,0)
ensanut_2020$HCT_fin <- ensanut_2020$HCT_biochem
ensanut_2020$HCT_biochem %>% table(useNA = "always"); ensanut_2020$HCT_fin %>% table(useNA = "always")
#2021
ensanut_2021$S.TC_fin <- ensanut_2021$valor_COLEST; ensanut_2021$S.TC_fin[ensanut_2021$san04<8]<-NA
ensanut_2021$HCT_biochem <- (with(ensanut_2021, ((S.TC_fin>=200)))) %>% ifelse(1,0)
ensanut_2021$HCT_fin <- ( 
  (ensanut_2021 %>% transmute(A0605A==1, S.TC_fin>=200) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2021$HCT_fin[is.na(ensanut_2021$A0605A) & is.na(ensanut_2021$S.TG_fin)] <- NA
ensanut_2021$HCT_biochem %>% table(useNA = "always"); ensanut_2021$HCT_fin %>% table(useNA = "always")
#2022
ensanut_2022$S.TC_fin <- ensanut_2022$valor_COLEST; ensanut_2022$S.TC_fin[ensanut_2022$san04<8]<-NA
ensanut_2022$HCT_biochem <- (with(ensanut_2022, ((S.TC_fin>=200)))) %>% ifelse(1,0)
ensanut_2022$HCT_fin <- ( 
  (ensanut_2022 %>% transmute(A0605A==1, S.TC_fin>=200) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2022$HCT_fin[is.na(ensanut_2022$A0605A) & is.na(ensanut_2022$S.TG_fin)] <- NA
ensanut_2022$HCT_biochem %>% table(useNA = "always"); ensanut_2022$HCT_fin %>% table(useNA = "always")

##-- HDL-C fasting correction --##
ensanut_2016$S.HDL_fin <- ensanut_2016$valor.COL_HDL; ensanut_2016$S.HDL_fin[ensanut_2016$sanvenh<8]<-NA
ensanut_2018$S.HDL_fin <- ensanut_2018$VALOR_COLEST; ensanut_2018$S.HDL_fin[ensanut_2018$P5_1.y<8]<-NA
ensanut_2020$S.HDL_fin <- ensanut_2020$valor.COLEST; ensanut_2020$S.HDL_fin[ensanut_2020$san04<8]<-NA
ensanut_2021$S.HDL_fin <- ensanut_2021$valor_COLEST; ensanut_2021$S.HDL_fin[ensanut_2021$san04<8]<-NA
ensanut_2022$S.HDL_fin <- ensanut_2022$valor_COLEST; ensanut_2022$S.HDL_fin[ensanut_2022$san04<8]<-NA

##-- Lipid lowering medications --##
ensanut_2016$TG_TX <- ensanut_2016$a610a; ensanut_2016$TC_TX <- ensanut_2016$a608a
ensanut_2018$TG_TX <- ensanut_2018$P6_7_1; ensanut_2018$TC_TX <- ensanut_2018$P6_5_1
ensanut_2021$TG_TX <- ensanut_2021$A0607A; ensanut_2021$TC_TX <- ensanut_2021$A0605A
ensanut_2022$TG_TX <- ensanut_2022$A0607A1; ensanut_2022$TC_TX <- ensanut_2022$A0605A

##-- Metabolic syndrome (IDF) --##
# Obesity + 2 of the following: High TG (or TX), Low HDL (or TX), High BP (or TX), High FPG
# We excluded High FPG due to this being the predictor -> thus: Obesity + 1 of the following
#2016
ensanut_2016$obesity <- ( #Central obesity (WC) or BMI ≥30
  (ensanut_2016 %>% transmute(ob_central, imc>=30) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2016$obesity[is.na(ensanut_2016$ob_central) & is.na(ensanut_2016$imc)] <- NA
ensanut_2016$highTG <- ( #TG ≥150 or treatment
  (ensanut_2016 %>% transmute(a610a==1, S.TG_fin>=150) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2016$highTG[is.na(ensanut_2016$a610a) & is.na(ensanut_2016$S.TG_fin)] <- NA
ensanut_2016$lowHDL <- ((ensanut_2016 %>% transmute( #HDL <50 (women)  or HDL <40 (men) or treatment
  a608a==1, ((sexo_2016==1 & S.HDL_fin<50) | (sexo_2016==0 & S.HDL_fin<40))) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2016$lowHDL[is.na(ensanut_2016$a608a) & is.na(ensanut_2016$S.HDL_fin)] <- NA
ensanut_2016$highBP <- ( #BP ≥130/85 or treatment
  (ensanut_2016 %>% transmute(a405.y==1, PAS_fin>=130, PAD_fin>=85) %>% apply(1, sum, na.rm=T))>0) %>% ifelse(
    1,0); ensanut_2016$highBP[is.na(ensanut_2016$a405.y) & is.na(ensanut_2016$PAS_fin) & is.na(ensanut_2016$PAD_fin)] <- NA
#MetS
ensanut_2016$MetS <- ((ensanut_2016 %>% select(highTG, lowHDL, highBP) %>% apply(1, sum, na.rm=T))>0) %>% ifelse(1,0) #Metabolic Syndrome (IDF)
ensanut_2016$MetS[with(ensanut_2016, (is.na(highTG) & is.na(lowHDL) & is.na(highBP)))] <- NA; ensanut_2016$MetS[ensanut_2016$obesity==0]<-0
ensanut_2016$MetS %>% table(useNA = "always")
ensanut_2016$MetS_n <- ensanut_2016 %>% select(obesity, highTG, lowHDL, highBP) %>% apply(1, sum, na.rm=T) #Number of criteria fulfilled
ensanut_2016$MetS_n[with(ensanut_2016, (is.na(highTG) & is.na(lowHDL) & is.na(highBP) & is.na(obesity)))] <- NA

#2018
ensanut_2018$obesity <- ( #Central obesity (WC) or BMI ≥30
  (ensanut_2018 %>% transmute(ob_central, imc_calculo_2018>=30) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2018$obesity[is.na(ensanut_2018$ob_central) & is.na(ensanut_2018$imc_calculo_2018)] <- NA
ensanut_2018$highTG <- ( #TG ≥150 or treatment
  (ensanut_2018 %>% transmute(P6_7_1==1, S.TG_fin>=150) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2018$highTG[is.na(ensanut_2018$P6_7_1) & is.na(ensanut_2018$S.TG_fin)] <- NA
ensanut_2018$lowHDL <- ((ensanut_2018 %>% transmute( #HDL <50 (women)  or HDL <40 (men) or treatment
  P6_5_1==1, ((sexo_2018==1 & S.HDL_fin<50) | (sexo_2018==0 & S.HDL_fin<40))) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2018$lowHDL[is.na(ensanut_2018$P6_5_1) & is.na(ensanut_2018$S.HDL_fin)] <- NA
ensanut_2018$highBP <- ( #BP ≥130/85 or treatment
  (ensanut_2018 %>% transmute(P4_4==1, PAS_fin>=130, PAD_fin>=85) %>% apply(1, sum, na.rm=T))>0) %>% ifelse(
    1,0); ensanut_2018$highBP[is.na(ensanut_2018$P4_4) & is.na(ensanut_2018$PAS_fin) & is.na(ensanut_2018$PAD_fin)] <- NA
#MetS
ensanut_2018$MetS <- ((ensanut_2018 %>% select(highTG, lowHDL, highBP) %>% apply(1, sum, na.rm=T))>0) %>% ifelse(1,0) #Metabolic Syndrome (IDF)
ensanut_2018$MetS[with(ensanut_2018, (is.na(highTG) & is.na(lowHDL) & is.na(highBP)))] <- NA; ensanut_2018$MetS[ensanut_2018$obesity==0]<-0
ensanut_2018$MetS %>% table(useNA = "always")
ensanut_2018$MetS_n <- ensanut_2018 %>% select(obesity, highTG, lowHDL, highBP) %>% apply(1, sum, na.rm=T) #Number of criteria fulfilled
ensanut_2018$MetS_n[with(ensanut_2018, (is.na(highTG) & is.na(lowHDL) & is.na(highBP) & is.na(obesity)))] <- NA

#2020
ensanut_2020$obesity <- ifelse(ensanut_2020$imc_calculo_2020>=30, 1,0) #Obesity (only BMI criteria)
ensanut_2020$highTG <- case_when(ensanut_2020$S.TG_fin>=150~1, ensanut_2020$S.TG_fin<150~0) #TG ≥150
ensanut_2020$lowHDL <- ifelse(with( #HDL <50 (women)  or HDL <40 (men)
  ensanut_2020, (sexo_2020==1 & S.HDL_fin<50) | (sexo_2020==0 & S.HDL_fin<40)),
  1,0); ensanut_2020$lowHDL[with(ensanut_2020, is.na(S.HDL_fin))] <- NA 
ensanut_2020$highBP <- ((ensanut_2020 %>% transmute(PAS_fin>=130, PAD_fin>=85) %>% apply(1, sum, na.rm=T))>0) %>% #BP ≥130/85
  ifelse(1,0); ensanut_2020$highBP[is.na(ensanut_2020$PAS_fin) & is.na(ensanut_2020$PAD_fin)] <- NA
#MetS
ensanut_2020$MetS <- ((ensanut_2020 %>% select(highTG, lowHDL, highBP) %>% apply(1, sum, na.rm=T))>0) %>% ifelse(1,0) #Metabolic Syndrome (IDF)
ensanut_2020$MetS[with(ensanut_2020, (is.na(highTG) & is.na(lowHDL) & is.na(highBP)))] <- NA; ensanut_2020$MetS[ensanut_2020$obesity==0]<-0
ensanut_2020$MetS %>% table(useNA = "always")
ensanut_2020$MetS_n <- ensanut_2020 %>% select(obesity, highTG, lowHDL, highBP) %>% apply(1, sum, na.rm=T) #Number of criteria fulfilled
ensanut_2020$MetS_n[with(ensanut_2020, (is.na(highTG) & is.na(lowHDL) & is.na(highBP) & is.na(obesity)))] <- NA

#2021
ensanut_2021$obesity <- ( #Central obesity (WC) or BMI ≥30
  (ensanut_2021 %>% transmute(ob_central, imc_calculo_2021>=30) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2021$obesity[is.na(ensanut_2021$ob_central) & is.na(ensanut_2021$imc_calculo_2021)] <- NA
ensanut_2021$highTG <- ( #TG ≥150 or treatment
  (ensanut_2021 %>% transmute(A0607A==1, S.TG_fin>=150) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2021$highTG[is.na(ensanut_2021$A0607A) & is.na(ensanut_2021$S.TG_fin)] <- NA
ensanut_2021$lowHDL <- ((ensanut_2021 %>% transmute( #HDL <50 (women)  or HDL <40 (men) or treatment
  A0605A==1, ((sexo_2021==1 & S.HDL_fin<50) | (sexo_2021==0 & S.HDL_fin<40))) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2021$lowHDL[is.na(ensanut_2021$A0605A) & is.na(ensanut_2021$S.HDL_fin)] <- NA
ensanut_2021$highBP <- ( #BP ≥130/85 or treatment
  (ensanut_2021 %>% transmute(a0404==1, PAS_fin>=130, PAD_fin>=85) %>% apply(1, sum, na.rm=T))>0) %>% ifelse(
    1,0); ensanut_2021$highBP[is.na(ensanut_2021$a0404) & is.na(ensanut_2021$PAS_fin) & is.na(ensanut_2021$PAD_fin)] <- NA
#MetS
ensanut_2021$MetS <- ((ensanut_2021 %>% select(highTG, lowHDL, highBP) %>% apply(1, sum, na.rm=T))>0) %>% ifelse(1,0) #Metabolic Syndrome (IDF)
ensanut_2021$MetS[with(ensanut_2021, (is.na(highTG) & is.na(lowHDL) & is.na(highBP)))] <- NA; ensanut_2021$MetS[ensanut_2021$obesity==0]<-0
ensanut_2021$MetS %>% table(useNA = "always")
ensanut_2021$MetS_n <- ensanut_2021 %>% select(obesity, highTG, lowHDL, highBP) %>% apply(1, sum, na.rm=T) #Number of criteria fulfilled
ensanut_2021$MetS_n[with(ensanut_2021, (is.na(highTG) & is.na(lowHDL) & is.na(highBP) & is.na(obesity)))] <- NA

#2022
ensanut_2022$obesity <- ( #Central obesity (WC) or BMI ≥30
  (ensanut_2022 %>% transmute(ob_central, imc_calculo_2022>=30) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2022$obesity[is.na(ensanut_2022$ob_central) & is.na(ensanut_2022$imc_calculo_2022)] <- NA
ensanut_2022$highTG <- ( #TG ≥150 or treatment
  (ensanut_2022 %>% transmute(A0607A1==1, S.TG_fin>=150) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2022$highTG[is.na(ensanut_2022$A0607A1) & is.na(ensanut_2022$S.TG_fin)] <- NA
ensanut_2022$lowHDL <- ((ensanut_2022 %>% transmute( #HDL <50 (women)  or HDL <40 (men) or treatment
  A0605A==1, ((sexo_2022==1 & S.HDL_fin<50) | (sexo_2022==0 & S.HDL_fin<40))) %>% apply(1, sum, na.rm=T))>0) %>%
  ifelse(1,0); ensanut_2022$lowHDL[is.na(ensanut_2022$A0605A) & is.na(ensanut_2022$S.HDL_fin)] <- NA
ensanut_2022$highBP <- ( #BP ≥130/85 or treatment
  (ensanut_2022 %>% transmute(a0404==1, PAS_fin>=130, PAD_fin>=85) %>% apply(1, sum, na.rm=T))>0) %>% ifelse(
    1,0); ensanut_2022$highBP[is.na(ensanut_2022$a0404) & is.na(ensanut_2022$PAS_fin) & is.na(ensanut_2022$PAD_fin)] <- NA
#MetS
ensanut_2022$MetS <- ((ensanut_2022 %>% select(highTG, lowHDL, highBP) %>% apply(1, sum, na.rm=T))>0) %>% ifelse(1,0) #Metabolic Syndrome (IDF)
ensanut_2022$MetS[with(ensanut_2022, (is.na(highTG) & is.na(lowHDL) & is.na(highBP)))] <- NA; ensanut_2022$MetS[ensanut_2022$obesity==0]<-0
ensanut_2022$MetS %>% table(useNA = "always")
ensanut_2022$MetS_n <- ensanut_2022 %>% select(obesity, highTG, lowHDL, highBP) %>% apply(1, sum, na.rm=T) #Number of criteria fulfilled
ensanut_2022$MetS_n[with(ensanut_2022, (is.na(highTG) & is.na(lowHDL) & is.na(highBP) & is.na(obesity)))] <- NA


##-- Chronic kidney disease --##
#2016
ensanut_2016$CKD_fin <- ifelse(ensanut_2016$a605c==1,1,0)
ensanut_2016$CKD_fin %>% table(useNA = "always")
#2018
ensanut_2018$CKD_fin <- ifelse(ensanut_2018$P6_1_3.x==1,1,0)
ensanut_2018$CKD_fin %>% table(useNA = "always")
#2021
ensanut_2021$CKD_fin <- ifelse(ensanut_2021$a0601c==1,1,0)
ensanut_2021$CKD_fin %>% table(useNA = "always")
#2022
ensanut_2022$CKD_fin <- ifelse(ensanut_2022$a0601c==1,1,0)
ensanut_2022$CKD_fin %>% table(useNA = "always")

## --- # --- # CARDIOVASCULAR DISEASES # --- # --- ##
##-- Myocardial infarction --##
#2016
ensanut_2016$a502a[ensanut_2016$a502a==3]<-NA
ensanut_2016$AMI_fin <- ifelse(ensanut_2016$a502a==1,1,0)
ensanut_2016$AMI_fin %>% table(useNA = "always")
#2018
ensanut_2018$AMI_fin <- ifelse(ensanut_2018$P5_2_1==1,1,0)
ensanut_2018$AMI_fin %>% table(useNA = "always")
#2021
ensanut_2021$a0502a <- replace_na(ensanut_2021$a0502a, 2)
ensanut_2021$AMI_fin <- ifelse(ensanut_2021$a0502a==1,1,0)
ensanut_2021$AMI_fin %>% table(useNA = "always")
#2022
ensanut_2022$a0502a <- replace_na(ensanut_2022$a0502a, 2)
ensanut_2022$AMI_fin <- ifelse(ensanut_2022$a0502a==1,1,0)
ensanut_2022$AMI_fin %>% table(useNA = "always")

##-- Angina --##
#2016
ensanut_2016$ANG_fin <- ifelse(ensanut_2016$a605b==1,1,0)
ensanut_2016$ANG_fin %>% table(useNA = "always")
#2018
ensanut_2018$ANG_fin <- ifelse(ensanut_2018$P5_2_2==1,1,0)
ensanut_2018$ANG_fin %>% table(useNA = "always")
#2021
ensanut_2021$a0502b <- replace_na(ensanut_2021$a0502b, 2)
ensanut_2021$ANG_fin <- ifelse(ensanut_2021$a0502b==1,1,0)
ensanut_2021$ANG_fin %>% table(useNA = "always")
#2022
ensanut_2022$a0502b <- replace_na(ensanut_2022$a0502b, 2)
ensanut_2022$ANG_fin <- ifelse(ensanut_2022$a0502b==1,1,0)
ensanut_2022$ANG_fin %>% table(useNA = "always")

##-- Heart failure --##
#2016
ensanut_2016$a502c[ensanut_2016$a502c==3]<-NA
ensanut_2016$HFA_fin <- ifelse(ensanut_2016$a502c==1,1,0)
ensanut_2016$HFA_fin %>% table(useNA = "always")
#2018
ensanut_2018$HFA_fin <- ifelse(ensanut_2018$P5_2_3==1,1,0)
ensanut_2018$HFA_fin %>% table(useNA = "always")
#2021
ensanut_2021$a0502c <- replace_na(ensanut_2021$a0502c, 2)
ensanut_2021$HFA_fin <- ifelse(ensanut_2021$a0502c==1,1,0)
ensanut_2021$HFA_fin %>% table(useNA = "always")
#2022
ensanut_2022$a0502c <- replace_na(ensanut_2022$a0502c, 2)
ensanut_2022$HFA_fin <- ifelse(ensanut_2022$a0502c==1,1,0)
ensanut_2022$HFA_fin %>% table(useNA = "always")

##-- Stroke --##
#2016
ensanut_2016$a611[ensanut_2016$a611==3]<-NA
ensanut_2016$EVC_fin <- ifelse(ensanut_2016$a611==1,1,0)
ensanut_2016$EVC_fin %>% table(useNA = "always")
#2018
ensanut_2018$P5_6 <- replace_na(ensanut_2018$P5_6, 2)
ensanut_2018$P5_6[ensanut_2018$P5_6==9]<-NA
ensanut_2018$EVC_fin <- ifelse(ensanut_2018$P5_6==1,1,0)
ensanut_2018$EVC_fin %>% table(useNA = "always")
#2021
ensanut_2021$a0506[ensanut_2021$a0506==9]<-NA
ensanut_2021$EVC_fin <- ifelse(ensanut_2021$a0506==1,1,0)
ensanut_2021$EVC_fin %>% table(useNA = "always")
#2022
ensanut_2022$a0502d[ensanut_2022$a0502d==9]<-NA
ensanut_2022$EVC_fin <- ifelse(ensanut_2022$a0502d==1,1,0)
ensanut_2022$EVC_fin %>% table(useNA = "always")

##-- Cardiovascular disease (MI | Heart failure | Stroke) --##
#2016
ensanut_2016$CVD_fin <- ((ensanut_2016 %>% select(
  AMI_fin, HFA_fin, EVC_fin) %>% apply(1, sum, na.rm=T))==0) %>% ifelse(0,1)
ensanut_2016$CVD_fin[with(ensanut_2016, is.na(a502a)&is.na(a502c)&is.na(a611))] <- NA
ensanut_2016$CVD_fin %>% table(useNA = "always") %>% prop.table()
#2018
ensanut_2018$CVD_fin <- ((ensanut_2018 %>% select(
  AMI_fin, HFA_fin, EVC_fin) %>% apply(1, sum, na.rm=T))==0) %>% ifelse(0,1)
ensanut_2018$CVD_fin[with(ensanut_2018, is.na(P5_2_1)&is.na(P5_2_3)&is.na(P5_6))] <- NA
ensanut_2018$CVD_fin %>% table(useNA = "always") %>% prop.table()
#2020
ensanut_2020$CVD_fin <- with(ensanut_2020,ifelse((H0902A.x=="4" | H0902B=="4" | H0902C=="4" | H0902D=="4")==T,1,0))
ensanut_2020$CVD_fin[with(ensanut_2020, is.na(H0902A.x)&is.na(H0902B)&is.na(H0902C)&is.na(H0902D))] <- NA
ensanut_2020$CVD_fin %>% table(useNA = "always") %>% prop.table()
#2021
ensanut_2021$CVD_fin <- ((ensanut_2021 %>% select(
  AMI_fin, HFA_fin, EVC_fin) %>% apply(1, sum, na.rm=T))==0) %>% ifelse(0,1)
ensanut_2021$CVD_fin[with(ensanut_2021, is.na(a0502a)&is.na(a0502c)&is.na(a0506))] <- NA
ensanut_2021$CVD_fin %>% table(useNA = "always") %>% prop.table()
#2022
ensanut_2022$CVD_fin <- ((ensanut_2022 %>% select(
  AMI_fin, HFA_fin, EVC_fin) %>% apply(1, sum, na.rm=T))==0) %>% ifelse(0,1)
ensanut_2022$CVD_fin[with(ensanut_2022, is.na(a0502a)&is.na(a0502c)&is.na(a0502d))] <- NA
ensanut_2022$CVD_fin %>% table(useNA = "always") %>% prop.table()


##-- LDL elevation --##
#2016
ensanut_2016$LDL_fin <- ensanut_2016$valor.COL_LDL; ensanut_2016$LDL_fin[ensanut_2016$sanvenh<8]<-NA
ensanut_2016$High.LDL <- (with(ensanut_2016, ((LDL_fin>=130)))) %>% ifelse(1,0)
ensanut_2016$High.LDL %>% table(useNA = "always"); ensanut_2016$High.LDL %>% table() %>% prop.table
#2018
ensanut_2018$LDL_fin <- ensanut_2018$VALOR_COL_LDL; ensanut_2018$LDL_fin[ensanut_2018$P5_1.y<8]<-NA
ensanut_2018$High.LDL <- (with(ensanut_2018, ((LDL_fin>=130)))) %>% ifelse(1,0)
ensanut_2018$High.LDL %>% table(useNA = "always"); ensanut_2018$High.LDL %>% table() %>% prop.table
#2020
ensanut_2020$LDL_fin <- ensanut_2020$valor.COL_LDL; ensanut_2020$LDL_fin[ensanut_2020$san04<8]<-NA
ensanut_2020$High.LDL <- (with(ensanut_2020, ((LDL_fin>=130)))) %>% ifelse(1,0)
ensanut_2020$High.LDL %>% table(useNA = "always"); ensanut_2020$High.LDL %>% table() %>% prop.table
#2021
ensanut_2021$LDL_fin <- ensanut_2021$valor_COL_LDL; ensanut_2021$LDL_fin[ensanut_2021$san04<8]<-NA
ensanut_2021$High.LDL <- (with(ensanut_2021, ((LDL_fin>=130)))) %>% ifelse(1,0)
ensanut_2021$High.LDL %>% table(useNA = "always"); ensanut_2021$High.LDL %>% table() %>% prop.table
#2022
ensanut_2022$LDL_fin <- ensanut_2022$valor_COL_LDL; ensanut_2022$LDL_fin[ensanut_2022$san04<8]<-NA
ensanut_2022$High.LDL <- (with(ensanut_2022, ((LDL_fin>=130)))) %>% ifelse(1,0)
ensanut_2022$High.LDL %>% table(useNA = "always"); ensanut_2022$High.LDL %>% table() %>% prop.table

##-- CRP elevation --##
#Joint cutoff
ensanut_2021$PCR_fin <- ensanut_2021$valor_PROTCREAC; ensanut_2021$PCR_fin[ensanut_2021$san04<8]<-NA
ensanut_2022$PCR_fin <- ensanut_2022$valor_PCR; ensanut_2022$PCR_fin[ensanut_2022$san04<8]<-NA
CRP.cutoff <- rbind(ensanut_2021 %>% select(PCR_fin), ensanut_2022 %>% select(PCR_fin)) %>%
  select(PCR_fin) %>% apply(2, quantile, 0.8, na.rm=T)
#2021
ensanut_2021$High.PCR <- (with(ensanut_2021, ((PCR_fin>= CRP.cutoff) ))) %>% ifelse(1,0)
ensanut_2021$High.PCR %>% table(useNA = "always"); ensanut_2021$High.PCR %>% table() %>% prop.table
#2022
ensanut_2022$High.PCR <- (with(ensanut_2022, ((PCR_fin>= CRP.cutoff)))) %>% ifelse(1,0)
ensanut_2022$High.PCR %>% table(useNA = "always"); ensanut_2022$High.PCR %>% table() %>% prop.table

##-- GFR reduction --##
#2016
ensanut_2016$Cr_fin <- ensanut_2016$valor.CREAT; ensanut_2016$Cr_fin[ensanut_2016$sanvenh<8]<-NA
ensanut_2016$Sex_GFR <- ifelse(ensanut_2016$sexo_2016==1,0,1); ensanut_2016$Eth_GFR <- 0
ensanut_2016$GFR <- with(ensanut_2016, nephro::CKDEpi.creat(Cr_fin, Sex_GFR, edad.x, Eth_GFR))
ensanut_2016$CKD_lab <- (ensanut_2016$GFR<60) %>% ifelse(1,0); ensanut_2016$CKD_lab%>%table()%>%prop.table
#2018
ensanut_2018$Cr_fin <- ensanut_2018$VALOR_CREAT; ensanut_2018$Cr_fin[ensanut_2018$P5_1.y<8]<-NA
ensanut_2018$Sex_GFR <- ifelse(ensanut_2018$sexo_2018==1,0,1); ensanut_2018$Eth_GFR <- 0
ensanut_2018$GFR <- with(ensanut_2018, nephro::CKDEpi.creat(Cr_fin, Sex_GFR, EDAD.x, Eth_GFR))
ensanut_2018$CKD_lab <- (ensanut_2018$GFR<60) %>% ifelse(1,0); ensanut_2018$CKD_lab%>%table()%>%prop.table
#2020
ensanut_2020$Cr_fin <- ensanut_2020$valor.CREAT; ensanut_2020$Cr_fin[ensanut_2020$san04<8]<-NA
ensanut_2020$Sex_GFR <- ifelse(ensanut_2020$sexo_2020==1,0,1); ensanut_2020$Eth_GFR <- 0
ensanut_2020$GFR <- with(ensanut_2020, nephro::CKDEpi.creat(Cr_fin, Sex_GFR, edad_num, Eth_GFR))
ensanut_2020$CKD_lab <- (ensanut_2020$GFR<60) %>% ifelse(1,0); ensanut_2020$CKD_lab%>%table()%>%prop.table
#2021
ensanut_2021$Cr_fin <- ensanut_2021$valor_CREAT; ensanut_2021$Cr_fin[ensanut_2021$san04<8]<-NA
ensanut_2021$Sex_GFR <- ifelse(ensanut_2021$sexo_2021==1,0,1); ensanut_2021$Eth_GFR <- 0
ensanut_2021$GFR <- with(ensanut_2021, nephro::CKDEpi.creat(Cr_fin, Sex_GFR, edad_num, Eth_GFR))
ensanut_2021$CKD_lab <- (ensanut_2021$GFR<60) %>% ifelse(1,0); ensanut_2021$CKD_lab%>%table()%>%prop.table
#2022
ensanut_2022$Cr_fin <- ensanut_2022$valor_CREAT; ensanut_2022$Cr_fin[ensanut_2022$san04<8]<-NA
ensanut_2022$Sex_GFR <- ifelse(ensanut_2022$sexo_2022==1,0,1); ensanut_2022$Eth_GFR <- 0
ensanut_2022$GFR <- with(ensanut_2022, nephro::CKDEpi.creat(Cr_fin, Sex_GFR, edad_num, Eth_GFR))
ensanut_2022$CKD_lab <- (ensanut_2022$GFR<60) %>% ifelse(1,0); ensanut_2022$CKD_lab%>%table()%>%prop.table



#### Prediabetes definitions  ####
# pred_comp1: FPG ≥100 (ADA)
# pred_comp2: HbA1c ≥5.7 (ADA)
# pred_comp3: Any ADA criteria
# pred_comp4: Both ADA criteria
# pred_comp5: FPG ≥110 (WHO)
# pred_comp6: HbA1c ≥6.0 (IEC)

# 2016
ensanut_2016$pred_comp1 <- ensanut_2016$prediabetes_ifg_2016; ensanut_2016$pred_comp1[ensanut_2016$diabetes_fin==1] <- 2
ensanut_2016$pred_comp2 <- ensanut_2016$prediabetes_hb1ac_2016; ensanut_2016$pred_comp2[ensanut_2016$diabetes_fin==1] <- 2
ensanut_2016$pred_comp3 <- ensanut_2016$prediabetes_allcriteria_2016; ensanut_2016$pred_comp3[ensanut_2016$diabetes_fin==1] <- 2
ensanut_2016$pred_comp4 <- ensanut_2016$prediabetes.prev_ifg_hb1ac_2016; ensanut_2016$pred_comp4[ensanut_2016$diabetes_fin==1] <- 2
ensanut_2016$pred_comp5 <- ensanut_2016$ifg_who; ensanut_2016$pred_comp5[ensanut_2016$diabetes_fin==1] <- 2
ensanut_2016$pred_comp6 <- ensanut_2016$a1c_iec; ensanut_2016$pred_comp6[ensanut_2016$diabetes_fin==1] <- 2
# 2018
ensanut_2018$pred_comp1 <- ensanut_2018$prediabetes_ifg_2018; ensanut_2018$pred_comp1[ensanut_2018$diabetes_fin==1] <- 2
ensanut_2018$pred_comp2 <- ensanut_2018$prediabetes_hb1ac_2018; ensanut_2018$pred_comp2[ensanut_2018$diabetes_fin==1] <- 2
ensanut_2018$pred_comp3 <- ensanut_2018$prediabetes_allcriteria_2018; ensanut_2018$pred_comp3[ensanut_2018$diabetes_fin==1] <- 2
ensanut_2018$pred_comp4 <- ensanut_2018$prediabetes.prev_ifg_hb1ac_2018; ensanut_2018$pred_comp4[ensanut_2018$diabetes_fin==1] <- 2
ensanut_2018$pred_comp5 <- ensanut_2018$ifg_who; ensanut_2018$pred_comp5[ensanut_2018$diabetes_fin==1] <- 2
ensanut_2018$pred_comp6 <- ensanut_2018$a1c_iec; ensanut_2018$pred_comp6[ensanut_2018$diabetes_fin==1] <- 2
# 2020
ensanut_2020$pred_comp1 <- ensanut_2020$prediabetes_ifg_2020; ensanut_2020$pred_comp1[ensanut_2020$diabetes_fin==1] <- 2
ensanut_2020$pred_comp2 <- ensanut_2020$prediabetes_hb1ac_2020; ensanut_2020$pred_comp2[ensanut_2020$diabetes_fin==1] <- 2
ensanut_2020$pred_comp3 <- ensanut_2020$prediabetes_allcriteria_2020; ensanut_2020$pred_comp3[ensanut_2020$diabetes_fin==1] <- 2
ensanut_2020$pred_comp4 <- ensanut_2020$prediabetes.prev_ifg_hb1ac_2020; ensanut_2020$pred_comp4[ensanut_2020$diabetes_fin==1] <- 2
ensanut_2020$pred_comp5 <- ensanut_2020$ifg_who; ensanut_2020$pred_comp5[ensanut_2020$diabetes_fin==1] <- 2
ensanut_2020$pred_comp6 <- ensanut_2020$a1c_iec; ensanut_2020$pred_comp6[ensanut_2020$diabetes_fin==1] <- 2
# 2021
ensanut_2021$pred_comp1 <- ensanut_2021$prediabetes_ifg_2021; ensanut_2021$pred_comp1[ensanut_2021$diabetes_fin==1] <- 2
ensanut_2021$pred_comp2 <- ensanut_2021$prediabetes_hb1ac_2021; ensanut_2021$pred_comp2[ensanut_2021$diabetes_fin==1] <- 2
ensanut_2021$pred_comp3 <- ensanut_2021$prediabetes_allcriteria_2021; ensanut_2021$pred_comp3[ensanut_2021$diabetes_fin==1] <- 2
ensanut_2021$pred_comp4 <- ensanut_2021$prediabetes.prev_ifg_hb1ac_2021; ensanut_2021$pred_comp4[ensanut_2021$diabetes_fin==1] <- 2
ensanut_2021$pred_comp5 <- ensanut_2021$ifg_who; ensanut_2021$pred_comp5[ensanut_2021$diabetes_fin==1] <- 2
ensanut_2021$pred_comp6 <- ensanut_2021$a1c_iec; ensanut_2021$pred_comp6[ensanut_2021$diabetes_fin==1] <- 2
# 2022
ensanut_2022$pred_comp1 <- ensanut_2022$prediabetes_ifg_2022; ensanut_2022$pred_comp1[ensanut_2022$diabetes_fin==1] <- 2
ensanut_2022$pred_comp2 <- ensanut_2022$prediabetes_hb1ac_2022; ensanut_2022$pred_comp2[ensanut_2022$diabetes_fin==1] <- 2
ensanut_2022$pred_comp3 <- ensanut_2022$prediabetes_allcriteria_2022; ensanut_2022$pred_comp3[ensanut_2022$diabetes_fin==1] <- 2
ensanut_2022$pred_comp4 <- ensanut_2022$prediabetes.prev_ifg_hb1ac_2022; ensanut_2022$pred_comp4[ensanut_2022$diabetes_fin==1] <- 2
ensanut_2022$pred_comp5 <- ensanut_2022$ifg_who; ensanut_2022$pred_comp5[ensanut_2022$diabetes_fin==1] <- 2
ensanut_2022$pred_comp6 <- ensanut_2022$a1c_iec; ensanut_2022$pred_comp6[ensanut_2022$diabetes_fin==1] <- 2


#### Pooled standard weights  ####
#ENSANUT 2016
ens16_fin.2 <- ensanut_2016 %>% mutate(
  "pond.lab"=ponde_f_vv, "estr.lab"=est_var, "sampleid"=paste(entidad_2016, as.numeric(estr.lab))) %>%
  filter(!is.na(pond.lab),!is.na(estr.lab)); wts <- with(ens16_fin.2, tapply(pond.lab, sampleid, sum)) %>%
  data.frame() %>% `names<-`("wt") %>% rownames_to_column("ids")
t1<-as.data.frame(table(ens16_fin.2$sampleid)); wts2<-data.frame(ids=wts$ids, sumwt=wts$wt, jn=t1$Freq)
ens16_fin.2 <- merge(ens16_fin.2, wts2, by.x="sampleid", by.y="ids", all.x=T)
ens16_fin.2$swts <- ens16_fin.2$pond.lab*(ens16_fin.2$jn/ens16_fin.2$sumwt)

#ENSANUT 2018
ens18_fin.2 <- ensanut_2018 %>% mutate(
  "pond.lab"=ponderador_glucosa, "estr.lab"=ESTRATO.y, "sampleid"=paste(entidad_2018, as.numeric(estr.lab))) %>%
  filter(!is.na(pond.lab),!is.na(estr.lab)); wts <- with(ens18_fin.2, tapply(pond.lab, sampleid, sum)) %>%
  data.frame() %>% `names<-`("wt") %>% rownames_to_column("ids")
t1<-as.data.frame(table(ens18_fin.2$sampleid)); wts2<-data.frame(ids=wts$ids, sumwt=wts$wt, jn=t1$Freq)
ens18_fin.2 <- merge(ens18_fin.2, wts2, by.x="sampleid", by.y="ids", all.x=T)
ens18_fin.2$swts <- ens18_fin.2$pond.lab*(ens18_fin.2$jn/ens18_fin.2$sumwt)

#ENSANUT 2020
ens20_fin.2 <- ensanut_2020 %>% mutate(
  "pond.lab"=ponde_g20.y, "estr.lab"=est_sel.x, "sampleid"=paste(entidad_2020, as.numeric(estr.lab))) %>%
  filter(!estr.lab%in%(which((estr.lab%>%table)<2)%>%names%>%as.numeric)) %>% 
  filter(!is.na(pond.lab),!is.na(estr.lab)); wts <- with(ens20_fin.2, tapply(pond.lab, sampleid, sum)) %>%
  data.frame() %>% `names<-`("wt") %>% rownames_to_column("ids")
t1<-as.data.frame(table(ens20_fin.2$sampleid)); wts2<-data.frame(ids=wts$ids, sumwt=wts$wt, jn=t1$Freq)
ens20_fin.2 <- merge(ens20_fin.2, wts2, by.x="sampleid", by.y="ids", all.x=T)
ens20_fin.2$swts <- ens20_fin.2$pond.lab*(ens20_fin.2$jn/ens20_fin.2$sumwt)

#ENSANUT 2021
ens21_fin.2 <- ensanut_2021 %>% mutate(
  "pond.lab"=ponde_g, "estr.lab"=est_sel.lab, "sampleid"=paste(entidad_2021, as.numeric(estr.lab))) %>%
  filter(!estr.lab%in%(which((estr.lab%>%table)<2)%>%names%>%as.numeric)) %>%
  filter(!is.na(pond.lab),!is.na(estr.lab)); wts <- with(ens21_fin.2, tapply(pond.lab, sampleid, sum)) %>%
  data.frame() %>% `names<-`("wt") %>% rownames_to_column("ids")
t1<-as.data.frame(table(ens21_fin.2$sampleid)); wts2<-data.frame(ids=wts$ids, sumwt=wts$wt, jn=t1$Freq)
ens21_fin.2 <- merge(ens21_fin.2, wts2, by.x="sampleid", by.y="ids", all.x=T)
ens21_fin.2$swts <- ens21_fin.2$pond.lab*(ens21_fin.2$jn/ens21_fin.2$sumwt)

#Ensanut 2022
ens22_fin.2 <- ensanut_2022 %>% mutate(
  "pond.lab"=ponde_v, "estr.lab"=est_sel.x.x.x, "sampleid"=paste(entidad_2022, as.numeric(estr.lab))) %>%
  filter(!is.na(pond.lab),!is.na(estr.lab)); wts <- with(ens22_fin.2, tapply(pond.lab, sampleid, sum)) %>%
  data.frame() %>% `names<-`("wt") %>% rownames_to_column("ids")
t1<-as.data.frame(table(ens22_fin.2$sampleid)); wts2<-data.frame(ids=wts$ids, sumwt=wts$wt, jn=t1$Freq)
ens22_fin.2 <- merge(ens22_fin.2, wts2, by.x="sampleid", by.y="ids", all.x=T)
ens22_fin.2$swts <- ens22_fin.2$pond.lab*(ens22_fin.2$jn/ens22_fin.2$sumwt)

list(ens16_fin.2,ens18_fin.2,ens20_fin.2,ens21_fin.2,ens22_fin.2) %>% lapply(nrow)

#### Logistic - Arterial Hypertension ####
# pred_comp1: FPG ≥100 (ADA)
rbind(ens16_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp1, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp1, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp1, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp1, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp1, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTA.PRED1

# pred_comp2: HbA1c ≥5.7 (ADA)
rbind(ens16_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp2, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp2, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp2, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp2, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp2, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTA.PRED2

# pred_comp3: Any ADA criteria
rbind(ens16_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp3, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp3, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp3, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp3, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp3, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTA.PRED3

# pred_comp4: Both ADA criteria
rbind(ens16_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp4, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp4, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp4, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp4, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp4, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTA.PRED4

# pred_comp5: FPG ≥110 (WHO)
rbind(ens16_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp5, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp5, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp5, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp5, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp5, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTA.PRED5

# pred_comp6: HbA1c ≥6.0 (IEC)
rbind(ens16_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp6, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp6, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp6, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp6, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HTA_fin, "pdm"=pred_comp6, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTA.PRED6

R.HTA <- paste0("R.HTA.PRED",1:6) %>% lapply(get) %>% sapply(R.getOR2) %>% t %>% data.frame("Labs"=1:6) %>% 
  mutate("Outcome"=1, "Labs"=ordered(
    Labs, levels=6:1, labels=c("HbA1c (IEC)", "IFG (WHO)", "Both (ADA)", "Any (ADA)", "HbA1c (ADA)", "IFG (ADA)")))


#### Logistic - Hyper Cholesterolemia ####
# pred_comp1: FPG ≥100 (ADA) 
rbind(ens16_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp1, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp1, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp1, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp1, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp1, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTC.PRED1
# pred_comp2: HbA1c ≥5.7 (ADA)
rbind(ens16_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp2, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp2, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp2, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp2, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp2, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTC.PRED2

# pred_comp3: Any ADA criteria
rbind(ens16_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp3, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp3, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp3, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp3, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp3, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTC.PRED3

# pred_comp4: Both ADA criteria
rbind(ens16_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp4, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp4, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp4, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp4, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp4, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTC.PRED4

# pred_comp5: FPG ≥110 (WHO)
rbind(ens16_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp5, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp5, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp5, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp5, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp5, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTC.PRED5

# pred_comp6: HbA1c ≥6.0 (IEC)
rbind(ens16_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp6, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp6, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp6, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp6, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HCT_fin, "pdm"=pred_comp6, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTC.PRED6

R.HTC <- paste0("R.HTC.PRED",1:6) %>% lapply(get) %>% sapply(R.getOR2) %>% t %>% data.frame("Labs"=1:6) %>% 
  mutate("Outcome"=2, "Labs"=ordered(
    Labs, levels=6:1, labels=c("HbA1c (IEC)", "IFG (WHO)", "Both (ADA)", "Any (ADA)", "HbA1c (ADA)", "IFG (ADA)")))


#### Logistic - Hyper Triglyceridemia ####
# pred_comp1: FPG ≥100 (ADA) 
rbind(ens16_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp1, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp1, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp1, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp1, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp1, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTG.PRED1

# pred_comp2: HbA1c ≥5.7 (ADA)
rbind(ens16_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp2, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp2, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp2, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp2, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp2, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTG.PRED2

# pred_comp3: Any ADA criteria
rbind(ens16_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp3, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp3, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp3, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp3, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp3, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTG.PRED3

# pred_comp4: Both ADA criteria
rbind(ens16_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp4, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp4, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp4, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp4, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp4, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTG.PRED4

# pred_comp5: FPG ≥110 (WHO)
rbind(ens16_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp5, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp5, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp5, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp5, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp5, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTG.PRED5

# pred_comp6: HbA1c ≥6.0 (IEC)
rbind(ens16_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp6, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp6, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp6, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp6, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=HTG_fin, "pdm"=pred_comp6, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.HTG.PRED6

R.HTG <- paste0("R.HTG.PRED",1:6) %>% lapply(get) %>% sapply(R.getOR2) %>% t %>% data.frame("Labs"=1:6) %>% 
  mutate("Outcome"=3, "Labs"=ordered(
    Labs, levels=6:1, labels=c("HbA1c (IEC)", "IFG (WHO)", "Both (ADA)", "Any (ADA)", "HbA1c (ADA)", "IFG (ADA)")))


#### Logistic - Insulin Resistance ####
# pred_comp1: FPG ≥100 (ADA) 
rbind(ens16_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp1, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp1, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp1, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp1, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp1, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
   glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.IR.PRED1

# pred_comp2: HbA1c ≥5.7 (ADA)
rbind(ens16_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp2, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp2, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp2, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp2, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp2, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
   glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.IR.PRED2

# pred_comp3: Any ADA criteria
rbind(ens16_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp3, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp3, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp3, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp3, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp3, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
   glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.IR.PRED3

# pred_comp4: Both ADA criteria
rbind(ens16_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp4, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp4, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp4, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp4, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp4, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
   glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.IR.PRED4

# pred_comp5: FPG ≥110 (WHO)
rbind(ens16_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp5, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp5, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp5, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp5, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp5, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
   glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.IR.PRED5

# pred_comp6: HbA1c ≥6.0 (IEC)
rbind(ens16_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp6, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp6, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp6, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp6, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=homa_cat2, "pdm"=pred_comp6, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
   glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.IR.PRED6

R.IR <- paste0("R.IR.PRED",1:6) %>% lapply(get) %>% sapply(R.getOR2) %>% t %>% data.frame("Labs"=1:6) %>% 
  mutate("Outcome"=4, "Labs"=ordered(Labs, levels=6:1, labels=c(
    "HbA1c (IEC)", "IFG (WHO)", "Both (ADA)", "Any (ADA)", "HbA1c (ADA)", "IFG (ADA)")))


#### Logistic - Metabolic Syndrome ####
# pred_comp1: FPG ≥100 (ADA) 
rbind(ens16_fin.2 %>% select("out"=MetS, "pdm"=pred_comp1, "sex"=sexo_2016, YEAR, swts, "age"=edad.x),
      ens18_fin.2 %>% select("out"=MetS, "pdm"=pred_comp1, "sex"=sexo_2018, YEAR, swts, "age"=EDAD.x),
      ens20_fin.2 %>% select("out"=MetS, "pdm"=pred_comp1, "sex"=sexo_2020, YEAR, swts, "age"=edad_num),
      ens21_fin.2 %>% select("out"=MetS, "pdm"=pred_comp1, "sex"=sexo_2021, YEAR, swts, "age"=edad_num),
      ens22_fin.2 %>% select("out"=MetS, "pdm"=pred_comp1, "sex"=sexo_2022, YEAR, swts, "age"=edad_num)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+(1|YEAR), family = binomial(), weights = swts) -> R.METS.PRED1

# pred_comp2: HbA1c ≥5.7 (ADA)
rbind(ens16_fin.2 %>% select("out"=MetS, "pdm"=pred_comp2, "sex"=sexo_2016, YEAR, swts, "age"=edad.x),
      ens18_fin.2 %>% select("out"=MetS, "pdm"=pred_comp2, "sex"=sexo_2018, YEAR, swts, "age"=EDAD.x),
      ens20_fin.2 %>% select("out"=MetS, "pdm"=pred_comp2, "sex"=sexo_2020, YEAR, swts, "age"=edad_num),
      ens21_fin.2 %>% select("out"=MetS, "pdm"=pred_comp2, "sex"=sexo_2021, YEAR, swts, "age"=edad_num),
      ens22_fin.2 %>% select("out"=MetS, "pdm"=pred_comp2, "sex"=sexo_2022, YEAR, swts, "age"=edad_num)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+(1|YEAR), family = binomial(), weights = swts) -> R.METS.PRED2

# pred_comp3: Any ADA criteria
rbind(ens16_fin.2 %>% select("out"=MetS, "pdm"=pred_comp3, "sex"=sexo_2016, YEAR, swts, "age"=edad.x),
      ens18_fin.2 %>% select("out"=MetS, "pdm"=pred_comp3, "sex"=sexo_2018, YEAR, swts, "age"=EDAD.x),
      ens20_fin.2 %>% select("out"=MetS, "pdm"=pred_comp3, "sex"=sexo_2020, YEAR, swts, "age"=edad_num),
      ens21_fin.2 %>% select("out"=MetS, "pdm"=pred_comp3, "sex"=sexo_2021, YEAR, swts, "age"=edad_num),
      ens22_fin.2 %>% select("out"=MetS, "pdm"=pred_comp3, "sex"=sexo_2022, YEAR, swts, "age"=edad_num)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+(1|YEAR), family = binomial(), weights = swts) -> R.METS.PRED3

# pred_comp4: Both ADA criteria
rbind(ens16_fin.2 %>% select("out"=MetS, "pdm"=pred_comp4, "sex"=sexo_2016, YEAR, swts, "age"=edad.x),
      ens18_fin.2 %>% select("out"=MetS, "pdm"=pred_comp4, "sex"=sexo_2018, YEAR, swts, "age"=EDAD.x),
      ens20_fin.2 %>% select("out"=MetS, "pdm"=pred_comp4, "sex"=sexo_2020, YEAR, swts, "age"=edad_num),
      ens21_fin.2 %>% select("out"=MetS, "pdm"=pred_comp4, "sex"=sexo_2021, YEAR, swts, "age"=edad_num),
      ens22_fin.2 %>% select("out"=MetS, "pdm"=pred_comp4, "sex"=sexo_2022, YEAR, swts, "age"=edad_num)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+(1|YEAR), family = binomial(), weights = swts) -> R.METS.PRED4

# pred_comp5: FPG ≥110 (WHO)
rbind(ens16_fin.2 %>% select("out"=MetS, "pdm"=pred_comp5, "sex"=sexo_2016, YEAR, swts, "age"=edad.x),
      ens18_fin.2 %>% select("out"=MetS, "pdm"=pred_comp5, "sex"=sexo_2018, YEAR, swts, "age"=EDAD.x),
      ens20_fin.2 %>% select("out"=MetS, "pdm"=pred_comp5, "sex"=sexo_2020, YEAR, swts, "age"=edad_num),
      ens21_fin.2 %>% select("out"=MetS, "pdm"=pred_comp5, "sex"=sexo_2021, YEAR, swts, "age"=edad_num),
      ens22_fin.2 %>% select("out"=MetS, "pdm"=pred_comp5, "sex"=sexo_2022, YEAR, swts, "age"=edad_num)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+(1|YEAR), family = binomial(), weights = swts) -> R.METS.PRED5

# pred_comp6: HbA1c ≥6.0 (IEC)
rbind(ens16_fin.2 %>% select("out"=MetS, "pdm"=pred_comp6, "sex"=sexo_2016, YEAR, swts, "age"=edad.x),
      ens18_fin.2 %>% select("out"=MetS, "pdm"=pred_comp6, "sex"=sexo_2018, YEAR, swts, "age"=EDAD.x),
      ens20_fin.2 %>% select("out"=MetS, "pdm"=pred_comp6, "sex"=sexo_2020, YEAR, swts, "age"=edad_num),
      ens21_fin.2 %>% select("out"=MetS, "pdm"=pred_comp6, "sex"=sexo_2021, YEAR, swts, "age"=edad_num),
      ens22_fin.2 %>% select("out"=MetS, "pdm"=pred_comp6, "sex"=sexo_2022, YEAR, swts, "age"=edad_num)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+(1|YEAR), family = binomial(), weights = swts) -> R.METS.PRED6

R.METS <- paste0("R.METS.PRED",1:6) %>% lapply(get) %>% sapply(R.getOR2) %>% t %>% data.frame("Labs"=1:6) %>% 
  mutate("Outcome"=5, "Labs"=ordered(
    Labs, levels=6:1, labels=c("HbA1c (IEC)", "IFG (WHO)", "Both (ADA)", "Any (ADA)", "HbA1c (ADA)", "IFG (ADA)")))


#### Logistic - CVD: AMI|HF|Stroke ####
# pred_comp1: FPG ≥100 (ADA) 
rbind(ens16_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp1, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp1, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp1, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020), 
      ens21_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp1, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp1, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial()) -> R.CVD.PRED1

# pred_comp2: HbA1c ≥5.7 (ADA)
rbind(ens16_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp2, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp2, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp2, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp2, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp2, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.CVD.PRED2

# pred_comp3: Any ADA criteria
rbind(ens16_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp3, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp3, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp3, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp3, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp3, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.CVD.PRED3

# pred_comp4: Both ADA criteria
rbind(ens16_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp4, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp4, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp4, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp4, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp4, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.CVD.PRED4

# pred_comp5: FPG ≥110 (WHO)
rbind(ens16_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp5, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp5, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp5, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp5, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp5, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.CVD.PRED5

# pred_comp6: HbA1c ≥6.0 (IEC)
rbind(ens16_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp6, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp6, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp6, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp6, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=CVD_fin, "pdm"=pred_comp6, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.CVD.PRED6

R.CVD <- paste0("R.CVD.PRED",1:6) %>% lapply(get) %>% sapply(R.getOR2) %>% t %>% data.frame("Labs"=1:6) %>% 
  mutate("Outcome"=6, "Labs"=ordered(
    Labs, levels=6:1, labels=c("HbA1c (IEC)", "IFG (WHO)", "Both (ADA)", "Any (ADA)", "HbA1c (ADA)", "IFG (ADA)")))


#### Logistic - High LDL ####
# pred_comp1: FPG ≥100 (ADA) 
rbind(ens16_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp1, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp1, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp1, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020), 
      ens21_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp1, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp1, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.LDL.PRED1

# pred_comp2: HbA1c ≥5.7 (ADA)
rbind(ens16_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp2, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp2, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp2, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp2, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp2, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.LDL.PRED2

# pred_comp3: Any ADA criteria
rbind(ens16_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp3, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp3, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp3, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp3, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp3, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.LDL.PRED3

# pred_comp4: Both ADA criteria
rbind(ens16_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp4, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp4, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp4, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp4, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp4, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.LDL.PRED4

# pred_comp5: FPG ≥110 (WHO)
rbind(ens16_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp5, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp5, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp5, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp5, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp5, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.LDL.PRED5

# pred_comp6: HbA1c ≥6.0 (IEC)
rbind(ens16_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp6, "sex"=sexo_2016, "age"=edad.x, YEAR, swts, "bmi"=imc),
      ens18_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp6, "sex"=sexo_2018, "age"=EDAD.x, YEAR, swts, "bmi"=imc_calculo_2018),
      ens20_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp6, "sex"=sexo_2020, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2020),
      ens21_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp6, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=High.LDL, "pdm"=pred_comp6, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.LDL.PRED6

R.LDL <- paste0("R.LDL.PRED",1:6) %>% lapply(get) %>% sapply(R.getOR2) %>% t %>% data.frame("Labs"=1:6) %>% 
  mutate("Outcome"=6, "Labs"=ordered(
    Labs, levels=6:1, labels=c("HbA1c (IEC)", "IFG (WHO)", "Both (ADA)", "Any (ADA)", "HbA1c (ADA)", "IFG (ADA)")))


#### Logistic - High CRP ####
# pred_comp1: FPG ≥100 (ADA) 
rbind(ens21_fin.2 %>% select("out"=High.PCR, "pdm"=pred_comp1, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=High.PCR, "pdm"=pred_comp1, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.PCR.PRED1

# pred_comp2: HbA1c ≥5.7 (ADA)
rbind(ens21_fin.2 %>% select("out"=High.PCR, "pdm"=pred_comp2, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=High.PCR, "pdm"=pred_comp2, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.PCR.PRED2

# pred_comp3: Any ADA criteria
rbind(ens21_fin.2 %>% select("out"=High.PCR, "pdm"=pred_comp3, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=High.PCR, "pdm"=pred_comp3, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.PCR.PRED3

# pred_comp4: Both ADA criteria
rbind(ens21_fin.2 %>% select("out"=High.PCR, "pdm"=pred_comp4, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=High.PCR, "pdm"=pred_comp4, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.PCR.PRED4

# pred_comp5: FPG ≥110 (WHO)
rbind(ens21_fin.2 %>% select("out"=High.PCR, "pdm"=pred_comp5, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=High.PCR, "pdm"=pred_comp5, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.PCR.PRED5

# pred_comp6: HbA1c ≥6.0 (IEC)
rbind(ens21_fin.2 %>% select("out"=High.PCR, "pdm"=pred_comp6, "sex"=sexo_2021, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2021),
      ens22_fin.2 %>% select("out"=High.PCR, "pdm"=pred_comp6, "sex"=sexo_2022, "age"=edad_num, YEAR, swts, "bmi"=imc_calculo_2022)) %>% 
  mutate("pdm"=factor(pdm, levels=0:2, labels=c("Normal", "Prediabetes", "Diabetes"))) %>% drop_na %>%
  glmer(formula = out~pdm+factor(sex)+scale(age)+scale(bmi)+(1|YEAR), family = binomial(), weights = swts) -> R.PCR.PRED6

R.PCR <- paste0("R.PCR.PRED",1:6) %>% lapply(get) %>% sapply(R.getOR2) %>% t %>% data.frame("Labs"=1:6) %>% 
  mutate("Outcome"=7, "Labs"=ordered(
    Labs, levels=6:1, labels=c("HbA1c (IEC)", "IFG (WHO)", "Both (ADA)", "Any (ADA)", "HbA1c (ADA)", "IFG (ADA)")))


#### Poisson - Prediabetes change by year ####
## Mexican population by year ##
conapo1 <- readxl::read_xlsx("Bases/base_municipios_final_datos_01.xlsx")
conapo2 <- readxl::read_xlsx("Bases/base_municipios_final_datos_02.xlsx")
conapo <- rbind(conapo1, conapo2)
conapo_fin <- conapo %>% filter(!EDAD_QUIN %in% c("pobm_00_04", "pobm_05_09", "pobm_10_14", "pobm_15_19")) %>%
  transmute("YEAR"=`A—O`, POB) %>% group_by(YEAR) %>% summarise(N=sum(POB)) %>% filter(YEAR%in%c(2016,2018,2020,2021,2022))
#Offset 1: total population ≥20 years
ens16_fin.2$POB <- conapo_fin$N[conapo_fin$YEAR==2016]
ens18_fin.2$POB <- conapo_fin$N[conapo_fin$YEAR==2018]
ens20_fin.2$POB <- conapo_fin$N[conapo_fin$YEAR==2020]
ens21_fin.2$POB <- conapo_fin$N[conapo_fin$YEAR==2021]
ens22_fin.2$POB <- conapo_fin$N[conapo_fin$YEAR==2022]
#Offset 2: total diabetes-free population ≥20 years
ens16_fin.2$POB2 <- conapo_fin$N[conapo_fin$YEAR==2016] - svytotal(~diabetes_fin,ensanut_2016_survey,na.rm = T)[1]
ens18_fin.2$POB2 <- conapo_fin$N[conapo_fin$YEAR==2018] - svytotal(~diabetes_fin,ensanut_2018_survey,na.rm = T)[1]
ens20_fin.2$POB2 <- conapo_fin$N[conapo_fin$YEAR==2020] - svytotal(~diabetes_fin,ensanut_2020_survey,na.rm = T)[1]
ens21_fin.2$POB2 <- conapo_fin$N[conapo_fin$YEAR==2021] - svytotal(~diabetes_fin,ensanut_2021_survey,na.rm = T)[1]
ens22_fin.2$POB2 <- conapo_fin$N[conapo_fin$YEAR==2022] - svytotal(~diabetes_fin,ensanut_2022_survey,na.rm = T)[1]

## Diabetes trends ##
rbind(ens16_fin.2 %>% transmute("dm"=diabetes_fin, "yrs"=YEAR-2016, "ent"=entidad, POB, swts), #Total
      ens18_fin.2 %>% transmute("dm"=diabetes_fin, "yrs"=YEAR-2016, "ent"=entidad, POB, swts),
      ens20_fin.2 %>% transmute("dm"=diabetes_fin, "yrs"=YEAR-2016, "ent"=entidad, POB, swts),
      ens21_fin.2 %>% transmute("dm"=diabetes_fin, "yrs"=YEAR-2016, "ent"=entidad, POB, swts),
      ens22_fin.2 %>% transmute("dm"=diabetes_fin, "yrs"=YEAR-2016, "ent"=entidad, POB, swts)) %>% drop_na %>%
  glm(formula = dm~yrs+offset(log(POB)), family = poisson(), weights = swts) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[2,c(1,4)] -> coef.00; model -> poisson.00
rbind(ens16_fin.2 %>% transmute("dm"=previous_diabetes_2016, "yrs"=YEAR-2016, "ent"=entidad, POB, swts), #Diagnosed
      ens18_fin.2 %>% transmute("dm"=previous_diabetes_2018, "yrs"=YEAR-2016, "ent"=entidad, POB, swts),
      ens20_fin.2 %>% transmute("dm"=previous_diabetes_2020, "yrs"=YEAR-2016, "ent"=entidad, POB, swts),
      ens21_fin.2 %>% transmute("dm"=previous_diabetes_2021, "yrs"=YEAR-2016, "ent"=entidad, POB, swts),
      ens22_fin.2 %>% transmute("dm"=previous_diabetes_2022, "yrs"=YEAR-2016, "ent"=entidad, POB, swts)) %>% drop_na %>%
  glm(formula = dm~yrs+offset(log(POB)), family = poisson(), weights = swts) -> poisson.00_d
rbind(ens16_fin.2 %>% transmute("dm"=undx_diabetes, "yrs"=YEAR-2016, "ent"=entidad, POB, swts), #Undiagnosed
      ens18_fin.2 %>% transmute("dm"=undx_diabetes, "yrs"=YEAR-2016, "ent"=entidad, POB, swts),
      ens20_fin.2 %>% transmute("dm"=undx_diabetes, "yrs"=YEAR-2016, "ent"=entidad, POB, swts),
      ens21_fin.2 %>% transmute("dm"=undx_diabetes, "yrs"=YEAR-2016, "ent"=entidad, POB, swts),
      ens22_fin.2 %>% transmute("dm"=undx_diabetes, "yrs"=YEAR-2016, "ent"=entidad, POB, swts)) %>% drop_na %>%
  glm(formula = dm~yrs+offset(log(POB)), family = poisson(), weights = swts) -> poisson.00_u

pois0_summ <- list(poisson.00, poisson.00_d, poisson.00_u) %>%
  lapply(jtools::summ, confint=T, digits=3) %>% lapply(with, coeftable)
pois0_summ[[1]][2,1:3] %>% exp; pois0_summ[[2]][2,1:3] %>% exp; pois0_summ[[3]][2,1:3] %>% exp


## Prediabetes trends ##
rbind(ens16_fin.2 %>% transmute("pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts), #ADA-IFG
      ens18_fin.2 %>% transmute("pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~yrs+offset(log(POB2)), family = poisson(), weights = swts) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[2,c(1,4)] -> coef.0B1; model -> poisson.01
rbind(ens16_fin.2 %>% transmute("pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts), #ADA-A1C
      ens18_fin.2 %>% transmute("pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~yrs+offset(log(POB2)), family = poisson(), weights = swts) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[2,c(1,4)] -> coef.0B2; model -> poisson.02
rbind(ens16_fin.2 %>% transmute("pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts), #ADA-ANY
      ens18_fin.2 %>% transmute("pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~yrs+offset(log(POB2)), family = poisson(), weights = swts) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[2,c(1,4)] -> coef.0B3; model -> poisson.03
rbind(ens16_fin.2 %>% transmute("pdm"=pred_comp4, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts), #ADA-BOTH
      ens18_fin.2 %>% transmute("pdm"=pred_comp4, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("pdm"=pred_comp4, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("pdm"=pred_comp4, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("pdm"=pred_comp4, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~yrs+offset(log(POB2)), family = poisson(), weights = swts) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[2,c(1,4)] -> coef.0B4; model -> poisson.04
rbind(ens16_fin.2 %>% transmute("pdm"=pred_comp5, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts), #WHO-IFG
      ens18_fin.2 %>% transmute("pdm"=pred_comp5, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("pdm"=pred_comp5, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("pdm"=pred_comp5, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("pdm"=pred_comp5, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~yrs+offset(log(POB2)), family = poisson(), weights = swts) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[2,c(1,4)] -> coef.0B5; model -> poisson.05
rbind(ens16_fin.2 %>% transmute("pdm"=pred_comp6, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts), #IEC-A1C
      ens18_fin.2 %>% transmute("pdm"=pred_comp6, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("pdm"=pred_comp6, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("pdm"=pred_comp6, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("pdm"=pred_comp6, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~yrs+offset(log(POB2)), family = poisson(), weights = swts) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[2,c(1,4)] -> coef.0B6; model -> poisson.06

V2 <- c(summary(poisson.00)$coefficients[2,1], summary(poisson.01)$coefficients[2,1],summary(poisson.02)$coefficients[2,1],summary(poisson.03)$coefficients[2,1],
        summary(poisson.04)$coefficients[2,1],summary(poisson.05)$coefficients[2,1],summary(poisson.06)$coefficients[2,1])
V3 <- rbind(confint(poisson.00)[2,],confint(poisson.01)[2,],confint(poisson.02)[2,],confint(poisson.03)[2,],
        confint(poisson.04)[2,],confint(poisson.05)[2,],confint(poisson.06)[2,])
data.frame(c("Diabetes","ADA-IFG","ADA-A1C","ADA-ANY","ADA-BOTH","WHO-IFG","IEC-A1C"), exp(V2) %>% round(4)) %>%
  cbind(exp(V3) %>% round(4)) %>% `names<-`(c("Definition","RR","Lower","Upper"))

## Age ##
rbind(ens16_fin.2 %>% transmute("fac"=edad_cat_2016 %>% as.factor, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=edad_cat_2018 %>% as.factor, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=edad_cat_2020 %>% as.factor, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=edad_cat_2021 %>% as.factor, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=edad_cat_2022 %>% as.factor, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[6,c(1,4)] -> coef.A1; model -> poisson.A1
rbind(ens16_fin.2 %>% transmute("fac"=edad_cat_2016 %>% as.factor, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=edad_cat_2018 %>% as.factor, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=edad_cat_2020 %>% as.factor, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=edad_cat_2021 %>% as.factor, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=edad_cat_2022 %>% as.factor, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[6,c(1,4)] -> coef.A2; model -> poisson.A2
rbind(ens16_fin.2 %>% transmute("fac"=edad_cat_2016 %>% as.factor, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=edad_cat_2018 %>% as.factor, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=edad_cat_2020 %>% as.factor, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=edad_cat_2021 %>% as.factor, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=edad_cat_2022 %>% as.factor, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[6,c(1,4)] -> coef.A3; model -> poisson.A3

## Sex ##
rbind(ens16_fin.2 %>% transmute("fac"=sexo_2016, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=sexo_2018, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=sexo_2020, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=sexo_2021, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=sexo_2022, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.B1; model -> poisson.B1
rbind(ens16_fin.2 %>% transmute("fac"=sexo_2016, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=sexo_2018, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=sexo_2020, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=sexo_2021, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=sexo_2022, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.B2; model -> poisson.B2
rbind(ens16_fin.2 %>% transmute("fac"=sexo_2016, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=sexo_2018, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=sexo_2020, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=sexo_2021, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=sexo_2022, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.B3; model -> poisson.B3

## BMI ##
rbind(ens16_fin.2 %>% transmute("fac"=imc_cat_2016 %>% as.factor, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=imc_cat_2018 %>% as.factor, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=imc_cat_2020 %>% as.factor, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=imc_cat_2021 %>% as.factor, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=imc_cat_2022 %>% as.factor, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[6,c(1,4)] -> coef.C1; model -> poisson.C1
rbind(ens16_fin.2 %>% transmute("fac"=imc_cat_2016 %>% as.factor, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=imc_cat_2018 %>% as.factor, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=imc_cat_2020 %>% as.factor, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=imc_cat_2021 %>% as.factor, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=imc_cat_2022 %>% as.factor, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[6,c(1,4)] -> coef.C2; model -> poisson.C2
rbind(ens16_fin.2 %>% transmute("fac"=imc_cat_2016 %>% as.factor, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=imc_cat_2018 %>% as.factor, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=imc_cat_2020 %>% as.factor, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=imc_cat_2021 %>% as.factor, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=imc_cat_2022 %>% as.factor, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[6,c(1,4)] -> coef.C3; model -> poisson.C3

## WC ##
rbind(ens16_fin.2 %>% transmute("fac"=ob_central, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=ob_central, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=ob_central, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=ob_central, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.D1; model -> poisson.D1
rbind(ens16_fin.2 %>% transmute("fac"=ob_central, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=ob_central, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=ob_central, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=ob_central, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.D2; model -> poisson.D2
rbind(ens16_fin.2 %>% transmute("fac"=ob_central, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=ob_central, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=ob_central, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=ob_central, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.D3; model -> poisson.D3

## Smoking ##
rbind(ens16_fin.2 %>% transmute("fac"=smoking %>% as.factor, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=smoking %>% as.factor, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=smoking %>% as.factor, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=smoking %>% as.factor, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[6,c(1,4)] -> coef.E1; model -> poisson.E1
rbind(ens16_fin.2 %>% transmute("fac"=smoking %>% as.factor, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=smoking %>% as.factor, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=smoking %>% as.factor, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=smoking %>% as.factor, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[6,c(1,4)] -> coef.E2; model -> poisson.E2
rbind(ens16_fin.2 %>% transmute("fac"=smoking %>% as.factor, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=smoking %>% as.factor, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=smoking %>% as.factor, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=smoking %>% as.factor, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[6,c(1,4)] -> coef.E3; model -> poisson.E3

## Indigenous identity ##
rbind(ens16_fin.2 %>% transmute("fac"=lengua_indigena, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=lengua_indigena, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=lengua_indigena, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=lengua_indigena, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=lengua_indigena, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.F1; model -> poisson.F1
rbind(ens16_fin.2 %>% transmute("fac"=lengua_indigena, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=lengua_indigena, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=lengua_indigena, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=lengua_indigena, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=lengua_indigena, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.F2; model -> poisson.F2
rbind(ens16_fin.2 %>% transmute("fac"=lengua_indigena, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=lengua_indigena, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=lengua_indigena, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=lengua_indigena, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=lengua_indigena, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.F3; model -> poisson.F3

## DISLI ##
rbind(ens16_fin.2 %>% transmute("fac"=DISLI_cat2, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=regiones_2016, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=DISLI_cat2, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=regiones_2018, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=DISLI_cat2, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=regiones_2020, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=DISLI_cat2, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=regiones_2021, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=DISLI_cat2, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=regiones_2022, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.G1; model -> poisson.G1
rbind(ens16_fin.2 %>% transmute("fac"=DISLI_cat2, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=regiones_2016, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=DISLI_cat2, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=regiones_2018, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=DISLI_cat2, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=regiones_2020, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=DISLI_cat2, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=regiones_2021, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=DISLI_cat2, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=regiones_2022, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.G2; model -> poisson.G2
rbind(ens16_fin.2 %>% transmute("fac"=DISLI_cat2, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=regiones_2016, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=DISLI_cat2, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=regiones_2018, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=DISLI_cat2, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=regiones_2020, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=DISLI_cat2, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=regiones_2021, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=DISLI_cat2, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=regiones_2022, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.G3; model -> poisson.G3

## Area ##
rbind(ens16_fin.2 %>% transmute("fac"=area2_2016, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=area_2018, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=area_2020, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=area2_2021, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=area2_2022, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.H1; model -> poisson.H1
rbind(ens16_fin.2 %>% transmute("fac"=area2_2016, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=area_2018, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=area_2020, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=area2_2021, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=area2_2022, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.H2; model -> poisson.H2
rbind(ens16_fin.2 %>% transmute("fac"=area2_2016, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=area_2018, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens20_fin.2 %>% transmute("fac"=area_2020, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=area2_2021, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=area2_2022, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.H3; model -> poisson.H3

## Family history ##
rbind(ens16_fin.2 %>% transmute("fac"=T2D.AHF, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=T2D.AHF, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=T2D.AHF, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=T2D.AHF, "pdm"=pred_comp1, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.I1; model -> poisson.I1
rbind(ens16_fin.2 %>% transmute("fac"=T2D.AHF, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=T2D.AHF, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=T2D.AHF, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=T2D.AHF, "pdm"=pred_comp2, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.I2; model -> poisson.I2
rbind(ens16_fin.2 %>% transmute("fac"=T2D.AHF, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens18_fin.2 %>% transmute("fac"=T2D.AHF, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens21_fin.2 %>% transmute("fac"=T2D.AHF, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts),
      ens22_fin.2 %>% transmute("fac"=T2D.AHF, "pdm"=pred_comp3, "yrs"=YEAR-2016, "ent"=entidad, POB2, swts)) %>%
  filter(pdm!=2) %>% drop_na %>% glm(formula = pdm~fac*yrs+offset(log(POB2)), family = poisson()) -> model; model %>%
  summary %>% coefficients() -> mcoef; mcoef[4,c(1,4)] -> coef.I3; model -> poisson.I3

#P-values for the interaction effect
paste0("Interaction\np = ",format(coef.A1[2], scientific = F, digits = 3) %>% toupper) -> PRL.A1
paste0("Interaction\np = ",format(coef.A2[2], scientific = F, digits = 3) %>% toupper) -> PRL.A2
paste0("Interaction\np = ",format(coef.A3[2], scientific = T, digits = 3) %>% toupper) -> PRL.A3

paste0("Interaction\np = ",format(coef.B1[2], scientific = F, digits = 3) %>% toupper) -> PRL.B1
paste0("Interaction\np = ",format(coef.B2[2], scientific = F, digits = 3) %>% toupper) -> PRL.B2
paste0("Interaction\np = ",format(coef.B3[2], scientific = F, digits = 3) %>% toupper) -> PRL.B3

paste0("Interaction\np = ",format(coef.C1[2], scientific = F, digits = 3) %>% toupper) -> PRL.C1
paste0("Interaction\np = ",format(coef.C2[2], scientific = T, digits = 3) %>% toupper) -> PRL.C2
paste0("Interaction\np = ",format(coef.C3[2], scientific = F, digits = 2) %>% toupper) -> PRL.C3

paste0("Interaction\np = ",format(coef.D1[2], scientific = F, digits = 3) %>% toupper) -> PRL.D1
paste0("Interaction\np = ",format(coef.D2[2], scientific = T, digits = 3) %>% toupper) -> PRL.D2
paste0("Interaction\np = ",format(coef.D3[2], scientific = F, digits = 2) %>% toupper) -> PRL.D3

paste0("Interaction\np = ",format(coef.E1[2], scientific = F, digits = 3) %>% toupper) -> PRL.E1
paste0("Interaction\np = ",format(coef.E2[2], scientific = F, digits = 3) %>% toupper) -> PRL.E2
paste0("Interaction\np = ",format(coef.E3[2], scientific = F, digits = 3) %>% toupper) -> PRL.E3

paste0("Interaction\np = ",format(coef.F1[2], scientific = F, digits = 2) %>% toupper) -> PRL.F1
paste0("Interaction\np = ",format(coef.F2[2], scientific = F, digits = 3) %>% toupper) -> PRL.F2
paste0("Interaction\np = ",format(coef.F3[2], scientific = F, digits = 3) %>% toupper) -> PRL.F3

paste0("Interaction\np = ",format(coef.G1[2], scientific = F, digits = 3) %>% toupper) -> PRL.G1
paste0("Interaction\np = ",format(coef.G2[2], scientific = F, digits = 3) %>% toupper) -> PRL.G2
paste0("Interaction\np = ",format(coef.G3[2], scientific = F, digits = 3) %>% toupper) -> PRL.G3

paste0("Interaction\np = ",format(coef.H1[2], scientific = F, digits = 3) %>% toupper) -> PRL.H1
paste0("Interaction\np = ",format(coef.H2[2], scientific = T, digits = 3) %>% toupper) -> PRL.H2
paste0("Interaction\np = ",format(coef.H3[2], scientific = T, digits = 3) %>% toupper) -> PRL.H3

paste0("Interaction\np = ",format(coef.I1[2], scientific = F, digits = 3) %>% toupper) -> PRL.I1
paste0("Interaction\np = ",format(coef.I2[2], scientific = T, digits = 3) %>% toupper) -> PRL.I2
paste0("Interaction\np = ",format(coef.I3[2], scientific = T, digits = 3) %>% toupper) -> PRL.I3

####------------------------------ FIGURES -----------------------------#### ----####
#### Figure 1: Disturbances in glucose metabolism ####
prev<-rbind(diabetes, prediabetes) %>% mutate(
  "prop"=round(prop,1),"IC95"=round(IC95,1), "lIC95"=round(lIC95,1),"uIC95"=round(uIC95,1),
  "cluster"=ordered(cluster, levels=c("Diabetes", "All", "None"), labels=c("Diabetes",  "Prediabetes", "None")))

fig1A <- rbind(diabetes, prediabetes) %>% 
  mutate("cluster"=ordered(cluster, levels=c("Diabetes", "All", "None"),
                           labels=c("Diabetes",  "Prediabetes", "None"))) %>% 
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=cluster,colour=cluster, linetype=cluster)) + 
  geom_line(size=1.5) + geom_point(size=2)+
  geom_errorbar(aes(ymin=lIC95, ymax=uIC95), width=.2,
                position=position_dodge(0.01), linetype=1)+ theme_pubclean()+
  ylab("Weighted prevalence (%)")+ xlab("ENSANUT cycle") + labs(colour="", linetype="") +
  scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) + ylim(8,34) +
  scale_color_manual(values=c("red4", "#552586")) + scale_linetype_manual(values=c(1,6)) +
  ggtitle("Disturbances in glucose metabolism") + theme(legend.position = "bottom") +
  theme(plot.title = element_text(size=15, face="bold", hjust=0.5, vjust=0)) +
  geom_text(data=prev %>% filter(cluster=="Prediabetes"), aes(label = paste0(
    prop[1:5],"%","\n","(",lIC95[1:5],"-",uIC95[1:5],")"), vjust=-1*((IC95/2)+0.7)), size=2.75)+
  geom_text(data=prev %>% filter(cluster=="Diabetes"), aes(
    label = paste0(prop,"%","\n","(",lIC95,"-",uIC95,")"), vjust=(IC95/2)+1.6),size=2.75)

fig1B <- rbind(prev_ifg_hba1c %>% mutate(x="C"), prev_disc_hba1c %>% mutate(x="B"),
               prev_disc_ifg %>% mutate(x="A")) %>%  select(YEAR, x, prop, IC95) %>%
  mutate(x=ordered(x, labels=c("IFG-Normal A1c","NFG-High A1c","IFG-High A1c"))) %>% 
  group_by(YEAR) %>% arrange(desc(x)) %>%
  mutate(pos = cumsum(prop), upper=pos+IC95, lower=pos-IC95) %>% ungroup() %>%
  ggplot(aes(x=YEAR, y=prop, fill=x), xLabels=NA) +
  geom_bar(stat="identity", color="black", linetype=2) + labs(fill="") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .1, col = "black") +
  theme_pubclean() + scale_x_discrete(limits = c("2016","2018","2020","2021","2022")) +
  xlab("ENSANUT cycle") + ylab ("Weighted prevalence (%)") +
  scale_fill_manual(values=c("#B589D6", "#804FB3", "#552586")) + ylim(0,32) +
  ggtitle("Discordance in prediabetes diagnosis") + theme(legend.position = "bottom") +
  theme(plot.title = element_text(size=15, face="bold", hjust=0.5, vjust=0)) +
  geom_text(aes(label = paste0(round(prop,1),"%","  (",round(prop,1)-round(IC95,1),"-", round(prop,1)+round(IC95,1),")")),
            size=2.7, col="white", position = position_stack(vjust = .5))

## HOMA2-IR ##
homa_16<- ensanut_2016 %>% dplyr::select(ID, HOMA2IR, HOMA2B, any_glucose) %>%
  mutate(year="2016"); names(homa_16)<-c("ID", "HOMA2IR", "HOMA2B","pred", "year")
homa_18<- ensanut_2018 %>% dplyr::select(ID, HOMA2IR, HOMA2B,any_glucose) %>%
  mutate(year="2018"); names(homa_18)<-c("ID", "HOMA2IR", "HOMA2B","pred","year")
homa_20<- ensanut_2020 %>% dplyr::select(ID, HOMA2IR, HOMA2B,any_glucose) %>%
  mutate(year="2020"); names(homa_20)<-c("ID", "HOMA2IR", "HOMA2B","pred","year")
homa_21<- ensanut_2021 %>% dplyr::select(ID, HOMA2IR, HOMA2B,any_glucose) %>%
  mutate(year="2021"); names(homa_21)<-c("ID", "HOMA2IR", "HOMA2B","pred","year")
homa_22<- ensanut_2022 %>% dplyr::select(ID, HOMA2IR, HOMA2B,any_glucose) %>%
  mutate(year="2022"); names(homa_22)<-c("ID", "HOMA2IR", "HOMA2B","pred","year")

homa_fin<-rbind(homa_16,homa_18, homa_20, homa_21, homa_22)
homa_fin[homa_fin=="Unable to calculate"]<-NA; homa_fin$HOMA2IR[homa_fin$HOMA2IR==0]<-NA
homa_fin$HOMA2IR<-as.numeric(homa_fin$HOMA2IR); homa_fin$HOMA2B<-as.numeric(homa_fin$HOMA2B)
labs<-c("Normoglycemic", "IFG-Normal HbA1c","NFG-High HbA1c","IFG-High HbA1c", "Diabetes")
homa_fin <- homa_fin %>% mutate("pred"=factor(pred, labels=labs), "homa_cat"=ifelse(HOMA2IR>=2.5,1,0))

fig1C <- homa_fin %>% mutate(homa_cat=factor(homa_cat, labels=c("No-IR", "IR"))) %>%
  drop_na() %>% group_by(year, pred, homa_cat) %>% summarise(n=n()) %>%
  ggplot(aes(fill=homa_cat, y=n, x=year)) + geom_bar(position="fill", stat="identity") +
  facet_wrap(~pred, ncol=5) + scale_y_continuous(labels = scales::percent) +
  labs(fill="Status") + ylab("Percent overall (%)") + xlab("ENSANUT cycle") +
  theme_pubclean() + theme(legend.position = "bottom") +
  scale_fill_manual(values=c("#C1549C", "#F0DF9D")) +
  ggtitle("Insulin resistance status") + theme(legend.position = "right") +
  theme(plot.title = element_text(size=15, face="bold", hjust=0.5, vjust=0)) +
  theme(strip.text = element_text(face = "bold.italic"))

figure1<-ggarrange(fig1A,fig1B, ncol=2, nrow=1, labels= letters[1:2]) %>%
  ggarrange(fig1C, ncol=1, nrow=2, labels=c("", "c"), heights = c(1, 0.75))
ggsave(figure1, file="Figures/Figure1.pdf", bg="transparent",
       width=32, height=25, units=c("cm"), dpi=600, limitsize = FALSE)

#### Figure 2: Modifiers of the prediabetes prevalence ####
#Figure limits
lim1 <- c(all_age$lIC95, all_sex$lIC95, all_imc$lIC95, all_obc$lIC95,
          all_inl$lIC95, all_aru$lIC95, all_sli$lIC95, all_ahf$lIC95) %>% min %>% floor()
lim2 <- c(all_age$uIC95, all_sex$uIC95, all_imc$uIC95, all_obc$uIC95,
          all_inl$uIC95, all_aru$uIC95, all_sli$uIC95, all_ahf$uIC95) %>% max %>% ceiling()

##-- Age --##
all_age %>% mutate("A"=c("Any criteria")) %>% 
  mutate(age_cat=factor(age_cat, labels = c("20-39", "40-59", ">=60"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(age_cat),colour=factor(age_cat))) + geom_line(size=1.5) +
  geom_point(size=2) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ggtitle("Age group (years)") +
  theme_pubclean() + scale_color_manual(values=c("#B498BF", "#CA4E46", "#353134")) + ylim(lim1, lim2) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5)) -> f3A.1
#annotate(geom = "text", x = 2018.5, y=35, label=PRL.A3)

##-- Sex --##
all_sex %>% mutate("A"=c("Any criteria")) %>% 
  mutate(sex=factor(sex, labels = c("Men", "Women"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(sex),colour=factor(sex))) + geom_line(size=1.5) +
  geom_point(size=2) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ggtitle("Sex") +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) + ylim(lim1, lim2) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5)) -> f3B.1

##-- BMI --##
all_imc %>% mutate("A"=c("Any criteria")) %>% 
  mutate(imc_cat=factor(imc_cat, labels = c("Normal\nweight", "Overweight", "Obesity"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(imc_cat),colour=factor(imc_cat))) + geom_line(size=1.5) +
  geom_point(size=2) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ggtitle("BMI category") +
  theme_pubclean() + scale_color_manual(values=c("#B498BF", "#CA4E46", "#353134")) + ylim(lim1, lim2) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5),
        legend.text = element_text(hjust=0.5)) -> f3C.1

##-- Central obesity --##
all_obc %>% mutate("A"=c("Any criteria")) %>% 
  mutate(ob_central=factor(ob_central, labels = c("Normal", "Central\nobesity"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(ob_central),colour=factor(ob_central))) + geom_line(size=1.5) +
  geom_point(size=2) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ggtitle("Waist circumference") +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) + ylim(lim1, lim2) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5),
        legend.text = element_text(hjust=0.5)) -> f3D.1

##-- Smoking --##
all_smo %>% mutate("A"=c("Any criteria")) %>% 
  mutate(smoking=factor(smoking, labels = c("Never\nsmoker", "Former\nsmoker", "Current\nsmoker"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(smoking),colour=factor(smoking))) + geom_line(size=1.5) +
  geom_point(size=2) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted Prevalence, (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ggtitle("Smoking status") +
  theme_pubclean() + scale_color_manual(values=c("#B498BF", "#CA4E46", "#353134")) + ylim(lim1, lim2) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5),
        legend.text = element_text(hjust=0.5)) -> f3E.1

##-- Indigenous language --##
all_inl %>% mutate("A"=c("Any criteria")) %>% 
  mutate(inl=factor(inl, labels = c("No", "Yes"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(inl),colour=factor(inl))) + geom_line(size=1.5) +
  geom_point(size=2) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ggtitle("Indigenous identity") +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) + ylim(lim1, lim2) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5)) -> f3F.1

##-- Density Independent Social Lag Index --##
all_sli %>% mutate("A"=c("Any criteria")) %>% 
  mutate(SLI=factor(SLI, labels = c("Low/Middle", "High"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(SLI),colour=factor(SLI))) + geom_line(size=1.5) +
  geom_point(size=2) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ggtitle("DISLI category") +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) + ylim(lim1, lim2) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5)) -> f3G.1

##-- Area --##
all_aru %>% mutate("A"=c("Any criteria")) %>% 
  mutate(area=factor(area, labels = c("Rural", "Urban"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(area),colour=factor(area))) + geom_line(size=1.5) +
  geom_point(size=2) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ggtitle("Area") +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) + ylim(lim1, lim2) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5)) -> f3H.1

##-- Family History --##
all_ahf %>% mutate("A"=c("Any criteria")) %>% 
  mutate(Fam=factor(T2D.AHF, labels = c("No", "Yes"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(Fam),colour=factor(Fam))) + geom_line(size=1.5) +
  geom_point(size=2) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ggtitle("Family history of diabetes") +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) + ylim(lim1, lim2) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5)) -> f3I.1

##-- Joint figure --##
#AGE (A), SEX (B), BMI (C), WC (D), SMO (E), IND (F), SLI(G), AREA(H)
figure2 <- ggarrange(f3A.1,f3B.1,f3C.1,f3D.1,"","","","",f3I.1,f3F.1,f3G.1,f3H.1, ncol=4, nrow=3,
                     labels=c(letters[1:4], rep("",4), letters[5:8]), heights=c(1,0.05,1))

ggsave(figure2, file="Figures/Figure2.pdf", bg="transparent", width=36, height=20.2,
       units=c("cm"), dpi=600, limitsize = FALSE)


#### Figure 3: Prediabetes and cardiometabolic outcomes ####
###ENSANUT 2016
svyby(~a1c_iec, by=~YEAR, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(a1c_iec*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster) -> iec_2016 
###ENSANUT 2018
svyby(~a1c_iec, by=~YEAR, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(a1c_iec*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster) -> iec_2018
###ENSANUT 2020
svyby(~a1c_iec, by=~YEAR, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(a1c_iec*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster) -> iec_2020
###ENSANUT 2021
svyby(~a1c_iec, by=~YEAR, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(a1c_iec*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster) -> iec_2021
###ENSANUT 2022
svyby(~a1c_iec, by=~YEAR, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(a1c_iec*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="HbA1c") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster) -> iec_2022
a1c_iec<-rbind(iec_2016,iec_2018,iec_2020,iec_2021,iec_2022) %>%
  as.data.frame(); a1c_iec$YEAR<-rownames(a1c_iec); a1c_iec

###ENSANUT 2016
svyby(~ifg_who, by=~YEAR, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(ifg_who*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster) -> who_2016
###ENSANUT 2018
svyby(~ifg_who, by=~YEAR, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(ifg_who*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster) -> who_2018
###ENSANUT 2020
svyby(~ifg_who, by=~YEAR, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(ifg_who*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster) -> who_2020
###ENSANUT 2021
svyby(~ifg_who, by=~YEAR, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(ifg_who*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster) -> who_2021
###ENSANUT 2022
svyby(~ifg_who, by=~YEAR, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(ifg_who*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IFG") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster) -> who_2022
ifg_who<-rbind(who_2016,who_2018,who_2020,who_2021,who_2022) %>%
  as.data.frame(); ifg_who$YEAR<-rownames(ifg_who); ifg_who

### Plots
rbind(prev_ifg, prev_hba1c, prev_ifg_hba1c, prediabetes) %>% mutate(
  cl2 = ordered(cluster, levels=c("IFG", "High HbA1c", "All", "High HbA1c-IFG"),
                labels=c("IFG", "HbA1c", "Any", "Both")), "A"=c("ADA criteria")) %>% 
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=cl2,colour=cl2)) + geom_line(size=1.5) + geom_point(size=2, shape=19) +
  facet_wrap(~A) + ylim(0, 32) + geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  theme_pubclean()+ ylab("Weighted prevalence (%)")+xlab("ENSANUT cycle") + labs(colour="") +
  scale_color_manual(values=c("#B589D6", "#804FB3", "#5D2E8D", "#2C1B3D")) + scale_x_continuous(breaks = c(2016,2018,2020,2021,2022)) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic")) -> fig2A.1
rbind(ifg_who, a1c_iec) %>% mutate(cl2 = ordered(
  cluster, levels=c("IFG", "HbA1c"), labels=c("IFG (WHO)", "HbA1c (IEC)")), "A"=c("Other criteria")) %>% 
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=cl2,colour=cl2)) +geom_line(size=1.5) + geom_point(size=2, shape=15) +
  facet_wrap(~A) + ylim(0, 32) + geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  theme_pubclean()+ ylab("Weighted prevalence (%)")+xlab("ENSANUT cycle")+ labs(colour="")+
  scale_color_manual(values=c("#E25E5E", "#9C2222")) + scale_x_continuous(breaks = c(2016,2018,2020,2021,2022)) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic")) -> fig2A.2
fig2A <- ggarrange(fig2A.1, fig2A.2, ncol=1, nrow=2, labels=letters[1:2])

rbind(R.HTA, R.HTC, R.HTG, R.IR, R.PCR, R.METS, R.CVD) %>% mutate("Outcome"=ordered(
  Outcome, levels=c(6:5,7,4:1), labels= c("Cardiovascular\nDisease","Metabolic\nSyndrome", "High CRP", "Insulin\nResistance",
                                 "Hyper-\ntriglyceridemia","Hyper-\ncholesterolemia", "Arterial\nHypertension"))) %>%
  ggplot(aes(x=OR, y=Outcome, color=Labs, shape=Labs)) + geom_vline(xintercept=1, linetype=2) +
  geom_pointrangeh(aes(xmin = L95, xmax = U95), position = position_dodge(width = 0.5), size = 0.5) +
  labs(y=NULL, x="Odds Ratio (log-scale)", color="Prediabetes\ndefinition", shape="Prediabetes\ndefinition") +
  theme_pubclean() + theme(legend.position = "bottom") + scale_x_log10() +
  scale_shape_manual(values=c(15,15,19,19,19,19), guide=guide_legend(reverse = T,  byrow=T, ncol=3)) +
  scale_color_manual(values=c("#9C2222", "#E25E5E", "#2C1B3D","#5D2E8D", "#804FB3", "#B589D6"),
                     guide=guide_legend(reverse = T, byrow=T, ncol=3)) + theme(legend.title.align = 0.5) -> fig2B

figure3<-ggarrange(fig2A, "", fig2B, ncol=3, nrow=1, widths = c(7,0.075,10), labels = c("","","c"))

ggsave(figure3, file="Figures/Figure3.pdf", bg="transparent",
       width=40*0.75, height=25*0.75, units=c("cm"), dpi=600, limitsize = FALSE)

#### Supp Figure 1: Participant selection ####
#Adults ≥20 years old
ensanut_2016 %>% nrow
ensanut_2018 %>% nrow
ensanut_2020 %>% nrow
ensanut_2021 %>% nrow
ensanut_2022 %>% nrow
nrow(ensanut_2016)+nrow(ensanut_2018)+nrow(ensanut_2020)+nrow(ensanut_2021)+nrow(ensanut_2022)
#Laboratory subset
ensanut_2016_2 %>% nrow
ensanut_2018_2 %>% nrow
ensanut_2020_2 %>% nrow
ensanut_2021_2 %>% nrow
ensanut_2022_2 %>% nrow
nrow(ensanut_2016_2)+nrow(ensanut_2018_2)+nrow(ensanut_2020_2)+nrow(ensanut_2021_2)+nrow(ensanut_2022_2)
#Survey design
ensanut_2016_survey %>% nrow
ensanut_2018_survey %>% nrow
ensanut_2020_survey %>% nrow
ensanut_2021_survey %>% nrow
ensanut_2022_survey %>% nrow
nrow(ensanut_2016_survey)+nrow(ensanut_2018_survey)+nrow(ensanut_2020_survey)+nrow(ensanut_2021_survey)+nrow(ensanut_2022_survey)

ensanut_2016_survey[["variables"]][["ponde_f_vv"]] %>% sum
ensanut_2018_survey[["variables"]][["ponderador_glucosa"]] %>% sum
ensanut_2020_survey[["variables"]][["ponde_g20.y"]] %>% sum
ensanut_2021_survey[["variables"]][["ponde_vv"]] %>% sum
ensanut_2022_survey[["variables"]][["ponde_v"]] %>% sum



#### Supp Figure 2: Missing data by year ####
names_supp1A <- c("Age","Sex","Indigenous","Area","BMI","WC","SBP","DBP","T2D Family Hx","Hx Diabetes", "Hx High BP", "Hx CVD")
names_supp1B <- c("Glucose","HbA1c","Insulin","Triglycerides","Total choleserol","HDL-cholesterol")

rbind(cbind("Year"=paste0(2016, " (n=", nrow(ensanut_2016), ")"), "Miss"=(ensanut_2016 %>% select(
  edad.x, sexo_2016, lengua_indigena, area2_2016, imc, cintura, PAS_fin, PAD_fin, T2D.AHF,
  previous_diabetes_2016, a401.x, CVD_fin) %>%  apply(2,is.na) %>%
    apply(2,sum)/nrow(ensanut_2016))) %>% cbind("Var"=names_supp1A),
  cbind("Year"=paste0(2018, " (n=", nrow(ensanut_2018), ")"), "Miss"=(ensanut_2018 %>% select(
    EDAD.x, sexo_2018, lengua_indigena, area_2018, imc_calculo_2018, waist_2018,
    PAS_fin, PAD_fin, T2D.AHF, previous_diabetes_2018, P4_1.x, CVD_fin) %>%  apply(2,is.na) %>%
      apply(2,sum)/nrow(ensanut_2018))) %>% cbind("Var"=names_supp1A),
  cbind("Year"=paste0(2020, " (n=", nrow(ensanut_2020), ")"), "Miss"=(ensanut_2020 %>% transmute(
    edad_num, sexo_2020, lengua_indigena, area_2020, imc_calculo_2020, "Waist"=NA,
    PAS_fin, PAD_fin, "T2D.AHF"=NA, previous_diabetes_2020, "HTA"=HX_HBP, CVD_fin) %>%  apply(2,is.na) %>%
      apply(2,sum)/nrow(ensanut_2020))) %>% cbind("Var"=names_supp1A),
  cbind("Year"=paste0(2021, " (n=", nrow(ensanut_2021), ")"), "Miss"=(ensanut_2021 %>% select(
    edad_num, sexo_2021, lengua_indigena, area2_2021, imc_calculo_2021, waist_2021,
    PAS_fin, PAD_fin, T2D.AHF, previous_diabetes_2021, a0401.x, CVD_fin) %>%  apply(2,is.na) %>%
      apply(2,sum)/nrow(ensanut_2021))) %>% cbind("Var"=names_supp1A),
  cbind("Year"=paste0(2022, " (n=", nrow(ensanut_2022), ")"), "Miss"=(ensanut_2022 %>% select(
    edad_num, sexo_2022, lengua_indigena, area2_2022, imc_calculo_2022, waist_2022,
    PAS_fin, PAD_fin, T2D.AHF, previous_diabetes_2022, a0401, CVD_fin) %>%  apply(2,is.na) %>%
      apply(2,sum)/nrow(ensanut_2022))) %>% cbind("Var"=names_supp1A)) %>% as.data.frame %>%
  mutate("Var"=ordered(Var,rev(names_supp1A),rev(names_supp1A)), "Miss"=as.numeric(Miss),
         "Lab"= (Miss*100) %>% round(1) %>% paste0("%")) %>% mutate(
           "Miss2"=ifelse(Miss<1, Miss, 0)*100, "Lab"=ifelse(Miss<1, Lab, "NA")) %>% 
  ggplot(aes(x=Var, y=Miss2)) + geom_col(position="dodge2", fill="#9C2222") + coord_flip(ylim=c(0,80)) +
  facet_wrap(~Year, ncol=5, nrow=1) + geom_text(aes(label=Lab, y=Miss2+1), size=3, fontface="bold.italic", hjust=0) +
  labs(y="Percentage of missingness", x=NULL, title="Adult subset") + theme_pubclean() +
  theme(strip.text = element_text(face = "bold.italic"), legend.position = "none",
        plot.title = element_text(hjust=0.5, face="bold"), strip.background = element_blank()) -> SuppF2A

rbind(cbind("Year"=paste0(2016, " (n=", nrow(ensanut_2016_2), ")"), "Miss"=(ensanut_2016_2 %>% select(
  valor.GLU_SUERO, valor.HB1AC, valor.INSULINA, valor.TRIG,
  valor.COLEST, valor.COL_HDL) %>%  apply(2,is.na) %>%
    apply(2,sum)/nrow(ensanut_2016_2))) %>% cbind("Var"=names_supp1B),
  cbind("Year"=paste0(2018, " (n=", nrow(ensanut_2018_2), ")"), "Miss"=(ensanut_2018_2 %>% select(
    VALOR_GLU_SUERO, VALOR_HB1AC, VALOR_INSULINA, VALOR_TRIG,
    VALOR_COLEST, VALOR_COL_HDL) %>%  apply(2,is.na) %>%
      apply(2,sum)/nrow(ensanut_2018_2))) %>% cbind("Var"=names_supp1B),
  cbind("Year"=paste0(2020, " (n=", nrow(ensanut_2020_2), ")"), "Miss"=(ensanut_2020_2 %>% transmute(
    valor.GLU_SUERO, HB1AC.Valor, valor.INSULINA, valor.TRIG,
    valor.COLEST, valor.COL_HDL) %>%  apply(2,is.na) %>%
      apply(2,sum)/nrow(ensanut_2020_2))) %>% cbind("Var"=names_supp1B),
  cbind("Year"=paste0(2021, " (n=", nrow(ensanut_2021_2), ")"), "Miss"=(ensanut_2021_2 %>% select(
    valor_GLU_SUERO, valor_HB1AC, valor_INSULINA, valor_TRIG,
    valor_COLEST, valor_COL_HDL) %>%  apply(2,is.na) %>%
      apply(2,sum)/nrow(ensanut_2021_2))) %>% cbind("Var"=names_supp1B),
  cbind("Year"=paste0(2022, " (n=", nrow(ensanut_2022_2), ")"), "Miss"=(ensanut_2022_2 %>% select(
    valor_GLU_SUERO, valor_HB1AC, valor_INSULINA, valor_TRIG,
    valor_COLEST, valor_COL_HDL) %>%  apply(2,is.na) %>%
      apply(2,sum)/nrow(ensanut_2022_2))) %>% cbind("Var"=names_supp1B)) %>% as.data.frame %>% 
  mutate("Var"=ordered(Var,rev(names_supp1B),rev(names_supp1B)), "Miss"=as.numeric(Miss),
         "Lab"= (Miss*100) %>% round(1) %>% paste0("%")) %>% mutate(
           "Miss2"=ifelse(Miss<1, Miss, 0)*100, "Lab"=ifelse(Miss<1, Lab, "NA")) %>% 
  ggplot(aes(x=Var, y=Miss2)) + geom_col(position="dodge2", fill="#9C2222") + coord_flip(ylim=c(0,0.0475*100)) +
  facet_wrap(~Year, ncol=5, nrow=1) + geom_text(aes(label=Lab, y=Miss2+0.001*100), size=3, fontface="bold.italic", hjust=0) +
  labs(y="Percentage of missingness", x=NULL, title="Laboratory subset") + theme_pubclean() +
  theme(strip.text = element_text(face = "bold.italic"), legend.position = "none",
        plot.title = element_text(hjust=0.5, face="bold"), strip.background = element_blank()) -> SuppF2B

supp_fig2 <- ggarrange(SuppF2A, SuppF2B, nrow=2, ncol=1, labels=letters[1:2])
ggsave(supp_fig2, file="Figures/Supp_fig2.pdf", bg="transparent",
       width=29, height=20.5, units=c("cm"), dpi=600, limitsize = FALSE)



#### Supp Figure 3: Modifiers (HbA1c and IFG) ####
#Figure limits
lim1 <- c(hba1c_age$lIC95, hba1c_sex$lIC95, hba1c_imc$lIC95, hba1c_obc$lIC95,
          hba1c_inl$lIC95, hba1c_aru$lIC95, hba1c_sli$lIC95, hba1c_ahf$lIC95,
          ifg_age$lIC95, ifg_sex$lIC95, ifg_imc$lIC95, ifg_obc$lIC95,
          ifg_inl$lIC95, ifg_aru$lIC95, ifg_sli$lIC95, ifg_ahf$lIC95) %>% min %>% floor()
lim2 <- c(hba1c_age$uIC95, hba1c_sex$uIC95, hba1c_imc$uIC95, hba1c_obc$uIC95,
          hba1c_inl$uIC95, hba1c_aru$uIC95, hba1c_sli$uIC95, hba1c_ahf$uIC95,
          ifg_age$uIC95, ifg_sex$uIC95, ifg_imc$uIC95, ifg_obc$uIC95,
          ifg_inl$uIC95, ifg_aru$uIC95, ifg_sli$uIC95, ifg_ahf$uIC95) %>% max %>% ceiling()

##-- Age --##
hba1c_age %>% mutate("A"=c("HbA1c criteria")) %>% 
  mutate(age_cat=factor(age_cat, labels = c("20-39", "40-59", ">=60"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(age_cat),colour=factor(age_cat))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#B498BF", "#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3A.2
ifg_age %>% mutate("A"=c("IFG criteria")) %>% 
  mutate(age_cat=factor(age_cat, labels = c("20-39", "40-59", ">=60"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(age_cat),colour=factor(age_cat))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#B498BF", "#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3A.3

##-- Sex --##
hba1c_sex %>% mutate("A"=c("HbA1c criteria")) %>% 
  mutate(sex=factor(sex, labels = c("Men", "Women"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(sex),colour=factor(sex))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3B.2
ifg_sex %>% mutate("A"=c("IFG criteria")) %>% 
  mutate(sex=factor(sex, labels = c("Men", "Women"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(sex),colour=factor(sex))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3B.3

##-- BMI --##
hba1c_imc %>% mutate("A"=c("HbA1c criteria")) %>% 
  mutate(imc_cat=factor(imc_cat, labels = c("Normal\nweight", "Overweight", "Obesity"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(imc_cat),colour=factor(imc_cat))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#B498BF", "#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3C.2
ifg_imc %>% mutate("A"=c("IFG criteria")) %>% 
  mutate(imc_cat=factor(imc_cat, labels = c("Normal\nweight", "Overweight", "Obesity"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(imc_cat),colour=factor(imc_cat))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#B498BF", "#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3C.3

##-- Central obesity --##
hba1c_obc %>% mutate("A"=c("HbA1c criteria")) %>% 
  mutate(ob_central=factor(ob_central, labels = c("Normal", "Central\nobesity"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(ob_central),colour=factor(ob_central))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3D.2
ifg_obc %>% mutate("A"=c("IFG criteria")) %>% 
  mutate(ob_central=factor(ob_central, labels = c("Normal", "Central\nobesity"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(ob_central),colour=factor(ob_central))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3D.3

##-- Smoking --##
hba1c_smo %>% mutate("A"=c("HbA1c criteria")) %>% 
  mutate(smoking=factor(smoking, labels = c("Never\nsmoker", "Former\nsmoker", "Current\nsmoker"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(smoking),colour=factor(smoking))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted Prevalence, (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#B498BF", "#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3E.2
ifg_smo %>% mutate("A"=c("IFG criteria")) %>% 
  mutate(smoking=factor(smoking, labels = c("Never\nsmoker", "Former\nsmoker", "Current\nsmoker"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(smoking),colour=factor(smoking))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted Prevalence, (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#B498BF", "#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3E.3

##-- Indigenous language --##
hba1c_inl %>% mutate("A"=c("HbA1c criteria")) %>% 
  mutate(inl=factor(inl, labels = c("No", "Yes"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(inl),colour=factor(inl))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3F.2
ifg_inl %>% mutate("A"=c("IFG criteria")) %>% 
  mutate(inl=factor(inl, labels = c("No", "Yes"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(inl),colour=factor(inl))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3F.3

##-- Social Lag Index --##
hba1c_sli %>% mutate("A"=c("HbA1c criteria")) %>% 
  mutate(SLI=factor(SLI, labels = c("Low/Middle", "High"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(SLI),colour=factor(SLI))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3G.2
ifg_sli %>% mutate("A"=c("IFG criteria")) %>% 
  mutate(SLI=factor(SLI, labels = c("Low/Middle", "High"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(SLI),colour=factor(SLI))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3G.3

##-- Area --##
hba1c_aru %>% mutate("A"=c("HbA1c criteria")) %>% 
  mutate(area=factor(area, labels = c("Rural", "Urban"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(area),colour=factor(area))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3H.2
ifg_aru %>% mutate("A"=c("IFG criteria")) %>% 
  mutate(area=factor(area, labels = c("Rural", "Urban"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(area),colour=factor(area))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3H.3

##-- Family History --##
hba1c_ahf %>% mutate("A"=c("HbA1c criteria")) %>% 
  mutate(Fam=factor(T2D.AHF, labels = c("No", "Yes"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(Fam),colour=factor(Fam))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3I.2
ifg_ahf %>% mutate("A"=c("IFG criteria")) %>% 
  mutate(Fam=factor(T2D.AHF, labels = c("No", "Yes"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(Fam),colour=factor(Fam))) + geom_line(size=1.5) +
  geom_point(size=2) + facet_wrap(~A) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ylim(lim1,lim2) +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold.italic"),
        legend.text = element_text(hjust=0.5)) -> f3I.3

##-- Joint figure --##
#AGE (A), SEX (B), BMI (C), WC (D), SMO (E), IND (F), SLI (G), AREA (H), FHX (I)
FS.3A <- ggarrange(f3A.2, f3A.3, nrow=1, ncol=2, common.legend = T, legend = "bottom") %>% 
  annotate_figure(top = text_grob("Prediabetes prevalence by age group (years)", color = "black", face = "bold", size = 14))
FS.3B <- ggarrange(f3B.2, f3B.3, nrow=1, ncol=2, common.legend = T, legend = "bottom") %>% 
  annotate_figure(top = text_grob("Prediabetes prevalence by sex", color = "black", face = "bold", size = 14))
FS.3C <- ggarrange(f3C.2, f3C.3, nrow=1, ncol=2, common.legend = T, legend = "bottom") %>% 
  annotate_figure(top = text_grob("Prediabetes prevalence by BMI category", color = "black", face = "bold", size = 14))
FS.3D <- ggarrange(f3D.2, f3D.3, nrow=1, ncol=2, common.legend = T, legend = "bottom") %>% 
  annotate_figure(top=text_grob("Prediabetes prevalence by waist circumference", color="black", face="bold", size=14))
FS.3E <- ggarrange(f3E.2, f3E.3, nrow=1, ncol=2, common.legend = T, legend = "bottom") %>% 
  annotate_figure(top=text_grob("Prediabetes prevalence by smoking status", color="black", face="bold", size=14))
FS.3F <- ggarrange(f3F.2, f3F.3, nrow=1, ncol=2, common.legend = T, legend = "bottom") %>% 
  annotate_figure(top = text_grob("Prediabetes prevalence by indigenous identity", color = "black", face = "bold", size = 14))
FS.3G <- ggarrange(f3G.2, f3G.3, nrow=1, ncol=2, common.legend = T, legend = "bottom") %>% 
  annotate_figure(top = text_grob("Prediabetes prevalence by DISLI category", color = "black", face = "bold", size = 14))
FS.3H <- ggarrange(f3H.2, f3H.3, nrow=1, ncol=2, common.legend = T, legend = "bottom") %>% 
  annotate_figure(top = text_grob("Prediabetes prevalence by area", color = "black", face = "bold", size = 14))
FS.3I <- ggarrange(f3I.2, f3I.3, nrow=1, ncol=2, common.legend = T, legend = "bottom") %>% 
  annotate_figure(top = text_grob("Prediabetes prevalence by family history of diabetes", color = "black", face = "bold", size = 14))

Supp_fig3 <- ggarrange(FS.3A,FS.3B,"","",FS.3C,FS.3D,"","",FS.3I,FS.3F,"","",FS.3G,FS.3H, ncol=2, nrow=7,
                       labels=c("a","b","","","c","d","","","e","f","","","g","h"), heights=c(1,1/20,1,1/20,1,1/20,1))
ggsave(Supp_fig3, file="Figures/Supp_fig3.pdf", bg="transparent", width=36, height=20.2*2,
       units=c("cm"), dpi=600, limitsize = FALSE)



#### Supp Figure 4: Indigenous population ####
##ENSANUT 2016
prediabetes_prev_2016<-svyby(
  ~homa_cat, by=~YEAR+lengua_indigena, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(homa_cat*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IR") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2018
prediabetes_prev_2018<-svyby(
  ~homa_cat, by=~YEAR+lengua_indigena, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(homa_cat*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IR") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2020
prediabetes_prev_2020<-svyby(
  ~homa_cat, by=~YEAR+lengua_indigena, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(homa_cat*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IR") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
##ENSANUT 2021
prediabetes_prev_2021<-svyby(
  ~homa_cat, by=~YEAR+lengua_indigena, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(homa_cat*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IR") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
##ENSANUT 2022
prediabetes_prev_2022<-svyby(
  ~homa_cat, by=~YEAR+lengua_indigena, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(homa_cat*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="IR") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
irp_inl <- rbind(prediabetes_prev_2016,prediabetes_prev_2018,prediabetes_prev_2020,
                 prediabetes_prev_2021,prediabetes_prev_2022)
irp_inl <- irp_inl %>% as.data.frame %>% mutate("inl"= c(rep(c(1,2),5)))

irp_inl %>% mutate("A"=c("Any criteria")) %>% 
  mutate(inl=factor(inl, labels = c("Non-indigenous", "Indigenous"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(inl),colour=factor(inl))) + geom_line(size=1.5) +
  geom_point(size=2) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ggtitle("Prevalence of insulin resistance") +
  theme_pubclean() + scale_color_manual(values=c("#CA4E46", "#353134")) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5)) -> suppF4.A

## HOMA2-IR ##
homa_16<- ensanut_2016 %>% dplyr::select(ID, HOMA2IR, HOMA2B, any_glucose, lengua_indigena) %>%
  mutate(year="2016"); names(homa_16)<-c("ID", "HOMA2IR", "HOMA2B","pred", "ind", "year")
homa_18<- ensanut_2018 %>% dplyr::select(ID, HOMA2IR, HOMA2B,any_glucose, lengua_indigena) %>%
  mutate(year="2018"); names(homa_18)<-c("ID", "HOMA2IR", "HOMA2B","pred", "ind", "year")
homa_20<- ensanut_2020 %>% dplyr::select(ID, HOMA2IR, HOMA2B,any_glucose, lengua_indigena) %>%
  mutate(year="2020"); names(homa_20)<-c("ID", "HOMA2IR", "HOMA2B","pred", "ind", "year")
homa_21<- ensanut_2021 %>% dplyr::select(ID, HOMA2IR, HOMA2B,any_glucose, lengua_indigena) %>%
  mutate(year="2021"); names(homa_21)<-c("ID", "HOMA2IR", "HOMA2B","pred", "ind", "year")
homa_22<- ensanut_2022 %>% dplyr::select(ID, HOMA2IR, HOMA2B,any_glucose, lengua_indigena) %>%
  mutate(year="2022"); names(homa_22)<-c("ID", "HOMA2IR", "HOMA2B","pred", "ind", "year")
homa_fin<-rbind(homa_16,homa_18, homa_20, homa_21, homa_22)
homa_fin[homa_fin=="Unable to calculate"]<-NA; homa_fin$HOMA2IR[homa_fin$HOMA2IR==0]<-NA
homa_fin$HOMA2IR<-as.numeric(homa_fin$HOMA2IR); homa_fin$HOMA2B<-as.numeric(homa_fin$HOMA2B)
labs<-c("Normoglycemic", "IFG-Normal HbA1c","NFG-High HbA1c","IFG-High HbA1c", "Diabetes")
homa_fin <- homa_fin %>% mutate("pred"=factor(pred, labels=labs), "homa_cat"=ifelse(HOMA2IR>=2.5,1,0))

np0 <- homa_fin %>% filter(ind==0) %>% drop_na %>% nrow; homa_fin %>%
  mutate(homa_cat=factor(homa_cat, labels=c("No-IR", "IR"))) %>% filter(ind==0) %>% 
  drop_na() %>% group_by(year, pred, homa_cat) %>% summarise(n=n()) %>%
  mutate(freq=n/sum(n)) %>% mutate(freq2=ifelse(homa_cat=="IR", round(freq*100,1), "")) %>% 
  ggplot(aes(fill=homa_cat, y=freq, x=year)) + geom_bar(position="fill", stat="identity") +
  facet_wrap(~pred, ncol=5) + scale_y_continuous(labels = scales::percent) +
  labs(fill="Status", y="Percent overall (%)", x="ENSANUT cycle", title=paste0("Non-indigenous (n = ", np0, ")")) +
  theme_pubclean() + theme(legend.position = "bottom") + scale_fill_manual(values=c("#C1549C", "#F0DF9D")) +
  geom_text(aes(label=freq2, y=freq, x=year), nudge_y = 0.04, color="black", fontface="bold.italic") +
  theme(legend.position = "right") + theme(plot.title = element_text(size=15, face="bold", hjust=0.5, vjust=0)) +
  theme(strip.text = element_text(face = "bold.italic")) -> suppF4.B1
np1 <- homa_fin %>% filter(ind==1) %>% drop_na %>% nrow; homa_fin %>%
  mutate(homa_cat=factor(homa_cat, labels=c("No-IR", "IR"))) %>% filter(ind==1) %>% 
  drop_na() %>% group_by(year, pred, homa_cat) %>% summarise(n=n()) %>%
  mutate(freq=n/sum(n)) %>% mutate(freq2=ifelse(homa_cat=="IR", round(freq*100,1), "")) %>% 
  ggplot(aes(fill=homa_cat, y=freq, x=year)) + geom_bar(position="fill", stat="identity") +
  facet_wrap(~pred, ncol=5) + scale_y_continuous(labels = scales::percent) +
  labs(fill="Status", y="Percent overall (%)", x="ENSANUT cycle", title=paste0("Indigenous (n = ", np1, ")")) +
  theme_pubclean() + theme(legend.position = "bottom") + scale_fill_manual(values=c("#C1549C", "#F0DF9D")) +
  geom_text(aes(label=freq2, y=freq, x=year), nudge_y = 0.04, color="black", fontface="bold.italic") +
  theme(legend.position = "right") + theme(plot.title = element_text(size=15, face="bold", hjust=0.5, vjust=0)) +
  theme(strip.text = element_text(face = "bold.italic")) -> suppF4.B2
suppF4.B <- ggarrange(suppF4.B1, suppF4.B2, nrow=2, ncol=1, common.legend = T, legend = "bottom")


all_indigenous <- rbind(
  ensanut_2016 %>% transmute("Age"=edad.x, "Sex"=sexo_2016, "DISLI"=DISLI, "BMI"=imc, "Waist"=cintura, "Area"=area2_2016,
                             "HOMA2IR"=HOMA2IR, "HOMA2B"=HOMA2B, "FPG"=valor.GLU_SUERO, "HbA1c"=valor.HB1AC, "Insulin"=valor.INSULINA,
                             "Triglycerides"=valor.TRIG, "Total_cholesterol"=valor.COLEST, "HDL_cholesterol"=valor.COL_HDL,
                             "YEAR"="2016", "Ind"=lengua_indigena %>% factor(levels=0:1, labels=c("Non-indigenous", "Indigenous"))),
  ensanut_2018 %>% transmute("Age"=EDAD.x, "Sex"=sexo_2018, "DISLI"=DISLI, "BMI"=imc_calculo_2018, "Waist"=waist_2018, "Area"=area_2018,
                             "HOMA2IR"=HOMA2IR, "HOMA2B"=HOMA2B, "FPG"=VALOR_GLU_SUERO, "HbA1c"=VALOR_HB1AC, "Insulin"=VALOR_INSULINA,
                             "Triglycerides"=VALOR_TRIG, "Total_cholesterol"=VALOR_COLEST, "HDL_cholesterol"=VALOR_COL_HDL,
                             "YEAR"="2018", "Ind"=lengua_indigena %>% factor(levels=0:1, labels=c("Non-indigenous", "Indigenous"))),
  ensanut_2020 %>% transmute("Age"=edad_num, "Sex"=sexo_2020, "DISLI"=DISLI, "BMI"=imc_calculo_2020, "Waist"=NA, "Area"=area_2020,
                             "HOMA2IR"=HOMA2IR, "HOMA2B"=HOMA2B, "FPG"=valor.GLU_SUERO, "HbA1c"=HB1AC.Valor, "Insulin"=valor.INSULINA,
                             "Triglycerides"=valor.TRIG, "Total_cholesterol"=valor.COLEST, "HDL_cholesterol"=valor.COL_HDL,
                             "YEAR"="2020", "Ind"=lengua_indigena %>% factor(levels=0:1, labels=c("Non-indigenous", "Indigenous"))),
  ensanut_2021 %>% transmute("Age"=edad_num, "Sex"=sexo_2021, "DISLI"=DISLI, "BMI"=imc_calculo_2021, "Waist"=waist_2021, "Area"=area2_2021,
                             "HOMA2IR"=HOMA2IR, "HOMA2B"=HOMA2B, "FPG"=valor_GLU_SUERO, "HbA1c"=valor_HB1AC, "Insulin"=valor_INSULINA,
                             "Triglycerides"=valor_TRIG, "Total_cholesterol"=valor_COLEST, "HDL_cholesterol"=valor_COL_HDL,
                             "YEAR"="2021", "Ind"=lengua_indigena %>% factor(levels=0:1, labels=c("Non-indigenous", "Indigenous"))),
  ensanut_2022 %>% transmute("Age"=edad_num, "Sex"=sexo_2022, "DISLI"=DISLI, "BMI"=imc_calculo_2022, "Waist"=waist_2022, "Area"=area2_2022,
                             "HOMA2IR"=HOMA2IR, "HOMA2B"=HOMA2B, "FPG"=valor_GLU_SUERO, "HbA1c"=valor_HB1AC, "Insulin"=valor_INSULINA,
                             "Triglycerides"=valor_TRIG, "Total_cholesterol"=valor_COLEST, "HDL_cholesterol"=valor_COL_HDL,
                             "YEAR"="2022", "Ind"=lengua_indigena %>% factor(levels=0:1, labels=c("Non-indigenous", "Indigenous")))) %>% filter(!is.na(Ind))

(all_indigenous %>% ggplot(aes(x=Ind, fill=Ind, y=Age))) %>% boxplot_indigenous -> suppF4.C1
(all_indigenous %>% ggplot(aes(x=Ind, fill=Ind, y=DISLI))) %>% boxplot_indigenous -> suppF4.C2
(all_indigenous %>% ggplot(aes(x=Ind, fill=Ind, y=BMI))) %>% boxplot_indigenous + ylim(0,100) -> suppF4.C3
(all_indigenous %>% ggplot(aes(x=Ind, fill=Ind, y=Waist))) %>% boxplot_indigenous + ylim(50,200) -> suppF4.C4
(all_indigenous %>% ggplot(aes(x=Ind, fill=Ind, y=HOMA2IR))) %>% boxplot_indigenous + ylim(0,10) -> suppF4.C5
(all_indigenous %>% ggplot(aes(x=Ind, fill=Ind, y=HOMA2B))) %>% boxplot_indigenous + ylim(20,200) -> suppF4.C6
(all_indigenous %>% ggplot(aes(x=Ind, fill=Ind, y=FPG))) %>% boxplot_indigenous + ylim(0,250) -> suppF4.C7
(all_indigenous %>% ggplot(aes(x=Ind, fill=Ind, y=HbA1c))) %>% boxplot_indigenous + ylim(3,15) -> suppF4.C8
(all_indigenous %>% ggplot(aes(x=Ind, fill=Ind, y=Insulin))) %>% boxplot_indigenous + ylim(0,55) -> suppF4.C9
(all_indigenous %>% ggplot(aes(x=Ind, fill=Ind, y=Triglycerides))) %>% boxplot_indigenous + ylim(0,1000) -> suppF4.C10
(all_indigenous %>% ggplot(aes(x=Ind, fill=Ind, y=Total_cholesterol))) %>% boxplot_indigenous + ylim(0,500) -> suppF4.C11
(all_indigenous %>% ggplot(aes(x=Ind, fill=Ind, y=HDL_cholesterol))) %>% boxplot_indigenous -> suppF4.C12

suppF4.C <- ggarrange(suppF4.C8, suppF4.C1, suppF4.C3, suppF4.C4, suppF4.C5, suppF4.C6, suppF4.C7, suppF4.C9,
                      suppF4.C10, suppF4.C11, suppF4.C12, suppF4.C2, nrow=2, ncol=6, common.legend = T, legend = "bottom")
suppF4 <- ggarrange(suppF4.A, "", suppF4.C, nrow=3, ncol=1, labels=c("a","","b"), heights = c(0.7,0.025,1))

ggsave(suppF4, file="Figures/Supp_fig4.pdf", bg="transparent", width=23.8*1.15, height=26.45*1.15,
       units=c("cm"), dpi=600, limitsize = FALSE)


####------------------------------ TABLES ------------------------------#### ----####
#### Main 1 (Table 2): Poisson regression (ADA-Any) ####
#Confidence intervals
list(poisson.00,poisson.01,poisson.02,poisson.03, poisson.04,poisson.05,poisson.06) %>%
  lapply(jtools::summ, confint=T, digits=3) %>% lapply(with, coeftable) -> OV.0
OV.0[[1]] -> OV.00; OV.0[[2]] -> OV.01; OV.0[[3]] -> OV.02; OV.0[[4]] -> OV.03; OV.0[[5]] -> OV.04; OV.0[[6]] -> OV.05; OV.0[[7]] -> OV.06

OV.00 <- OV.00[-1,c(1:3)] %>% exp %>% sprintf(fmt="%#.3f"); PCI.00 <- paste0("(RR ",OV.00[1],", 95%CI ",OV.00[2]," – ",OV.00[3],")")
OV.01 <- OV.01[-1,c(1:3)] %>% exp %>% sprintf(fmt="%#.3f"); PCI.01 <- paste0("(RR ",OV.01[1],", 95%CI ",OV.01[2]," – ",OV.01[3],")")
OV.02 <- OV.02[-1,c(1:3)] %>% exp %>% sprintf(fmt="%#.3f"); PCI.02 <- paste0("(RR ",OV.02[1],", 95%CI ",OV.02[2]," – ",OV.02[3],")")
OV.03 <- OV.03[-1,c(1:3)] %>% exp %>% sprintf(fmt="%#.3f"); PCI.03 <- paste0("(RR ",OV.03[1],", 95%CI ",OV.03[2]," – ",OV.03[3],")")
OV.04 <- OV.04[-1,c(1:3)] %>% exp %>% sprintf(fmt="%#.3f"); PCI.04 <- paste0("(RR ",OV.04[1],", 95%CI ",OV.04[2]," – ",OV.04[3],")")
OV.05 <- OV.05[-1,c(1:3)] %>% exp %>% sprintf(fmt="%#.3f"); PCI.05 <- paste0("(RR ",OV.05[1],", 95%CI ",OV.05[2]," – ",OV.05[3],")")
OV.06 <- OV.06[-1,c(1:3)] %>% exp %>% sprintf(fmt="%#.3f"); PCI.06 <- paste0("(RR ",OV.06[1],", 95%CI ",OV.06[2]," – ",OV.06[3],")")

#Change through the years (2016-2021)
PCI.00 #Diabetes
PCI.01 #ADA-IFG
PCI.02 #ADA-A1C
PCI.03 #ADA-ANY
PCI.04 #ADA-BOTH
PCI.05 #WHO-IFG
PCI.06 #IEC-A1C

PCI.A1 <- R.getIRR(poisson.A1); PCI.A2 <- R.getIRR(poisson.A2); PCI.A3 <- R.getIRR(poisson.A3)
PCI.B1 <- R.getIRR(poisson.B1); PCI.B2 <- R.getIRR(poisson.B2); PCI.B3 <- R.getIRR(poisson.B3)
PCI.C1 <- R.getIRR(poisson.C1); PCI.C2 <- R.getIRR(poisson.C2); PCI.C3 <- R.getIRR(poisson.C3)
PCI.D1 <- R.getIRR(poisson.D1); PCI.D2 <- R.getIRR(poisson.D2); PCI.D3 <- R.getIRR(poisson.D3)
PCI.I1 <- R.getIRR(poisson.I1); PCI.I2 <- R.getIRR(poisson.I2); PCI.I3 <- R.getIRR(poisson.I3)
PCI.F1 <- R.getIRR(poisson.F1); PCI.F2 <- R.getIRR(poisson.F2); PCI.F3 <- R.getIRR(poisson.F3)
PCI.G1 <- R.getIRR(poisson.G1); PCI.G2 <- R.getIRR(poisson.G2); PCI.G3 <- R.getIRR(poisson.G3)
PCI.H1 <- R.getIRR(poisson.H1); PCI.H2 <- R.getIRR(poisson.H2); PCI.H3 <- R.getIRR(poisson.H3)

c("20-39","40-59", ">=60", "ENSANUT year", "Year*Age 40-59", "Year*Age >=60") -> PLAB.A
c("Men","Women", "ENSANUT year", "Year*Women") -> PLAB.B
c("Normal weight","Overweight", "Obesity", "ENSANUT year", "Year*Overweight", "Year*Obesity") -> PLAB.C
c("Normal","Central obesity", "ENSANUT year", "Year*Obesity") -> PLAB.D
c("No","Yes", "ENSANUT year", "Year*Yes") -> PLAB.I
c("No","Yes", "ENSANUT year", "Year*Yes") -> PLAB.F
c("Low-Middle","High", "ENSANUT year", "Year*High") -> PLAB.G
c("Rural", "Urban", "ENSANUT year", "Year*Urban") -> PLAB.H

data.frame(
  c("Age (years)","","","","","", "Sex","","","", "BMI\ncategories","","","","","", "Waist\ncircumference","","","", "Family history\nof diabetes","","","",
    "Indigenous\nlanguage","","","", "DISLI","","","", "Area","","",""), c(PLAB.A, PLAB.B, PLAB.C, PLAB.D, PLAB.I, PLAB.F, PLAB.G, PLAB.H),
  c(PCI.A3, PCI.B3, PCI.C3, PCI.D3, PCI.I3, PCI.F3, PCI.G3, PCI.H3)) %>% `names<-`(c("Model", "Predictor", "RR (95% CI)")) %>%
  flextable(cwidth = c(2,1.5,1.5,2,1)) %>% align(align = "center",part = "all") %>% autofit() %>% save_as_docx(path="Tables/table2.docx")


#### Supp 1 (Table 1): Population characteristics ####
ensanut_2016$smoking3 <- factor(ensanut_2016$smoking, 0:2, c("Never smoker", "Former smoker", "Current smoker"))
ensanut_2018$smoking3 <- factor(ensanut_2018$smoking, 0:2, c("Never smoker", "Former smoker", "Current smoker"))
ensanut_2020$smoking3 <- factor(ensanut_2020$smoking, 0:2, c("Never smoker", "Former smoker", "Current smoker"))
ensanut_2021$smoking3 <- factor(ensanut_2021$smoking, 0:2, c("Never smoker", "Former smoker", "Current smoker"))
ensanut_2022$smoking3 <- factor(ensanut_2022$smoking, 0:2, c("Never smoker", "Former smoker", "Current smoker"))

ensanut_all <- rbind(
  ensanut_2016 %>% transmute(
    "Year"=2016, "Age"=edad.x, "Sex"=sexo_2016, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area2_2016, "AHF"=T2D.AHF, "BMI"=imc, "WC"=cintura,
    "Glu"=valor.GLU_SUERO, "A1c"=valor.HB1AC, "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=NA, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin),
  ensanut_2018 %>% transmute(
    "Year"=2018, "Age"=EDAD.x, "Sex"=sexo_2018, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area_2018, "AHF"=T2D.AHF, "BMI"=imc_calculo_2018, "WC"=waist_2018,
    "Glu"=VALOR_GLU_SUERO, "A1c"=VALOR_HB1AC,  "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=NA, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin),
  ensanut_2020 %>% transmute(
    "Year"=2020, "Age"=edad_num, "Sex"=sexo_2020, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area_2020, "AHF"=NA, "BMI"=imc_calculo_2020, "WC"=NA,
    "Glu"=valor.GLU_SUERO, "A1c"=HB1AC.Valor,  "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=NA, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin),
  ensanut_2021 %>% transmute(
    "Year"=2021, "Age"=edad_num, "Sex"=sexo_2021, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area2_2021, "AHF"=T2D.AHF, "BMI"=imc_calculo_2021, "WC"=waist_2021,
    "Glu"=valor_GLU_SUERO, "A1c"=valor_HB1AC,  "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=valor_PROTCREAC, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin),
  ensanut_2022 %>% transmute(
    "Year"=2022, "Age"=edad_num, "Sex"=sexo_2022, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area2_2022, "AHF"=T2D.AHF, "BMI"=imc_calculo_2022, "WC"=waist_2022,
    "Glu"=valor_GLU_SUERO, "A1c"=valor_HB1AC,  "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=valor_PCR, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin))

lab.in <- ensanut_all[-1] %>% names; lab.out <- c(
  "Age (years)", "Sex (women)", "Smoking status", "Indigenous language", "Urban area", "Family history of diabetes", "Body mass index (kg/m^2)",
  "Waist circumference (cm)", "Fasting plasma glucose (mg/dL)", "HbA1c (%)", "Total cholesterol (mg/dL)", "Triglycerides (mg/dL)", "CRP (mg/L)", "Hypertension",
  "High cholesterol or medication","High triglycerides or medication","Insulin resistance","Cardiovascular disease"); lab.all <- list(
    lab.in[1]~lab.out[1], lab.in[2]~lab.out[2], lab.in[3]~lab.out[3], lab.in[4]~lab.out[4], lab.in[5]~lab.out[5], lab.in[6]~lab.out[6],
    lab.in[7]~lab.out[7], lab.in[8]~lab.out[8], lab.in[9]~lab.out[9], lab.in[10]~lab.out[10], lab.in[11]~lab.out[11], lab.in[12]~lab.out[12],
    lab.in[13]~lab.out[13], lab.in[14]~lab.out[14], lab.in[15]~lab.out[15], lab.in[16]~lab.out[16],  lab.in[17]~lab.out[17],  lab.in[18]~lab.out[18])

gtsummary::tbl_summary(ensanut_all, by=Year, missing_text = "Missing", label = lab.all) %>% bold_labels() %>%
  modify_table_body(~.x %>% mutate(stat_1 = ifelse(stat_1%in%c("0 (NA%)", "NA (NA, NA)"), "-", stat_1),
                                   stat_2 = ifelse(stat_2%in%c("0 (NA%)", "NA (NA, NA)"), "-", stat_2),
                                   stat_3 = ifelse(stat_3%in%c("0 (NA%)", "NA (NA, NA)"), "-", stat_3))) %>% 
  as_flex_table() %>% align(align = "center", part = "all") %>% autofit() %>%
  save_as_docx(path="Tables/table1.docx", pr_section = prop_section(page_size = page_size(orient = "landscape")))

(ensanut_all$Sex %>% table %>% prop.table*100) %>% round(1)
(ensanut_all$Age %>% summary) %>% round(1)


#### Supp 2 (Table 7): Differences in population subsets ####
ensanut_2016_2<-ensanut_2016%>%filter(!is.na(ponde_f_vv), !is.na(est_var)) #F2
ensanut_2018_2 <- ensanut_2018 %>% filter(!is.na(ponderador_glucosa), !is.na(ESTRATO.y)) #F2
ensanut_2020_2 <- ensanut_2020 %>% filter(!is.na(ponde_g20.y), !is.na(est_sel.x)) %>% #F2
  filter(!est_sel.x%in%(which((est_sel.x%>%table)<2)%>%names%>%as.numeric))
ensanut_2021_2 <- ensanut_2021 %>% filter(!is.na(ponde_vv), !is.na(est_sel.lab)) %>% #F2
  filter(!est_sel.lab%in%(which((est_sel.lab%>%table)<2)%>%names%>%as.numeric)) 
ensanut_2022_2 <- ensanut_2022 %>% filter(!is.na(ponde_v), !is.na(est_sel.x.x.x)) #F2

tab.na.1 <- rbind(
  ensanut_2016 %>% transmute(
    "Sub"=1,"Year"=2016, "Age"=edad.x, "Sex"=sexo_2016, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area2_2016, "AHF"=T2D.AHF, "BMI"=imc, "WC"=cintura,
    "Glu"=valor.GLU_SUERO, "A1c"=valor.HB1AC, "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=NA, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin),
  ensanut_2018 %>% transmute(
    "Sub"=1,"Year"=2018, "Age"=EDAD.x, "Sex"=sexo_2018, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area_2018, "AHF"=T2D.AHF, "BMI"=imc_calculo_2018, "WC"=waist_2018,
    "Glu"=VALOR_GLU_SUERO, "A1c"=VALOR_HB1AC,  "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=NA, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin),
  ensanut_2020 %>% transmute(
    "Sub"=1,"Year"=2020, "Age"=edad_num, "Sex"=sexo_2020, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area_2020, "AHF"=NA, "BMI"=imc_calculo_2020, "WC"=NA,
    "Glu"=valor.GLU_SUERO, "A1c"=HB1AC.Valor,  "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=NA, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin),
  ensanut_2021 %>% transmute(
    "Sub"=1,"Year"=2021, "Age"=edad_num, "Sex"=sexo_2021, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area2_2021, "AHF"=T2D.AHF, "BMI"=imc_calculo_2021, "WC"=waist_2021,
    "Glu"=valor_GLU_SUERO, "A1c"=valor_HB1AC,  "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=valor_PROTCREAC, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin),
  ensanut_2022 %>% transmute(
    "Sub"=1,"Year"=2022, "Age"=edad_num, "Sex"=sexo_2022, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area2_2022, "AHF"=T2D.AHF, "BMI"=imc_calculo_2022, "WC"=waist_2022,
    "Glu"=valor_GLU_SUERO, "A1c"=valor_HB1AC,  "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=valor_PCR, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin))

tab.na.2 <- rbind(
  ensanut_2016_2 %>% transmute(
    "Sub"=2,"Year"=2016, "Age"=edad.x, "Sex"=sexo_2016, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area2_2016, "AHF"=T2D.AHF, "BMI"=imc, "WC"=cintura,
    "Glu"=valor.GLU_SUERO, "A1c"=valor.HB1AC, "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=NA, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin),
  ensanut_2018_2 %>% transmute(
    "Sub"=2,"Year"=2018, "Age"=EDAD.x, "Sex"=sexo_2018, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area_2018, "AHF"=T2D.AHF, "BMI"=imc_calculo_2018, "WC"=waist_2018,
    "Glu"=VALOR_GLU_SUERO, "A1c"=VALOR_HB1AC,  "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=NA, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin),
  ensanut_2020_2 %>% transmute(
    "Sub"=2,"Year"=2020, "Age"=edad_num, "Sex"=sexo_2020, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area_2020, "AHF"=NA, "BMI"=imc_calculo_2020, "WC"=NA,
    "Glu"=valor.GLU_SUERO, "A1c"=HB1AC.Valor,  "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=NA, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin),
  ensanut_2021_2 %>% transmute(
    "Sub"=2,"Year"=2021, "Age"=edad_num, "Sex"=sexo_2021, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area2_2021, "AHF"=T2D.AHF, "BMI"=imc_calculo_2021, "WC"=waist_2021,
    "Glu"=valor_GLU_SUERO, "A1c"=valor_HB1AC,  "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=valor_PROTCREAC, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin),
  ensanut_2022_2 %>% transmute(
    "Sub"=2,"Year"=2022, "Age"=edad_num, "Sex"=sexo_2022, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area2_2022, "AHF"=T2D.AHF, "BMI"=imc_calculo_2022, "WC"=waist_2022,
    "Glu"=valor_GLU_SUERO, "A1c"=valor_HB1AC,  "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=valor_PCR, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin))

tab.na.3 <- rbind(
  ensanut_2016_2 %>% filter((sanvenh>=8)&(a301.x!=2)&!is.na(valor.HB1AC)&!is.na(valor.GLU_SUERO)) %>% transmute(
    "Sub"=3,"Year"=2016, "Age"=edad.x, "Sex"=sexo_2016, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area2_2016, "AHF"=T2D.AHF, "BMI"=imc, "WC"=cintura,
    "Glu"=valor.GLU_SUERO, "A1c"=valor.HB1AC, "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=NA, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin),
  ensanut_2018_2 %>% filter((P5_1.y>=8)&(P3_1!=2)&!is.na(VALOR_HB1AC)&!is.na(VALOR_GLU_SUERO)) %>% transmute(
    "Sub"=3,"Year"=2018, "Age"=EDAD.x, "Sex"=sexo_2018, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area_2018, "AHF"=T2D.AHF, "BMI"=imc_calculo_2018, "WC"=waist_2018,
    "Glu"=VALOR_GLU_SUERO, "A1c"=VALOR_HB1AC,  "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=NA, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin),
  ensanut_2020_2 %>% filter((san04>=8)&!is.na(HB1AC.Valor)&!is.na(valor.GLU_SUERO)) %>% transmute(
    "Sub"=3,"Year"=2020, "Age"=edad_num, "Sex"=sexo_2020, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area_2020, "AHF"=NA, "BMI"=imc_calculo_2020, "WC"=NA,
    "Glu"=valor.GLU_SUERO, "A1c"=HB1AC.Valor,  "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=NA, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin),
  ensanut_2021_2 %>% filter((san04>=8)&(a0301.x!=2)&!is.na(valor_HB1AC)&!is.na(valor_GLU_SUERO)) %>% transmute(
    "Sub"=3,"Year"=2021, "Age"=edad_num, "Sex"=sexo_2021, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area2_2021, "AHF"=T2D.AHF, "BMI"=imc_calculo_2021, "WC"=waist_2021,
    "Glu"=valor_GLU_SUERO, "A1c"=valor_HB1AC,  "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=valor_PROTCREAC, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin),
  ensanut_2022_2 %>% filter((san04>=8)&(a0301!=2)&!is.na(valor_HB1AC)&!is.na(valor_GLU_SUERO)) %>% transmute(
    "Sub"=3,"Year"=2022, "Age"=edad_num, "Sex"=sexo_2022, "Smo"=smoking3, "Ind"=lengua_indigena, "Urb"=area2_2022, "AHF"=T2D.AHF, "BMI"=imc_calculo_2022, "WC"=waist_2022,
    "Glu"=valor_GLU_SUERO, "A1c"=valor_HB1AC,  "TC"=S.TC_fin, "TG"=S.TG_fin, "CRP"=valor_PCR, "HTN"=HTA_fin, "HCT"=HCT_fin, "HTG"=HTG_fin, "IR"=homa_cat, "CVD"=CVD_fin))

lab.in <- tab.na.1[-1] %>% names; lab.out <- c(
  "ENSANUT year", "Age (years)", "Sex (women)", "Smoking status", "Indigenous language", "Urban area", "Family history of diabetes", "Body mass index (kg/m^2)",
  "Waist circumference (cm)", "Fasting plasma glucose (mg/dL)", "HbA1c (%)", "Total cholesterol (mg/dL)", "Triglycerides (mg/dL)", "CRP (mg/L)", "Hypertension",
  "High cholesterol or medication","High triglycerides or medication","Insulin resistance","Cardiovascular disease"); lab.all <- list(
    lab.in[1]~lab.out[1], lab.in[2]~lab.out[2], lab.in[3]~lab.out[3], lab.in[4]~lab.out[4], lab.in[5]~lab.out[5], lab.in[6]~lab.out[6],
    lab.in[7]~lab.out[7], lab.in[8]~lab.out[8], lab.in[9]~lab.out[9], lab.in[10]~lab.out[10], lab.in[11]~lab.out[11], lab.in[12]~lab.out[12],
    lab.in[13]~lab.out[13], lab.in[14]~lab.out[14], lab.in[15]~lab.out[15], lab.in[16]~lab.out[16],  lab.in[17]~lab.out[17],  lab.in[18]~lab.out[18],  lab.in[19]~lab.out[19])


rbind(tab.na.1,tab.na.2,tab.na.3) %>% mutate("Sub"=ordered(Sub, 1:3, labels=c("Total sample","Laboratory sample", "Final sample"))) %>%
  gtsummary::tbl_summary(by=Sub, missing_text = "Missing", missing = "ifany", label = lab.all) %>% bold_labels() %>% 
  as_flex_table() %>% align(align = "center", part = "all") %>% autofit() %>%
  save_as_docx(path="Tables/table7.docx", pr_section = prop_section(page_size = page_size(orient = "landscape")))


#### Supp 3 (Table 5): Prevalence of each definition ####
rbind(
  rbind(diabetes[c(1,3:4)] %>% round(1), OV.00), rbind(prev_ifg[c(1,3:4)] %>% round(1), OV.01), rbind(prev_hba1c[c(1,3:4)] %>% round(1), OV.02),
  rbind(prediabetes[c(1,3:4)] %>% round(1), OV.03), rbind(prev_ifg_hba1c[c(1,3:4)] %>% round(1), OV.04),
  rbind(who_2016[c(1,3,4)] %>% round(1), who_2018[c(1,3,4)] %>% round(1), who_2020[c(1,3,4)] %>% round(1),who_2021[c(1,3,4)] %>% round(1), who_2022[c(1,3,4)] %>% round(1), OV.05),
  rbind(iec_2016[c(1,3,4)] %>% round(1), iec_2018[c(1,3,4)] %>% round(1), iec_2020[c(1,3,4)] %>% round(1), iec_2021[c(1,3,4)] %>% round(1), iec_2022[c(1,3,4)] %>% round(1),OV.06)
  ) %>% transmute( "D"=c(rep("Diabetes",6), rep("ADA-IFG",6), rep("ADA-A1c",6), rep("ADA-Any",6), rep("ADA-Both",6), rep("WHO-IFG",6), rep("IEC-A1c",6)),
                   "Y"=rep(c("2016","2018","2020","2021","2022","RR"), 7),"P"=prop,"C"=paste0(lIC95," – ",uIC95)) %>% `rownames<-`(NULL) %>%
  `colnames<-`(c("Definition","Year","Prevalence (%)", "95% CI")) %>% flextable(cwidth = c(2,1.5,1.5,2,1)) %>%
  align(align = "center",part = "all") %>% autofit() %>% save_as_docx(path="Tables/table5.docx")

#### Supp 4 (Table 3): Poisson regression (ADA-A1C) ####
data.frame(
  c("Age (years)","","","","","", "Sex","","","", "BMI\ncategories","","","","","", "Waist\ncircumference","","","", "Family history\nof diabetes","","","",
    "Indigenous\nlanguage","","","", "DISLI","","","", "Area","","",""), c(PLAB.A, PLAB.B, PLAB.C, PLAB.D, PLAB.I, PLAB.F, PLAB.G, PLAB.H),
  c(PCI.A2, PCI.B2, PCI.C2, PCI.D2, PCI.I2, PCI.F2, PCI.G2, PCI.H2)) %>% `names<-`(c("Model", "Predictor", "RR (95% CI)")) %>%
  flextable(cwidth = c(2,1.5,1.5,2,1)) %>% align(align = "center",part = "all") %>% autofit() %>% save_as_docx(path="Tables/table3.docx")

#### Supp 5 (Table 4): Poisson regression (ADA-IFG) ####
data.frame(
  c("Age (years)","","","","","", "Sex","","","", "BMI\ncategories","","","","","", "Waist\ncircumference","","","", "Family history\nof diabetes","","","",
    "Indigenous\nlanguage","","","", "DISLI","","","", "Area","","",""), c(PLAB.A, PLAB.B, PLAB.C, PLAB.D, PLAB.I, PLAB.F, PLAB.G, PLAB.H),
  c(PCI.A1, PCI.B1, PCI.C1, PCI.D1, PCI.I1, PCI.F1, PCI.G1, PCI.H1)) %>% `names<-`(c("Model", "Predictor", "RR (95% CI)")) %>%
  flextable(cwidth = c(2,1.5,1.5,2,1)) %>% align(align = "center",part = "all") %>% autofit() %>% save_as_docx(path="Tables/table4.docx")


#### Supp 6 (Table 6): Logistic regression ####
R.HTA2  <- R.HTA[1:3]  %>% rbind(R.HTA.PRED3  %>% R.getOR.DM2); R.HTA3  <- R.HTA2  %>% apply(2, sprintf, fmt="%#.3f")
R.HTC2  <- R.HTC[1:3]  %>% rbind(R.HTC.PRED3  %>% R.getOR.DM2); R.HTC3  <- R.HTC2  %>% apply(2, sprintf, fmt="%#.3f")
R.HTG2  <- R.HTG[1:3]  %>% rbind(R.HTG.PRED3  %>% R.getOR.DM2); R.HTG3  <- R.HTG2  %>% apply(2, sprintf, fmt="%#.3f")
R.IR2   <- R.IR[1:3]   %>% rbind(R.IR.PRED3   %>% R.getOR.DM2); R.IR3   <- R.IR2   %>% apply(2, sprintf, fmt="%#.3f")
R.PCR2  <- R.PCR[1:3]  %>% rbind(R.PCR.PRED3  %>% R.getOR.DM2); R.PCR3  <- R.PCR2  %>% apply(2, sprintf, fmt="%#.3f")
R.METS2 <- R.METS[1:3] %>% rbind(R.METS.PRED3 %>% R.getOR.DM2); R.METS3 <- R.METS2 %>% apply(2, sprintf, fmt="%#.3f")
R.CVD2  <- R.CVD[1:3]  %>% rbind(R.CVD.PRED3  %>% R.getOR.DM2); R.CVD3  <- R.CVD2  %>% apply(2, sprintf, fmt="%#.3f")

ST4.n <- c(R.HTA.PRED3@frame %>% nrow,R.HTC.PRED3@frame %>% nrow,R.HTG.PRED3@frame %>% nrow,R.IR.PRED3@frame %>% nrow,
           R.PCR.PRED3@frame %>% nrow,R.METS.PRED3@frame %>% nrow,R.CVD.PRED3@frame %>% nrow)
ST4.n2 <- paste0(c("Hypertension\nn=","Hypercholesterolemia\nn=","Hypertriglyceridemia\nn=",
                   "Insulin Resistance\nn=","High CRP\nn=","Metabolic Syndrome\nn=","Cardiovascular disease\nn="), ST4.n)
ST4.1 <- c(ST4.n2[1],"","","","","","", ST4.n2[2],"","","","","","", ST4.n2[3],"","","","","","",
           ST4.n2[4],"","","","","","", ST4.n2[5],"","","","","","", ST4.n2[6],"","","","","","",ST4.n2[7],"","","","","","")
ST4.2 <- rep(c("IFG (ADA)","HbA1c (ADA)","Any (ADA)","Both (ADA)","IFG (WHO)","HbA1c (IEC)","Diabetes mellitus"),7)
ST4.3 <- rbind(R.HTA3, R.HTC3, R.HTG3, R.IR3, R.PCR3, R.METS3, R.CVD3) %>% as.data.frame %>%  transmute("OR2"=paste0(OR," (",L95," – ",U95,")"))

cbind(ST4.1,ST4.2,ST4.3) %>% `names<-`(c("Outcome", "Predictor", "OR (95% CI)")) %>% flextable(cwidth = c(2,1.5,1.5,2,1)) %>%
  align(align = "center",part = "all") %>% autofit() %>% save_as_docx(path="Tables/table6.docx")




#### Drafts  ####




# Prevalence from all subgroups
rbind(prev_disc_ifg,prev_disc_hba1c,prev_ifg_hba1c,diabetes) %>% 
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=cluster,colour=cluster)) + 
  geom_line(size=1.5) + geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  theme_pubclean()+ ylab("Weighted Prevalence, (%)")+xlab("ENSANUT cycle")+ labs(colour="")+
  ggsci::scale_color_futurama()+ theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))+
  ggtitle("Disturbances in glucose metabolism")
rbind(prev_hba1c, prev_ifg, homa_cat) %>% 
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=cluster,colour=cluster)) + 
  geom_line(size=1.5) + geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  theme_pubclean()+ ylab("Weighted Prevalence, (%)")+ xlab("ENSANUT cycle")+ labs(colour="") +
  ggsci::scale_color_futurama()+ theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))+
  ggtitle("Disturbances in glucose metabolism")

# HOMA comparison
homa_2016_2<- ensanut_2016 %>% dplyr::select(ID, HOMA2IR, HOMA2B, any_glucose) %>%
  mutate(year="2016"); names(homa_2016_2)<-c("ID", "HOMA2IR", "HOMA2B","pred", "year")
homa_2018_2<- ensanut_2018 %>% dplyr::select(ID, HOMA2IR, HOMA2B,any_glucose) %>%
  mutate(year="2018"); names(homa_2018_2)<-c("ID", "HOMA2IR", "HOMA2B","pred","year")
homa_2020_2<- ensanut_2020 %>% dplyr::select(ID, HOMA2IR, HOMA2B,any_glucose) %>%
  mutate(year="2020"); names(homa_2020_2)<-c("ID", "HOMA2IR", "HOMA2B","pred","year")
homa_2021_2<- ensanut_2021 %>% dplyr::select(ID, HOMA2IR, HOMA2B,any_glucose) %>%
  mutate(year="2021"); names(homa_2021_2)<-c("ID", "HOMA2IR", "HOMA2B","pred","year")

homa_fin<-rbind(homa_2016_2,homa_2018_2, homa_2020_2, homa_2021_2)
homa_fin[homa_fin=="Unable to calculate"]<-NA; homa_fin$HOMA2IR[homa_fin$HOMA2IR==0]<-NA
homa_fin$HOMA2IR<-as.numeric(homa_fin$HOMA2IR); homa_fin$HOMA2B<-as.numeric(homa_fin$HOMA2B)
labs<-c("Normoglycemic", "IFG-Normal HbA1c", "High HbA1c-NFG", "IFG-High HbA1c", "Diabetes")
homa_fin$pred<-factor(homa_fin$pred, labels = labs); homa_fin$homa_cat<-ifelse(homa_fin$HOMA2IR>=2.5,1,0)

homa_fin %>% filter(pred!="Diabetes") %>%
  ggplot(aes(x=year, y=HOMA2IR, fill=pred))+geom_boxplot()+ facet_wrap(~pred)+theme_classic()+
  scale_y_log10()+ xlab("ENSANUT cycle")+ylab("HOMA2-IR")+theme(legend.position = "top")
homa_fin %>%filter(pred!="Diabetes") %>%
  ggplot(aes(x=year, y=HOMA2B, fill=pred))+ geom_boxplot()+ facet_wrap(~pred)+theme_classic()+
  xlab("ENSANUT cycle")+ylab("HOMA2-Beta (%)")+ theme(legend.position = "top")

# Figure 2???
homa<- homa_fin %>% drop_na() %>% group_by(year, pred) %>%
  summarise(n=n()); homa$pred<-ordered(
    homa$pred, levels=c("Normoglycemic", "IFG-Normal HbA1c", "High HbA1c-NFG", "IFG-High HbA1c", "Diabetes"),
    labels=c("Normoglycemic", "Prediabetes\nIFG", "Prediabetes\nHigh HbA1c", "Prediabetes\nBoth", "Diabetes"))
homa %>% drop_na() %>% 
  ggplot(aes(fill=pred, y=n, x=year)) + geom_bar(position="fill", stat="identity")+
  scale_y_continuous(labels = scales::percent)+  theme_classic() + labs(fill="Status")+
  ylab("Percent (%)")+xlab("ENSANUT year")+ theme(legend.position = "top", legend.text.align = 0.5)

# Alternate figure 1B
rbind(prev_ifg_hba1c %>% mutate(x="C"), prev_disc_hba1c %>% mutate(x="B"),
      prev_disc_ifg %>% mutate(x="A")) %>%  select(YEAR, x, prop, lIC95, uIC95) %>%
  mutate(x=ordered(x, labels=c("IFG-Normal HbA1c","NFG-High HbA1c","IFG-High HbA1c"))) %>% 
  ggplot(aes(x=YEAR, y=prop, fill=x), xLabels=NA) +
  geom_bar(stat="identity", color="black", linetype=2) + labs(fill="") +
  geom_errorbar(aes(ymin = lIC95, ymax = uIC95), width = .1, col = "black") +
  theme_pubclean() + scale_x_discrete(limits = c("2016","2018","2020","2021")) +
  xlab("ENSANUT cycle") + ylab ("Weighted prevalence (%)") + facet_wrap(~x, nrow = 3, scales = "fixed") +
  scale_fill_manual(values=c("#B589D6", "#804FB3", "#552586")) +
  theme(strip.text = element_text(face = "bold.italic")) +
  ggtitle("Prediabetes criteria") + theme(legend.position = "none") +
  theme(plot.title = element_text(size=15, face="bold", hjust=0.5, vjust=0))

#### HOMA stratification
#ENSANUT 2016
prediabetes_prevalence_2016_year<-svyby(~homa_cat, by=~YEAR+imc_cat_2016, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(homa_cat*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Insulin Resistance") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster, YEAR)
prediabetes_prevalence_2016_year
#ENSANUT 2018
prediabetes_prevalence_2018_year<-svyby(~homa_cat, by=~YEAR+imc_cat_2018, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(homa_cat*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Insulin Resistance") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
prediabetes_prevalence_2018_year
#ENSANUT 2020
prediabetes_prevalence_2020_year<-svyby(~homa_cat, by=~YEAR+imc_cat_2020, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(homa_cat*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Insulin Resistance") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
prediabetes_prevalence_2020_year
#ENSANUT 2021
prediabetes_prevalence_2021_year<-svyby(~homa_cat, by=~YEAR+imc_cat_2021, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(homa_cat*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95,cluster="Insulin Resistance") %>%
  dplyr::select(prop,IC95,lIC95,uIC95, cluster,YEAR)
prediabetes_prevalence_2021_year

homa_imc<-rbind(prediabetes_prevalence_2016_year,prediabetes_prevalence_2018_year,prediabetes_prevalence_2020_year,prediabetes_prevalence_2021_year)
homa_imc<-as.data.frame(homa_imc)
homa_imc$imc_cat<-c(rep(c(1,2,3),4))
homa_imc %>% 
  mutate(imc_cat=factor(imc_cat, labels = c("Normal weight", "Overweight", "Obesity"))) %>%
  ggplot(aes(x=as.numeric(YEAR), y=prop,group=factor(imc_cat),colour=factor(imc_cat))) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2,position=position_dodge(0.01))+
  theme_pubclean()+
  ylab("Weighted Prevalence, (%)")+
  xlab("ENSANUT cycle")+
  labs(colour="BMI category")+
  ggsci::scale_color_futurama()+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))+
  ggtitle("Insulin Resistance by HOMA-IR >2.5")


#### DEFINITIONS WITH PREVIOUS DIAGNOSIS
##-- Arterial hypertension (DX | TX | ≥140/90) --##
#2016
ensanut_2016$PAS_fin <- ensanut_2016 %>% select(sistol3.y, sistol4.y) %>% apply(1, mean, na.rm=T)
ensanut_2016$PAD_fin <- ensanut_2016 %>% select(diastol3.y, diastol4.y) %>% apply(1, mean, na.rm=T)
ensanut_2016$HTA_fin <- (with(ensanut_2016, (
  a401.x==1 | (a405.y==1&!is.na(a405.y)) | PAS_fin>=140&!is.na(PAS_fin) |
    PAD_fin>=90&!is.na(PAD_fin) ))) %>% ifelse(1,0); ensanut_2016$HTA_fin %>% table(useNA = "always")
#2018
ensanut_2018$P27_1_1[ensanut_2018$P27_1_1<80]<-NA; ensanut_2018$P27_2_1[ensanut_2018$P27_2_1<80]<-NA
ensanut_2018$P27_1_2[ensanut_2018$P27_1_2<49]<-NA; ensanut_2018$P27_2_2[ensanut_2018$P27_2_2<49]<-NA
ensanut_2018$P27_1_2[ensanut_2018$P27_1_2==222]<-NA
ensanut_2018$PAS_fin <- ensanut_2018 %>% select(P27_1_1, P27_2_1) %>% apply(1, mean, na.rm=T)
ensanut_2018$PAD_fin <- ensanut_2018 %>% select(P27_1_2, P27_2_2) %>% apply(1, mean, na.rm=T)
ensanut_2018$HTA_fin <- (with(ensanut_2018, (
  P4_1.x==1 | (P4_4==1&!is.na(P4_4)) | (PAS_fin>=140&!is.na(PAS_fin)) |
    (PAD_fin>=90&!is.na(PAD_fin)) ))) %>% ifelse(1,0); ensanut_2018$HTA_fin %>% table(useNA = "always")
#2020 (only ≥140/90)
ensanut_2020$an08_01s[ensanut_2020$an08_01s==999]<-NA; ensanut_2020$an08_01d[ensanut_2020$an08_01d==999]<-NA
ensanut_2020$an08_02s[ensanut_2020$an08_02s==999]<-NA; ensanut_2020$an08_02d[ensanut_2020$an08_02d==999]<-NA
ensanut_2020$an08_03s[ensanut_2020$an08_03s==999]<-NA; ensanut_2020$an08_03d[ensanut_2020$an08_03d==999]<-NA
ensanut_2020$PAS_fin <- ensanut_2020 %>% select(an08_01s, an08_02s, an08_03s) %>% apply(1, mean, na.rm=T)
ensanut_2020$PAD_fin <- ensanut_2020 %>% select(an08_01d, an08_02d, an08_03d) %>% apply(1, mean, na.rm=T)
ensanut_2020$HTA_fin <- (with(ensanut_2020, (PAS_fin>=140&!is.na(PAS_fin) | PAD_fin>=90&!is.na(PAD_fin)))) %>%
  ifelse(1,0); ensanut_2020$HTA_fin[with(ensanut_2020, is.na(PAS_fin)&is.na(PAD_fin))] <- NA
ensanut_2020$HTA_fin %>% table(useNA = "always")
#2021
ensanut_2021 <- ensanut_2021 %>% mutate("a0401.x"=ifelse(a0401==3,0,1))
ensanut_2021$an27_01s[ensanut_2021$an27_01s==999]<-NA; ensanut_2021$an27_01d[ensanut_2021$an27_01d==999]<-NA
ensanut_2021$an27_02s[ensanut_2021$an27_02s==999]<-NA; ensanut_2021$an27_02d[ensanut_2021$an27_02d==999]<-NA
ensanut_2021$an27_03s[ensanut_2021$an27_03s==999]<-NA; ensanut_2021$an27_03d[ensanut_2021$an27_03d==999]<-NA
ensanut_2021$PAS_fin <- ensanut_2021 %>% select(an27_01s, an27_02s, an27_03s) %>% apply(1, mean, na.rm=T)
ensanut_2021$PAD_fin <- ensanut_2021 %>% select(an27_01d, an27_02d, an27_03d) %>% apply(1, mean, na.rm=T)
ensanut_2021$HTA_fin <- (with(ensanut_2021, (
  a0401.x==1 | (a0404==1&!is.na(a0404)) | PAS_fin>=140&!is.na(PAS_fin) |
    PAD_fin>=90&!is.na(PAD_fin) ))) %>% ifelse(1,0); ensanut_2021$HTA_fin %>% table(useNA = "always")

##-- Hypercholesterolemia (DX | high CT) --##
#2016
ensanut_2016$a607.x[ensanut_2016$a607.x==3] <- NA; ensanut_2016$a607.x[ensanut_2016$a600a==2] <- 2
ensanut_2016$S.TC_fin <- ensanut_2016$valor.COLEST; ensanut_2016$S.TC_fin[ensanut_2016$sanvenh<8]<-NA
ensanut_2016$HCT_fin <- (with(ensanut_2016, (a607.x==1 | (S.TC_fin>=200&!is.na(S.TC_fin)) ))) %>%
  ifelse(1,0); ensanut_2016$HCT_fin %>% table(useNA = "always")
#2018
ensanut_2018$S.TC_fin <- ensanut_2018$VALOR_COLEST; ensanut_2018$S.TC_fin[ensanut_2018$P5_1.y<8]<-NA
ensanut_2018$HCT_fin <- (with(ensanut_2018, (P6_4==1 | (S.TC_fin>=200&!is.na(S.TC_fin)) ))) %>%
  ifelse(1,0); ensanut_2018$HCT_fin %>% table(useNA = "always")
#2020 (only labs)
ensanut_2020$S.TC_fin <- ensanut_2020$valor.COLEST; ensanut_2020$S.TC_fin[ensanut_2020$san04<8]<-NA
ensanut_2020$HCT_fin <- (with(ensanut_2020, ((S.TC_fin>=200&!is.na(S.TC_fin))))) %>%
  ifelse(1,0); ensanut_2020$HCT_fin[with(ensanut_2020, is.na(S.TC_fin))] <- NA
ensanut_2020$HCT_fin %>% table(useNA = "always")
#2021
ensanut_2021$S.TC_fin <- ensanut_2021$valor_COLEST; ensanut_2021$S.TC_fin[ensanut_2021$san04<8]<-NA
ensanut_2021$HCT_fin <- (with(ensanut_2021, (a0604==1 | (S.TC_fin>=200&!is.na(S.TC_fin)) ))) %>%
  ifelse(1,0); ensanut_2021$HCT_fin %>% table(useNA = "always")

##-- Hypertriglyceridemia (DX | high TG) --##
#2016
ensanut_2016$a609.x[ensanut_2016$a609.x==3] <- NA; ensanut_2016$a609.x[ensanut_2016$a600a==2] <- 2
ensanut_2016$S.TG_fin <- ensanut_2016$valor.TRIG; ensanut_2016$S.TG_fin[ensanut_2016$sanvenh<8]<-NA
ensanut_2016$HTG_fin <- (with(ensanut_2016, (a609.x==1 | (S.TG_fin>=150&!is.na(S.TG_fin)) ))) %>%
  ifelse(1,0); ensanut_2016$HTG_fin %>% table(useNA = "always")
#2018
ensanut_2018$S.TG_fin <- ensanut_2018$VALOR_TRIG; ensanut_2018$S.TG_fin[ensanut_2018$P5_1.y<8]<-NA
ensanut_2018$HTG_fin <- (with(ensanut_2018, (P6_6==1 | (S.TG_fin>=150&!is.na(S.TG_fin)) ))) %>%
  ifelse(1,0); ensanut_2018$HTG_fin %>% table(useNA = "always")
#2020 (only labs)
ensanut_2020$S.TG_fin <- ensanut_2020$valor.TRIG; ensanut_2020$S.TG_fin[ensanut_2020$san04<8]<-NA
ensanut_2020$HTG_fin <- (with(ensanut_2020, ((S.TG_fin>=150&!is.na(S.TG_fin))))) %>%
  ifelse(1,0); ensanut_2020$HTG_fin[with(ensanut_2020, is.na(S.TG_fin))] <- NA
ensanut_2020$HTG_fin %>% table(useNA = "always")
#2021
ensanut_2021$S.TG_fin <- ensanut_2021$valor_TRIG; ensanut_2021$S.TG_fin[ensanut_2021$san04<8]<-NA
ensanut_2021$HTG_fin <- (with(ensanut_2021, (a0606==1 | (S.TG_fin>=150&!is.na(S.TG_fin)) ))) %>%
  ifelse(1,0); ensanut_2021$HTG_fin %>% table(useNA = "always")



#### JOINT PLOTS FOR ANY, IFG AND HBA1C
fig3A <- ggarrange(f3A.1, f3A.2, f3A.3, nrow=1, ncol=3, common.legend = T, legend = "bottom") %>% 
  annotate_figure(top = text_grob("Prediabetes prevalence by age", color = "black", face = "bold", size = 14))
fig3B <- ggarrange(f3B.1, f3B.2, f3B.3, nrow=1, ncol=3, common.legend = T, legend = "bottom") %>% 
  annotate_figure(top = text_grob("Prediabetes prevalence by sex", color = "black", face = "bold", size = 14))
fig3C <- ggarrange(f3C.1, f3C.2, f3C.3, nrow=1, ncol=3, common.legend = T, legend = "bottom") %>% 
  annotate_figure(top = text_grob("Prediabetes prevalence by BMI", color = "black", face = "bold", size = 14))
fig3D <- ggarrange(f3D.1, f3D.2, f3D.3, nrow=1, ncol=3, common.legend = T, legend = "bottom") %>% 
  annotate_figure(top=text_grob("Prediabetes prevalence by waist circumference", color="black", face="bold", size=14))
fig3E <- ggarrange(f3E.1, f3E.2, f3E.3, nrow=1, ncol=3, common.legend = T, legend = "bottom") %>% 
  annotate_figure(top=text_grob("Prediabetes prevalence by smoking status", color="black", face="bold", size=14))
fig3F <- ggarrange(f3F.1, f3F.2, f3F.3, nrow=1, ncol=3, common.legend = T, legend = "bottom") %>% 
  annotate_figure(top = text_grob("Prediabetes prevalence by Indigenous identity", color = "black", face = "bold", size = 14))
fig3G <- ggarrange(f3G.1, f3G.2, f3G.3, nrow=1, ncol=3, common.legend = T, legend = "bottom") %>% 
  annotate_figure(top = text_grob("Prediabetes prevalence by DISLI category", color = "black", face = "bold", size = 14))
fig3H <- ggarrange(f3H.1, f3H.2, f3H.3, nrow=1, ncol=3, common.legend = T, legend = "bottom") %>% 
  annotate_figure(top = text_grob("Prediabetes prevalence by area", color = "black", face = "bold", size = 14))




#### Supp Fig 5 ####
##-- Geographical shapes --##
load("Bases/df_mx.rda")
geom_mx <-  sf::st_read(dsn="Bases/shapes", layer="areas_geoestadisticas_estatales")  %>%
  sf::st_transform(crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",stringsAsFactors=FALSE)
geom_mx<-rmapshaper::ms_simplify(geom_mx, keep = 0.01, keep_shapes = T)
geom_mx$CVE_ENT_NUM<-as.numeric(geom_mx$CVE_ENT)
geo_data<-spdep::poly2nb(geom_mx)
#Regions of Mexico (9)
geom_mx$region2<-geom_mx$CVE_ENT_NUM
geom_mx$region2[geom_mx$CVE_ENT_NUM %in% c(2,3,18,25,26)]<-1 #Pacific-North
geom_mx$region2[geom_mx$CVE_ENT_NUM %in% c(5,8,19,28)]<-2 #Frontier
geom_mx$region2[geom_mx$CVE_ENT_NUM %in% c(6,14,16)]<-3 #Pacific-Center
geom_mx$region2[geom_mx$CVE_ENT_NUM %in% c(1,10,11,22,24,32)]<-4 #Center-North
geom_mx$region2[geom_mx$CVE_ENT_NUM %in% c(13,29,30)]<-5 #Center
geom_mx$region2[geom_mx$CVE_ENT_NUM %in% c(9)]<-6 #Mexico City
geom_mx$region2[geom_mx$CVE_ENT_NUM %in% c(15)]<-7 #State of Mexico
geom_mx$region2[geom_mx$CVE_ENT_NUM %in% c(12,17,20,21)]<-8 #Pacific-South
geom_mx$region2[geom_mx$CVE_ENT_NUM %in% c(4,7,23,27,31)]<-9 #Peninsula


##-- Difference in prevalence (2020-2021) BY STATE --##
geom_mx$entidad<-geom_mx$CVE_ENT_NUM
## All criteria ##
all_sta2 <- all_sta %>% mutate(year2=YEAR-2016) %>% group_by(entidad) %>% filter(YEAR==2020|YEAR==2021)
all_sta2_diff <- all_sta2 %>% summarise(diff=prop[YEAR==2021]) %>% left_join(geom_mx,by="entidad") %>%
  left_join((SLI %>% select(DISLI,entidad)), by="entidad") %>% bi_class(x=DISLI, y=diff, style="quantile", dim=2)
## Both criteria ##
bth_sta2 <- bth_sta %>% mutate(year2=YEAR-2016) %>% group_by(entidad) %>% filter(YEAR==2020|YEAR==2021)
bth_sta2_diff <- bth_sta2 %>% summarise(diff=prop[YEAR==2021]) %>% left_join(geom_mx,by="entidad") %>%
  left_join((SLI %>% select(DISLI,entidad)), by="entidad") %>% bi_class(x=DISLI, y=diff, style="quantile", dim=2)
## High HbA1c ##
a1c_sta2 <- a1c_sta %>% mutate(year2=YEAR-2016) %>% group_by(entidad) %>% filter(YEAR==2020|YEAR==2021)
a1c_sta2_diff <- a1c_sta2 %>% summarise(diff=prop[YEAR==2021]) %>% left_join(geom_mx,by="entidad") %>%
  left_join((SLI %>% select(DISLI,entidad)), by="entidad") %>% bi_class(x=DISLI, y=diff, style="quantile", dim=2)
## High FPG ##
ifg_sta2 <- ifg_sta %>% mutate(year2=YEAR-2016) %>% group_by(entidad) %>% filter(YEAR==2020|YEAR==2021)
ifg_sta2_diff <- ifg_sta2 %>% summarise(diff=prop[YEAR==2021]) %>% left_join(geom_mx,by="entidad") %>%
  left_join((SLI %>% select(DISLI,entidad)), by="entidad") %>% bi_class(x=DISLI, y=diff, style="quantile", dim=2)


##-- Difference in prevalence (2020-2021) BY REGION (9) --##
## All criteria ##
all_reg2 <- all_reg9 %>% mutate(year2=YEAR-2016) %>% group_by(region2) %>% filter(YEAR==2020|YEAR==2021)
all_reg2_diff <- all_reg2 %>% summarise(diff=prop[YEAR==2021]) %>% left_join(geom_mx,by="region2") %>%
  left_join((SLI %>% select(DISLI,region2)), by="region2") %>% bi_class(x=DISLI, y=diff, style="quantile", dim=2)
allreg_fin <- all_sta2_diff; allreg_fin$diff[allreg_fin$region2==1]<-(all_reg2_diff$diff[all_reg2_diff$region2==1])[1]
allreg_fin$diff[allreg_fin$region2==2]<-(all_reg2_diff$diff[all_reg2_diff$region2==2])[1]
allreg_fin$diff[allreg_fin$region2==3]<-(all_reg2_diff$diff[all_reg2_diff$region2==3])[1]
allreg_fin$diff[allreg_fin$region2==4]<-(all_reg2_diff$diff[all_reg2_diff$region2==4])[1]
allreg_fin$diff[allreg_fin$region2==5]<-(all_reg2_diff$diff[all_reg2_diff$region2==5])[1]
allreg_fin$diff[allreg_fin$region2==6]<-(all_reg2_diff$diff[all_reg2_diff$region2==6])[1]
allreg_fin$diff[allreg_fin$region2==7]<-(all_reg2_diff$diff[all_reg2_diff$region2==7])[1]
allreg_fin$diff[allreg_fin$region2==8]<-(all_reg2_diff$diff[all_reg2_diff$region2==8])[1]
allreg_fin$diff[allreg_fin$region2==9]<-(all_reg2_diff$diff[all_reg2_diff$region2==9])[1]
## Both criteria ##
bth_reg2 <- bth_reg9 %>% mutate(year2=YEAR-2016) %>% group_by(region2) %>% filter(YEAR==2020|YEAR==2021)
bth_reg2_diff <- bth_reg2 %>% summarise(diff=prop[YEAR==2021]) %>% left_join(geom_mx,by="region2") %>%
  left_join((SLI %>% select(DISLI,region2)), by="region2") %>% bi_class(x=DISLI, y=diff, style="quantile", dim=2)
bthreg_fin <- bth_sta2_diff; bthreg_fin$diff[bthreg_fin$region2==1]<-(bth_reg2_diff$diff[bth_reg2_diff$region2==1])[1]
bthreg_fin$diff[bthreg_fin$region2==2]<-(bth_reg2_diff$diff[bth_reg2_diff$region2==2])[1]
bthreg_fin$diff[bthreg_fin$region2==3]<-(bth_reg2_diff$diff[bth_reg2_diff$region2==3])[1]
bthreg_fin$diff[bthreg_fin$region2==4]<-(bth_reg2_diff$diff[bth_reg2_diff$region2==4])[1]
bthreg_fin$diff[bthreg_fin$region2==5]<-(bth_reg2_diff$diff[bth_reg2_diff$region2==5])[1]
bthreg_fin$diff[bthreg_fin$region2==6]<-(bth_reg2_diff$diff[bth_reg2_diff$region2==6])[1]
bthreg_fin$diff[bthreg_fin$region2==7]<-(bth_reg2_diff$diff[bth_reg2_diff$region2==7])[1]
bthreg_fin$diff[bthreg_fin$region2==8]<-(bth_reg2_diff$diff[bth_reg2_diff$region2==8])[1]
bthreg_fin$diff[bthreg_fin$region2==9]<-(bth_reg2_diff$diff[bth_reg2_diff$region2==9])[1]
## High HbA1c ##
a1c_reg2 <- a1c_reg9 %>% mutate(year2=YEAR-2016) %>% group_by(region2) %>% filter(YEAR==2020|YEAR==2021)
a1c_reg2_diff <- a1c_reg2 %>% summarise(diff=prop[YEAR==2021]) %>% left_join(geom_mx,by="region2") %>%
  left_join((SLI %>% select(DISLI,region2)), by="region2") %>% bi_class(x=DISLI, y=diff, style="quantile", dim=2)
a1creg_fin <- a1c_sta2_diff; a1creg_fin$diff[a1creg_fin$region2==1]<-(a1c_reg2_diff$diff[a1c_reg2_diff$region2==1])[1]
a1creg_fin$diff[a1creg_fin$region2==2]<-(a1c_reg2_diff$diff[a1c_reg2_diff$region2==2])[1]
a1creg_fin$diff[a1creg_fin$region2==3]<-(a1c_reg2_diff$diff[a1c_reg2_diff$region2==3])[1]
a1creg_fin$diff[a1creg_fin$region2==4]<-(a1c_reg2_diff$diff[a1c_reg2_diff$region2==4])[1]
a1creg_fin$diff[a1creg_fin$region2==5]<-(a1c_reg2_diff$diff[a1c_reg2_diff$region2==5])[1]
a1creg_fin$diff[a1creg_fin$region2==6]<-(a1c_reg2_diff$diff[a1c_reg2_diff$region2==6])[1]
a1creg_fin$diff[a1creg_fin$region2==7]<-(a1c_reg2_diff$diff[a1c_reg2_diff$region2==7])[1]
a1creg_fin$diff[a1creg_fin$region2==8]<-(a1c_reg2_diff$diff[a1c_reg2_diff$region2==8])[1]
a1creg_fin$diff[a1creg_fin$region2==9]<-(a1c_reg2_diff$diff[a1c_reg2_diff$region2==9])[1]
## High FPG ##
ifg_reg2 <- ifg_reg9 %>% mutate(year2=YEAR-2016) %>% group_by(region2) %>% filter(YEAR==2020|YEAR==2021)
ifg_reg2_diff <- ifg_reg2 %>% summarise(diff=prop[YEAR==2021]) %>% left_join(geom_mx,by="region2") %>%
  left_join((SLI %>% select(DISLI,region2)), by="region2") %>% bi_class(x=DISLI, y=diff, style="quantile", dim=2)
ifgreg_fin <- ifg_sta2_diff; ifgreg_fin$diff[ifgreg_fin$region2==1]<-(ifg_reg2_diff$diff[ifg_reg2_diff$region2==1])[1]
ifgreg_fin$diff[ifgreg_fin$region2==2]<-(ifg_reg2_diff$diff[ifg_reg2_diff$region2==2])[1]
ifgreg_fin$diff[ifgreg_fin$region2==3]<-(ifg_reg2_diff$diff[ifg_reg2_diff$region2==3])[1]
ifgreg_fin$diff[ifgreg_fin$region2==4]<-(ifg_reg2_diff$diff[ifg_reg2_diff$region2==4])[1]
ifgreg_fin$diff[ifgreg_fin$region2==5]<-(ifg_reg2_diff$diff[ifg_reg2_diff$region2==5])[1]
ifgreg_fin$diff[ifgreg_fin$region2==6]<-(ifg_reg2_diff$diff[ifg_reg2_diff$region2==6])[1]
ifgreg_fin$diff[ifgreg_fin$region2==7]<-(ifg_reg2_diff$diff[ifg_reg2_diff$region2==7])[1]
ifgreg_fin$diff[ifgreg_fin$region2==8]<-(ifg_reg2_diff$diff[ifg_reg2_diff$region2==8])[1]
ifgreg_fin$diff[ifgreg_fin$region2==9]<-(ifg_reg2_diff$diff[ifg_reg2_diff$region2==9])[1]


##-- Global Lee's L --##
spdep::lee.test(ifgreg_fin$DISLI, ifgreg_fin$diff, nb2listw(geo_data),
                zero.policy = TRUE, alternative = "two.sided", na.action = na.omit)


##-- Plots --##
## By region ##
all_reg2_diff %>% ggplot() +
  geom_sf(mapping = aes(fill = bi_class, geometry=geometry), color = "white", size = 0.1, show.legend = F) +
  bi_scale_fill(pal = "GrPink", dim = 2) + labs(title = "Any criteria") + bi_theme() +
  theme(plot.title = element_text(size=15)) -> map; legend <- bi_legend(
    pal="GrPink", dim=2, xlab="DISLI", ylab="Change", size=12.5
  ); ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.7, .55, 0.25, 0.25) -> GGG_B.1
a1c_reg2_diff %>% ggplot() +
  geom_sf(mapping = aes(fill = bi_class, geometry=geometry), color = "white", size = 0.1, show.legend = F) +
  bi_scale_fill(pal = "GrPink", dim = 2) + labs(title = "HbA1c criteria") + bi_theme() +
  theme(plot.title = element_text(size=15)) -> map; legend <- bi_legend(
    pal="GrPink", dim=2, xlab="DISLI", ylab="Change", size=12.5
  ); ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.7, .55, 0.25, 0.25) -> GGG_B.2
ifg_reg2_diff %>% ggplot() +
  geom_sf(mapping = aes(fill = bi_class, geometry=geometry), color = "white", size = 0.1, show.legend = F) +
  bi_scale_fill(pal = "GrPink", dim = 2) + labs(title = "IFG criteria") + bi_theme() +
  theme(plot.title = element_text(size=15)) -> map; legend <- bi_legend(
    pal="GrPink", dim=2, xlab="DISLI", ylab="Change", size=12.5
  ); ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.7, .55, 0.25, 0.25) -> GGG_B.3
bth_reg2_diff %>% ggplot() +
  geom_sf(mapping = aes(fill = bi_class, geometry=geometry), color = "white", size = 0.1, show.legend = F) +
  bi_scale_fill(pal = "GrPink", dim = 2) + labs(title = "Both criteria") + bi_theme() +
  theme(plot.title = element_text(size=15)) -> map; legend <- bi_legend(
    pal="GrPink", dim=2, xlab="DISLI", ylab="Change", size=12.5
  ); ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.7, .55, 0.25, 0.25) -> GGG_B.4
ggarrange(GGG_B.3,GGG_B.2,GGG_B.1,GGG_B.4, nrow = 2, ncol = 2) %>% annotate_figure(
  top = text_grob("Change in prediabetes prevalence vs. DISLI by region (2020-2021)", color = "black", face = "bold", size = 14))

## By state ##
all_sta2_diff %>% ggplot() +
  geom_sf(mapping = aes(fill = bi_class, geometry=geometry), color = "white", size = 0.1, show.legend = t) +
  bi_scale_fill(pal = "GrPink", dim = 2) + labs(title = "Any criteria") + bi_theme() +
  theme(plot.title = element_text(size=13.5)) -> map; legend <- bi_legend(
    pal="GrPink", dim=2, xlab="DISLI", ylab="Change", size=12.5
  ); ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.7, .55, 0.25, 0.25) -> fig4A
a1c_sta2_diff %>% ggplot() +
  geom_sf(mapping = aes(fill = bi_class, geometry=geometry), color = "white", size = 0.1, show.legend = F) +
  bi_scale_fill(pal = "GrPink", dim = 2) + labs(title = "HbA1c criteria") + bi_theme() +
  theme(plot.title = element_text(size=13.5)) -> map; legend <- bi_legend(
    pal="GrPink", dim=2, xlab="DISLI", ylab="Change", size=12.5
  ); ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.7, .55, 0.25, 0.25) -> fig4B
ifg_sta2_diff %>% ggplot() +
  geom_sf(mapping = aes(fill = bi_class, geometry=geometry), color = "white", size = 0.1, show.legend = F) +
  bi_scale_fill(pal = "GrPink", dim = 2) + labs(title = "IFG criteria") + bi_theme() +
  theme(plot.title = element_text(size=13.5)) -> map; legend <- bi_legend(
    pal="GrPink", dim=2, xlab="DISLI", ylab="Change", size=12.5
  ); ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.7, .55, 0.25, 0.25) -> fig4C
bth_sta2_diff %>% ggplot() +
  geom_sf(mapping = aes(fill = bi_class, geometry=geometry), color = "white", size = 0.1, show.legend = F) +
  bi_scale_fill(pal = "GrPink", dim = 2) + labs(title = "Both criteria") + bi_theme() +
  theme(plot.title = element_text(size=13.5)) -> map; legend <- bi_legend(
    pal="GrPink", dim=2, xlab="DISLI", ylab="Change", size=12.5
  ); ggdraw() + draw_plot(map, 0, 0, 1, 1) + draw_plot(legend, 0.7, .55, 0.25, 0.25) -> fig4D

figure4 <- ggarrange(fig4C,fig4B,fig4A,fig4D, nrow = 2, ncol = 2) %>% annotate_figure(
  top = text_grob("\nChange in prediabetes prevalence vs. DISLI by state (2020-2021)", color = "black", face = "bold", size = 15))

ggsave(figure4, file="Figures/Figure4.pdf", bg="transparent",
       width=35, height=22, units=c("cm"), dpi=600, limitsize = FALSE)



#### Figure ADA ####
prev<-rbind(diabetes, prediabetes) %>% 
  mutate("cluster"=ordered(cluster, levels=c("Diabetes", "All", "None"),
                           labels=c("Diabetes",  "Prediabetes", "None"))) 
f1a<-ggplot(prev, aes(x=as.numeric(YEAR), y=prop,group=cluster,colour=cluster, linetype=cluster)) + 
  geom_line(size=1.5) + geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2,
                position=position_dodge(0.01), linetype=1)+ theme_pubclean()+
  ylab("Weighted prevalence (%)")+ xlab("ENSANUT cycle") + labs(colour="", linetype="") +
  scale_x_continuous(breaks = c(2016, 2018, 2020, 2021)) + ylim(8,32) +
  scale_color_manual(values=c("red4", "#552586")) + scale_linetype_manual(values=c(1,6)) +
  ggtitle("Disturbances in glucose metabolism") + theme(legend.position = "bottom") +
  theme(plot.title = element_text(size=15, face="bold", hjust=0.5, vjust=0))+
  geom_text(data=prev %>% filter(cluster=="Prediabetes"),
            aes(label = paste0(prop[1:4],"%","\n","(",prop[1:4]-IC95[1:4],"-",prop[1:4]+IC95[1:4],")")),
            size=3.0, vjust=-2.1)+
  geom_text(data=prev %>% filter(cluster=="Diabetes"),
            aes(label = paste0(prop,"%","\n","(",prop-IC95,"-",prop+IC95,")")),
            size=3.0, vjust=2.7)

f1b<-fig1B + geom_text(aes(label = paste0(prop,"%","\n","(",prop-IC95,"-",prop+IC95,")")),
                       size=3.0, col="white", position = position_stack(vjust = .5))

figabs1 <- ggarrange(f1a, f1b, ncol=1, nrow=2, labels=LETTERS[1:2])
figureabs<-ggarrange(figabs1, "", fig2B, ncol=3, nrow=1, widths = c(10,0.075,10), labels = c("","","C"))

ggsave(figureabs, file="Figures/Figure_ADA.pdf", bg="transparent",
       width=50, height=30, units=c("cm"), dpi=600, limitsize = FALSE)



