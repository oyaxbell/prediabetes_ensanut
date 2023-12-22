## ENSANUT database management ##
# Working directory
setwd("/Users/carlosfermin/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/Prediabetes ENSANUT")
#setwd("/Users/henry-macbook/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/research/Prediabetes ENSANUT")

#Load packages
pacman::p_load(readr, haven, tidyverse, data.table, readxl, survey)

#### DISLI------- ####
SLI <- readxl::read_xlsx("Bases/irs.xlsx") %>% mutate("entidad"=as.numeric(id)) %>% transmute(entidad, IRS, IRS_cat) %>% 
  mutate(dens=c(253.9,52.8,10.8,16.1,20.8,130.0,75.6,15.1, 6163.3,14.9,201.5,55.7,148.1,106.2,760.2,81.0, 404.1,44.4,90.2,
                44.1,191.9,202.6,41.6,46.2, 52.8,16.4,97.1,44.0,336.0,112.3,58.7,21.5), entidad2=as.factor(entidad)) %>%
  mutate("DISLI"=(lm(formula=IRS~dens) %>% residuals)) %>% mutate("DISLI_cat"=cut(DISLI, breaks = c(
    -Inf, quantile(x = DISLI, probs = c(.25,.5,.75)), Inf), labels=c("Very-Low", "Low", "Moderate", "High"))) %>%
  mutate("DISLI_cat2"=ifelse(DISLI_cat %in% c("Very-Low","Low", "Moderate"),0,1),
         "Region"=case_when(entidad%in% c(2,3,5,8,19,25,26,28)~1, entidad%in% c(1,6,10,11,14,16,18,24,32)~2,
                            entidad%in% c(9,13,15,17,21,22,29)~3, entidad%in% c(4,7,12,20,23,27,30,31)~4))

#### ENSA 2000--- ####
#Join data
salud_2000 <- read_csv2("Bases/ENSA_2000/Adultos.csv") %>%
  mutate("pond_sal"=p_adulto, "pond_ant"=p_adulto, "FOLIO_INT"=paste0(folio_v,integ))
social_2000 <- read_csv2("Bases/ENSA_2000/Hogar_hogar.csv") %>%
  mutate("est_final"=estrato)
social2_2000 <- read_csv2("Bases/ENSA_2000/Hogar_integrantes.csv") %>%
  mutate("pond_soc"=p_hogar, "FOLIO_INT" = paste0(folio_v, integ)) %>%
  subset(h204_1 == "A2") %>% rename(folio_vi=folio_v)

ens00 <- salud_2000 %>% merge(social2_2000, by="FOLIO_INT", all.x=T) %>%
  merge(social_2000, by="folio_v", all.x=T); n00 <- nrow(social2_2000)
rm(salud_2000, social_2000, social2_2000)

#Database cleaning
ens00_fin <- ens00 %>% mutate(
  # SOCIODEMOGRAPHIC
  "Sex" = case_when(sexo==1~1, sexo==2~0), # 0: Hombre, 1: Mujer
  "Age" = case_when(edada%in%1:110~edada), 
  "BMI" = a2901_2/(a2902_2/100)^2, 
  "Waist" = a2903_2,
  "WHtR" = as.numeric(a2903_2)/a2902_2,
  "Ob_central" = ifelse((Sex==1 & a2903_2 >= 80)|(Sex==0 & a2903_2 >= 90),1,0), # 0: No, 1: Sí
  "Year" = 2000,
  "Smoking" = case_when(a2101==0~0, a2103==2~1, a2103==1~2), # 0: Never, 1: Former, 2: Current
  "Area" = case_when(tam_loc.x==1~0, tam_loc.x==2~1), # 0: Rural, 1: Urban
  "State" = as.numeric(ent.x),
  "Region"=case_when(State%in% c(2,3,5,8,19,25,26,28)~1, State%in% c(1,6,10,11,14,16,18,24,32)~2, #1:North // 2:Central-West
                     State%in% c(9,13,15,17,21,22,29)~3, State%in% c(4,7,12,20,23,27,30,31)~4), #3:Central // 4:South-Southeast
  ##"Indigenous" = ifelse(), # 0: No, 1: Sí
  # LABORATORY (REVIEW)
  "Glucose" = a21003, # Capillary glucose
  "Fasting" = case_when(a21001==1~0, a21001==2~1), # 0: No Fasting, 1: Fasting
  # HEALTH HISTORY
  "HX_T2D" = case_when(a2502==1~1, a2502==2~0), # 0: No, 1: Sí
  "TX_T2D"= case_when(a2507==0~4, a2507==1~2, a2507==2~1, a2507==9~3), # 1:Insulin // 2:Pills // 3:Both // 4:None
  "HX_HBP" = case_when(a2603==1~1, a2603==2~0), # 0: No, 1: Sí
  "TX_HBP" = case_when(a2607==2~0, a2607==1~1), # 0: No, 1: Sí
  "HX_CKD" = case_when(a2702_3==1~1, a2702_3==2~0), # 0: No, 1: Sí
  "CKD_dialysis" = ifelse((a2703_1==3|a2703_2==3|a2703_3==3|a2703_4==3|a2703_5==3|a2703_6==3), 1,0) # 0: No, 1: Sí
  ## No AMI, HF or stroke variables
  )

#SBP/DBP
ens00_fin$a2904_2[ens00_fin$a2904_2 == 999] <- NA; ens00_fin$a2905_2[ens00_fin$a2905_2 == 999] <- NA
ens00_fin$a2904_3[ens00_fin$a2904_3 == 999] <- NA; ens00_fin$a2905_3[ens00_fin$a2905_3 == 999] <- NA
ens00_fin$SBP <- ens00_fin %>% select(a2904_2, a2905_2) %>% apply(1, mean, na.rm=T)
ens00_fin$DBP <- ens00_fin %>% select(a2904_3, a2905_3) %>% apply(1, mean, na.rm=T)

#Diabetes
ens00_fin$diabetes_previous <- with(ens00_fin, ifelse(HX_T2D==1|TX_T2D%in%1:3, 1, 0))
ens00_fin$diabetes_biochem <- with(ens00_fin, ifelse(Glucose>=126&Fasting==1, 1, 0))
ens00_fin$diabetes_undx <- with(ens00_fin, ifelse(diabetes_previous==0&diabetes_biochem==1, 1, 0))
ens00_fin$diabetes_fin <- with(ens00_fin, ifelse(diabetes_previous==1|diabetes_biochem==1, 1, 0))

# SLI and DISLI
ens00_fin$entidad <- ens00_fin$State; e16.SLI <- merge(ens00_fin, SLI, by="entidad", all.x = T)
ens00_fin$DISLI <- e16.SLI$DISLI; ens00_fin$DISLI_cat <- e16.SLI$DISLI_cat; ens00_fin$DISLI_cat2 <- e16.SLI$DISLI_cat2

setattr(ens00_fin$Sex, "description", "0:Men // 1:Women")
setattr(ens00_fin$Smoking, "description", "0:Never // 1:Former // 2:Current")
setattr(ens00_fin$Area, "description", "0:Rural // 1:Urban")
setattr(ens00_fin$TX_T2D, "description", "1:Insulin // 2:Pills // 3:Both // 4:Any")

#### ENSANUT 2006 ####
#Join data
social_2006 <- read_csv("Bases/ENSANUT_2006/Hogar_integrantes.csv") %>%
  mutate("pond_soc"=p_int, "est_final" = estrato, "FOLIO_INT"=paste0(folio, inth))
salud_2006 <- read_csv("Bases/ENSANUT_2006/Adultos.csv") %>%
  mutate("pond_sal"=p_adult, "pond_ant"=p_adult, "pond_ant"=p_adult,
         "FOLIO_INT"=paste0(folio, inth))
antro_2006 <- read_csv("Bases/ENSANUT_2006/antropometria_sr.csv") %>%
  mutate("FOLIO_INT"=paste0(FOLIO, inth))
sangre_2006 <- read_csv("Bases/ENSANUT_2006/BASE_CRONICAS _ENSANUT.csv") %>%
  mutate("FOLIO_INT"=paste0(folio, inth))

ens06 <- salud_2006 %>% merge(social_2006, by="FOLIO_INT", all.x=T) %>%
  merge(antro_2006, by="FOLIO_INT", all.x=T) %>%
  left_join(sangre_2006, by="FOLIO_INT"); n06 <- nrow(social_2006)

#Database cleaning
ens06_fin <- ens06 %>% mutate(
  # SOCIODEMOGRAPHIC
  "Sex" = case_when(sexo.x==1~0, sexo.x==2~1), # 0: Hombre, 1: Mujer
  "Age" = case_when(edad.x%in%1:110~edad.x), 
  "BMI" = PESON.x/(TALLAN.x/100)^2, 
  "Waist" = cintura.x,
  "WHtR" = as.numeric(cintura.x)/TALLAN.x,
  "Ob_central" = ifelse((Sex==1 & Waist >= 80)|(Sex==0 & Waist >= 90),1,0), # 0: No, 1: Sí
  "Year" = 2006,
  "Smoking" = case_when(a1301==3~0, a1302==2~1, a1302==1~2), # 0: Never, 1: Former, 2: Current
  "Area" = case_when(tam_loc.x==1~0, tam_loc.x==2~1, tam_loc.x==3~1), # 0: Rural, 1: Urban
  "State" = as.numeric(ent.x),
  "Region"=case_when(State%in% c(2,3,5,8,19,25,26,28)~1, State%in% c(1,6,10,11,14,16,18,24,32)~2, #1:North // 2:Central-West
                     State%in% c(9,13,15,17,21,22,29)~3, State%in% c(4,7,12,20,23,27,30,31)~4), #3:Central // 4:South-Southeast
  "Indigenous" = case_when(h211.x==2~0, h211.x==1~1), # 0: No, 1: Sí
  # LABORATORY
  "Glucose" = GLU,
  "HBA1C" = HBGlu,
  "Insulin" = Ins,
  "TG" = Trig,
  "CT" = Col,
  "HDL" = HDLC,
  "LDL" = LDL,
  "Fasting" = halim,
  # HEALTH HISTORY
  "HX_T2D" = case_when(a401==1~1, a401 %in% 2:3 ~ 0), # 0: No, 1: Sí
  "TX_T2D"= a406, # 1:Insulin // 2:Pills // 3:Both // 4:None
  "HX_HBP" = case_when(a501==1~1, a501==2~0), # 0: No, 1: Sí
  "TX_HBP" = case_when(a505==2~0, a505==1~1), # 0: No, 1: Sí
  "HX_CKD" = case_when(a701c==1~1, a701c==2~0), # 0: No, 1: Sí
  "CKD_dialysis" = ifelse((a703a==3|a703b==3|a703c==3|a703d==3|a703e==3|a703f==3|a703g==3|a703h==3|
                             a703a==5|a703b==5|a703c==5|a703d==5|a703e==5|a703f==5|a703g==5|a703h==5),
                          1,0), # 0: No, 1: Sí
  "HX_AMI" = case_when(a602a==1~1, a602a==2~0), 
  "HX_HFA" = case_when(a602c==1~1, a602c==2~0),
  "HX_EVC" = case_when(a708==1~1, a708 %in% c(2,9)~0)
)

#SBP/DBP
ens06_fin$sis1.x <- as.numeric(ens06_fin$sis1.x); ens06_fin$sis2.x <- as.numeric(ens06_fin$sis2.x)
ens06_fin$dias1.x <- as.numeric(ens06_fin$dias1.x); ens06_fin$dias2.x <- as.numeric(ens06_fin$dias2.x)
ens06_fin$sis1.x[ens06_fin$sis1.x==222]<-NA; ens06_fin$sis2.x[ens06_fin$sis2.x==222]<-NA
ens06_fin$dias1.x[ens06_fin$dias1.x==222]<-NA; ens06_fin$dias2.x[ens06_fin$dias2.x==222]<-NA
ens06_fin$SBP <- ens06_fin %>% select(sis1.x, sis2.x) %>% apply(1, mean, na.rm=T)
ens06_fin$DBP <- ens06_fin %>% select(dias1.x, dias2.x) %>% apply(1, mean, na.rm=T)

#CVD
ens06_fin$HX_CVD <- ((ens06_fin %>% select(HX_AMI, HX_HFA, HX_EVC) %>% apply(1, sum, na.rm=T))==0) %>%
  ifelse(0,1); ens06_fin$HX_CVD[with(ens06_fin, is.na(HX_AMI)&is.na(HX_HFA)&is.na(HX_EVC))] <- NA
#Diabetes
ens06_fin$diabetes_previous <- with(ens06_fin, ifelse(HX_T2D==1|TX_T2D%in%1:3, 1, 0))
ens06_fin$diabetes_biochem <- with(ens06_fin, ifelse(Glucose>=126&Fasting>=8, 1, 0))
ens06_fin$diabetes_undx <- with(ens06_fin, ifelse(diabetes_previous==0&diabetes_biochem==1, 1, 0))
ens06_fin$diabetes_fin <- with(ens06_fin, ifelse(diabetes_previous==1|diabetes_biochem==1, 1, 0))

## HOMA-IR
#ens06_fin <- ens06_fin %>% mutate("FOLIO_INT"=row_number()); homa_2016<-read.csv("Bases/ensanut.2016.HOMA.csv")
#homa_2016[homa_2016=="Unable to calculate"]<-NA; homa_2016$HOMA2IR[homa_2016$HOMA2IR==0]<-NA
#homa_2016$HOMA2IR<-as.numeric(homa_2016$HOMA2IR); homa_2016$HOMA2B<-as.numeric(homa_2016$HOMA2B)
#homa_2016$homa_cat<-ifelse(homa_2016$HOMA2IR>=2.5,1,0); ens06_fin<-ens06_fin%>% left_join(homa_2016, by="FOLIO_INT")

## SLI and DISLI
ens06_fin$entidad <- ens06_fin$State; e16.SLI <- merge(ens06_fin, SLI, by="entidad", all.x = T)
ens06_fin$DISLI <- e16.SLI$DISLI; ens06_fin$DISLI_cat <- e16.SLI$DISLI_cat; ens06_fin$DISLI_cat2 <- e16.SLI$DISLI_cat2

setattr(ens06_fin$Sex, "description", "0:Men // 1:Women")
setattr(ens06_fin$Smoking, "description", "0:Never // 1:Former // 2:Current")
setattr(ens06_fin$Area, "description", "0:Rural // 1:Urban")
setattr(ens06_fin$Indigenous, "description", "0:Non-indigenous // 1:Indigenous")
setattr(ens06_fin$TX_T2D, "description", "1:Insulin // 2:Pills // 3:Both // 4:Any")

#### ENSANUT 2012 ####
#Join data
social_2012 <- read_sav("Bases/ENSANUT_2012/Hogar_integrantes.sav") %>%
  mutate("pond_soc"=pondei, "est_final"=est_var, "FOLIO_INT"=paste0(folio, intp))
salud_2012 <- read_sav("Bases/ENSANUT_2012/Adultos.sav") %>%
  mutate("pond_sal"=pondef, "FOLIO_INT"=paste0(folio, intp))
antro_2012 <- read_sav("Bases/ENSANUT_2012/Antropometria.sav") %>%
  mutate("pond_ant"=pondef, "FOLIO_INT"=paste0(folio, intp))
sangre_2012 <- read_sav("Bases/ENSANUT_2012/Datos_sangre_capilar.sav") %>%
  mutate("pond_san"=pondef, "FOLIO_INT"=paste0(folio, intp))
hba1c_2012 <- read_sav("Bases/ENSANUT_2012/base hba1c_svy_diab.sav") %>%
  mutate("pond_a1c"=pondev3, "FOLIO_INT"=paste0(folio, intp))
glu_2012 <- read_sav("Bases/ENSANUT_2012/base glucosa.sav") %>%
  mutate("pond_glu"=PONDEV3, "FOLIO_INT"=paste0(folio, intp))

ens12 <- salud_2012 %>% merge(social_2012, by="FOLIO_INT", all.x=T) %>%
  merge(antro_2012, by="FOLIO_INT", all.x=T) %>%
  left_join(sangre_2012, by="FOLIO_INT") %>%
  left_join(hba1c_2012, by="FOLIO_INT") %>%
  left_join(glu_2012, by="FOLIO_INT"); n12 <- nrow(social_2012)
remove(social_2012,salud_2012,antro_2012,sangre_2012,hba1c_2012,glu_2012)

#Database cleaning
ens12_fin <- ens12 %>% mutate(
  #SOCIODEMOGRAPHIC
  "Sex"=case_when(sexo.x==1~0, sexo.x==2~1), #0:Men // 1:Women
  "Age"=case_when(edad.x%in%1:110~edad.x),
  "Height"=case_when(talla!=222.2~talla),
  "Weight"=case_when(peso!=222.2~peso),
  "Waist"=case_when(cintura!=222.2~cintura),
  "BMI"=Weight/(Height/100)^2,
  "WHtR"=Waist/Height,
  "Ob_central"=ifelse((Sex==1 & Waist>=80)|(Sex==0 & Waist>=90),1,0),
  "Year"=2012,
  "Smoking"=case_when(a1301==3~0, a1303a==8~1, a1303a%in%1:5~2), #0:Never // 1:Former // 2:Current
  "Smoking_intensity"=case_when(a1303a==8~0, a1303a%in%1:5~a1303a), #0:Non-smoker // 1:Daily // 2:Weekly // 3:Monthly // 4:Occasionally // 5:Yearly
  "Area"=case_when(est_urb.x==1~0, est_urb.x==2~1, est_urb.x==3~1), #0:Rural // 1:Urban
  "State"=entidad.x,
  "Region"=case_when(State%in% c(2,3,5,8,19,25,26,28)~1, State%in% c(1,6,10,11,14,16,18,24,32)~2, #1:North // 2:Central-West
                     State%in% c(9,13,15,17,21,22,29)~3, State%in% c(4,7,12,20,23,27,30,31)~4), #3:Central // 4:South-Southeast
  "Indigenous"=case_when(h212==2~0, h212==1~1), #0:Non-indigenous // 1:Indigenous,
  #LABORATORY DATA
  "Glucose"=glucosa,
  "HBA1C"=hba1c2,
  "Fasting"=sanvenh,
  #HEALTH HISTORY
  "HX_T2D"=case_when(a301==2~0, a303==1~0, a301==1~1), #0:No // 1:Yes
  "TX_T2D"=a307, #1:Insulin // 2:Pills // 3:Both // 4:None
  "HX_HBP"=case_when(a401==2~0, a401==1~1), #0:No // 1:Yes
  "TX_HBP"=case_when(a405==2~0, a405==1~1), #0:No // 1:Yes
  "HX_CKD"=NA,
  "CKD_dialysis"=NA,
  "HX_AMI"=case_when(a502a==2~0, a502a==1~1), #0:No // 1:Yes
  "HX_HFA"=case_when(a502c==2~0, a502c==1~1), #0:No // 1:Yes
  "HX_EVC"=case_when(a604==2~0, a604==1~1)) #0:No // 1:Yes

#SBP/DBP
ens12_fin$sistol1[ens12_fin$sistol1==999]<-NA; ens12_fin$sistol2[ens12_fin$sistol2==999]<-NA
ens12_fin$diastol1[ens12_fin$diastol1==999]<-NA; ens12_fin$diastol2[ens12_fin$diastol2==999]<-NA
ens12_fin$SBP <- ens12_fin %>% select(sistol1, sistol2) %>% apply(1, mean, na.rm=T)
ens12_fin$DBP <- ens12_fin %>% select(diastol1, diastol2) %>% apply(1, mean, na.rm=T)
#CVD
ens12_fin$HX_CVD <- ((ens12_fin %>% select(HX_AMI, HX_HFA, HX_EVC) %>% apply(1, sum, na.rm=T))==0) %>%
  ifelse(0,1); ens12_fin$HX_CVD[with(ens12_fin, is.na(HX_AMI)&is.na(HX_HFA)&is.na(HX_EVC))] <- NA
#Diabetes
ens12_fin$diabetes_previous <- with(ens12_fin, ifelse(HX_T2D==1|TX_T2D%in%1:3, 1, 0))
ens12_fin$diabetes_biochem <- with(ens12_fin, ifelse(Glucose>=126&Fasting>=8, 1, 0))
ens12_fin$diabetes_undx <- with(ens12_fin, ifelse(diabetes_previous==0&diabetes_biochem==1, 1, 0))
ens12_fin$diabetes_fin <- with(ens12_fin, ifelse(diabetes_previous==1|diabetes_biochem==1, 1, 0))
#SLI and DISLI
ens12_fin$entidad <- ens12_fin$State; e12.SLI <- merge(ens12_fin, SLI, by="entidad", all.x = T)
ens12_fin$DISLI <- e12.SLI$DISLI; ens12_fin$DISLI_cat <- e12.SLI$DISLI_cat; ens12_fin$DISLI_cat2 <- e12.SLI$DISLI_cat2

#Variable description
setattr(ens12_fin$Sex, "description", "0:Men // 1:Women")
setattr(ens12_fin$Smoking, "description", "0:Never // 1:Former // 2:Current")
list("labels","label","class") %>% lapply(setattr, x=ens12_fin$Smoking_intensity, value=NULL); setattr(
  ens12_fin$Smoking_intensity, "description", "#0:Non-smoker // 1:Daily // 2:Weekly // 3:Monthly // 4:Occasionally // 5:Yearly")
setattr(ens12_fin$Area, "description", "0:Rural // 1:Urban")
setattr(ens12_fin$Region, "description", "1:North // 2:Central-West // 3:Central // 4:South-Southeast")
setattr(ens12_fin$Indigenous, "description", "0:Non-indigenous // 1:Indigenous")
list("labels","label","class") %>% lapply(setattr, x=ens12_fin$TX_T2D, value=NULL); setattr(
  ens12_fin$TX_T2D, "description", "1:Insulin // 2:Pills // 3:Both // 4:Any")


#### ENSANUT 2016 ####
ens16 <- read_csv("Bases/ensanut2016.csv")
ens16_fin <- ens16 %>% mutate(
  #SOCIODEMOGRAPHIC
  "Sex"=case_when(sexo.x==1~0, sexo.x==2~1), #0:Men // 1:Women
  "Age"=case_when(edad.x%in%1:110~edad.x),
  "BMI"=peso/(talla/100)^2,
  "Waist"=cintura,
  "WHtR"=cintura/talla,
  "Ob_central"=ifelse((Sex==1 & Waist>=80)|(Sex==0 & Waist>=90),1,0),
  "Year"=2016,
  "Smoking"=case_when(a1301==3~0, a1301a==2~1, a1301a==1~2), #0:Never // 1:Former // 2:Current
  "Area"=case_when(rural==1~0, rural==2~1), #0:Rural // 1:Urban
  "State"=entidad.x,
  "Region"=case_when(State%in% c(2,3,5,8,19,25,26,28)~1, State%in% c(1,6,10,11,14,16,18,24,32)~2, #1:North // 2:Central-West
                     State%in% c(9,13,15,17,21,22,29)~3, State%in% c(4,7,12,20,23,27,30,31)~4), #3:Central // 4:South-Southeast
  "Indigenous"=case_when(h212==2~0, h212==1~1), #0:Non-indigenous // 1:Indigenous,
  #LABORATORY DATA
  "Glucose"=valor.GLU_SUERO,
  "HBA1C"=valor.HB1AC,
  "Insulin"=valor.INSULINA,
  "TG"=valor.TRIG,
  "CT"=valor.COLEST,
  "HDL"=valor.COL_HDL,
  "LDL"=valor.COL_LDL,
  "Albumin"=valor.ALBU,
  "Creatinine"=valor.CREAT,
  "Fasting"=sanvenh,
  #HEALTH HISTORY
  "HX_T2D"=case_when(a301.x%in%2:3~0, a301.x==1~1), #0:No // 1:Yes
  "TX_T2D"=a307, #1:Insulin // 2:Pills // 3:Both // 4:None
  "HX_HBP"=case_when(a401.x==2~0, a401.x==1~1), #0:No // 1:Yes
  "TX_HBP"=case_when(a405.y==2~0, a405.y==1~1), #0:No // 1:Yes
  "HX_CKD"=case_when(a605c==2~0, a605c==1~1), #0:No // 1:Yes
  "CKD_dialysis"=case_when((a606d==2|a605c==2)~0, a606d==1~1), #0:No // 1:Yes
  "HX_AMI"=case_when(a502a==2~0, a502a==1~1), #0:No // 1:Yes
  "HX_HFA"=case_when(a502c==2~0, a502c==1~1), #0:No // 1:Yes
  "HX_EVC"=case_when(a611==2~0, a611==1~1)) #0:No // 1:Yes

#SBP/DBP
ens16_fin$sistol3.x[ens16_fin$sistol3.x==999]<-NA; ens16_fin$sistol4.x[ens16_fin$sistol4.x==999]<-NA
ens16_fin$diastol3.x[ens16_fin$diastol3.x==999]<-NA; ens16_fin$diastol4.x[ens16_fin$diastol4.x==999]<-NA
ens16_fin$SBP <- ens16_fin %>% select(sistol3.y, sistol4.y) %>% apply(1, mean, na.rm=T)
ens16_fin$DBP <- ens16_fin %>% select(diastol3.y, diastol4.y) %>% apply(1, mean, na.rm=T)
#CVD
ens16_fin$HX_CVD <- ((ens16_fin %>% select(HX_AMI, HX_HFA, HX_EVC) %>% apply(1, sum, na.rm=T))==0) %>%
  ifelse(0,1); ens16_fin$HX_CVD[with(ens16_fin, is.na(HX_AMI)&is.na(HX_HFA)&is.na(HX_EVC))] <- NA
#Diabetes
ens16_fin$diabetes_previous <- with(ens16_fin, ifelse(HX_T2D==1|TX_T2D%in%1:3, 1, 0))
ens16_fin$diabetes_biochem <- with(ens16_fin, ifelse((Glucose>=126&Fasting>=8)|HBA1C>=6.5, 1, 0))
ens16_fin$diabetes_undx <- with(ens16_fin, ifelse(diabetes_previous==0&diabetes_biochem==1, 1, 0))
ens16_fin$diabetes_fin <- with(ens16_fin, ifelse(diabetes_previous==1|diabetes_biochem==1, 1, 0))
## HOMA-IR
ens16_fin <- ens16_fin %>% mutate("FOLIO_INT"=row_number()); homa_2016<-read.csv("Bases/ensanut.2016.HOMA.csv")
homa_2016[homa_2016=="Unable to calculate"]<-NA; homa_2016$HOMA2IR[homa_2016$HOMA2IR==0]<-NA
homa_2016$HOMA2IR<-as.numeric(homa_2016$HOMA2IR); homa_2016$HOMA2B<-as.numeric(homa_2016$HOMA2B)
homa_2016$homa_cat<-ifelse(homa_2016$HOMA2IR>=2.5,1,0); ens16_fin<-ens16_fin%>% left_join(homa_2016, by="FOLIO_INT")
## SLI and DISLI
ens16_fin$entidad <- ens16_fin$State; e16.SLI <- merge(ens16_fin, SLI, by="entidad", all.x = T)
ens16_fin$DISLI <- e16.SLI$DISLI; ens16_fin$DISLI_cat <- e16.SLI$DISLI_cat; ens16_fin$DISLI_cat2 <- e16.SLI$DISLI_cat2

setattr(ens16_fin$Sex, "description", "0:Men // 1:Women")
setattr(ens16_fin$Smoking, "description", "0:Never // 1:Former // 2:Current")
setattr(ens16_fin$Area, "description", "0:Rural // 1:Urban")
setattr(ens16_fin$Region, "description", "1:North // 2:Central-West // 3:Central // 4:South-Southeast")
setattr(ens16_fin$Indigenous, "description", "0:Non-indigenous // 1:Indigenous")
setattr(ens16_fin$TX_T2D, "description", "1:Insulin // 2:Pills // 3:Both // 4:Any")



#### ENSANUT 2018 ####
ens18 <- read_csv("Bases/ensanut2018.csv")
ens18.edu <- haven::read_sav("Bases/CS_RESIDENTES.sav")
ens18.edu$ID <- with(ens18.edu, paste0(UPM, VIV_SEL, HOGAR, NUMREN))
ens18$ID <- with(ens18, paste0(UPM.x, VIV_SEL.x, HOGAR.x, NUMREN.x))
ens18 <- merge(ens18, ens18.edu %>% select(ID, NIVEL), by="ID", all.x = T)

ens18_fin <- ens18 %>% mutate(
  #SOCIODEMOGRAPHIC
  "Sex"=case_when(SEXO.x==1~0, SEXO.x==2~1),
  "Age"=case_when(EDAD.x%in%1:110~EDAD.x),
  "Year"=2018,
  "Smoking"=case_when(P13_4==3~0, P13_2==3~1, P13_2%in%1:2~2), #0:Never // 1:Former // 2:Current
  "Area"=case_when(DOMINIO.x==1~1, DOMINIO.x==2~0), #0:Rural // 1:Urban
  "State"=ENT.x,
  "Region"=case_when(State%in% c(2,3,5,8,19,25,26,28)~1, State%in% c(1,6,10,11,14,16,18,24,32)~2, #1:North // 2:Central-West
                     State%in% c(9,13,15,17,21,22,29)~3, State%in% c(4,7,12,20,23,27,30,31)~4), #3:Central // 4:South-Southeast
  #LABORATORY DATA
  "Glucose"=VALOR_GLU_SUERO,
  "HBA1C"=VALOR_HB1AC,
  "Insulin"=VALOR_INSULINA,
  "TG"=VALOR_TRIG,
  "CT"=VALOR_COLEST,
  "HDL"=VALOR_COL_HDL,
  "LDL"=VALOR_COL_LDL,
  "Albumin"=VALOR_ALBUM,
  "Creatinine"=VALOR_CREAT,
  "Fasting"=P5_1.y,
  #HEALTH HISTORY
  "HX_T2D"=case_when(P3_1%in%2:3~0, P3_1==1~1), #0:No // 1:Yes
  "TX_T2D"=P3_8, #1:Insulin // 2:Pills // 3:Both // 4:None
  "HX_HBP"=case_when(P4_1.x==2~0, P4_1.x==1~1), #0:No // 1:Yes
  "TX_HBP"=case_when(P4_4==2~0, P4_4==1~1), #0:No // 1:Yes
  "HX_CKD"=case_when(P6_1_3.x==2~0, P6_1_3.x==1~1), #0:No // 1:Yes
  "CKD_dialysis"=case_when((P6_2_3==2&P6_2_4==2|P6_1_3.x==2)~0, (P6_2_3==1|P6_2_4==1)~1), #0:No // 1:Yes
  "HX_AMI"=case_when(P5_2_1==2~0, P5_2_1==1~1), #0:No // 1:Yes
  "HX_HFA"=case_when(P5_2_3==2~0, P5_2_3==1~1), #0:No // 1:Yes
  "HX_EVC"=case_when(P5_6==2~0, P5_6==1~1)) #0:No // 1:Yes

##Unified Height
ens18_fin$Height <- ens18_fin$TALLA4_1
ens18_fin$Height[ens18_fin$EDAD.x>=60&!is.na(ens18_fin$TALLA15_1)] <- with(ens18_fin, TALLA15_1[!is.na(TALLA15_1)])
##Unified Weight
ens18_fin$Weight <- ens18_fin$PESO1_1
ens18_fin$Weight[ens18_fin$EDAD.x>=60&!is.na(ens18_fin$PESO12_1)] <- with(ens18_fin, PESO12_1[!is.na(PESO12_1)])
##Unified Waist
ens18_fin$Waist <- ens18_fin$CIRCUNFERENCIA8_1
ens18_fin$Waist[ens18_fin$EDAD.x>=60&!is.na(ens18_fin$CINTURA21_1)] <- with(ens18_fin, CINTURA21_1[!is.na(CINTURA21_1)])
#BMI, WHtR, Central Obesity
ens18_fin <- ens18_fin %>% mutate("BMI"=Weight/(Height/100)^2, "WHtR"=Waist/Height,
                                  "Ob_central"=ifelse((Sex==1 & Waist>=80)|(Sex==0 & Waist>=90),1,0))

#SBP/DBP
ens18_fin$P27_1_1[ens18_fin$P27_1_1==222]<-NA; ens18_fin$P27_2_1[ens18_fin$P27_2_1==222]<-NA
ens18_fin$P27_1_2[ens18_fin$P27_1_2==222]<-NA; ens18_fin$P27_2_2[ens18_fin$P27_2_2==222]<-NA
ens18_fin$SBP <- ens18_fin %>% select(P27_1_1, P27_2_1) %>% apply(1, mean, na.rm=T)
ens18_fin$DBP <- ens18_fin %>% select(P27_1_2, P27_2_2) %>% apply(1, mean, na.rm=T)
#CVD
ens18_fin$HX_CVD <- ((ens18_fin %>% select(HX_AMI, HX_HFA, HX_EVC) %>% apply(1, sum, na.rm=T))==0) %>%
  ifelse(0,1); ens18_fin$HX_CVD[with(ens18_fin, is.na(HX_AMI)&is.na(HX_HFA)&is.na(HX_EVC))] <- NA
#Diabetes
ens18_fin$diabetes_previous <- with(ens18_fin, ifelse(HX_T2D==1|TX_T2D%in%1:3, 1, 0))
ens18_fin$diabetes_biochem <- with(ens18_fin, ifelse((Glucose>=126&Fasting>=8)|HBA1C>=6.5, 1, 0))
ens18_fin$diabetes_undx <- with(ens18_fin, ifelse(diabetes_previous==0&diabetes_biochem==1, 1, 0))
ens18_fin$diabetes_fin <- with(ens18_fin, ifelse(diabetes_previous==1|diabetes_biochem==1, 1, 0))
## HOMA-IR
ens18_fin$FOLIO_INT <- ens18_fin$id; homa_2018<-read.csv("Bases/ensanut.2018.HOMA.csv")
homa_2018[homa_2018=="Unable to calculate"]<-NA; homa_2018$HOMA2IR[homa_2018$HOMA2IR==0]<-NA
homa_2018$HOMA2IR<-as.numeric(homa_2018$HOMA2IR); homa_2018$HOMA2B<-as.numeric(homa_2018$HOMA2B)
homa_2018$homa_cat<-ifelse(homa_2018$HOMA2IR>=2.5,1,0); ens18_fin<-ens18_fin%>% left_join(homa_2018, by="FOLIO_INT")
## SLI and DISLI
ens18_fin$entidad <- ens18_fin$State; e18.SLI <- merge(ens18_fin, SLI, by="entidad", all.x = T)
ens18_fin$DISLI <- e18.SLI$DISLI; ens18_fin$DISLI_cat <- e18.SLI$DISLI_cat; ens18_fin$DISLI_cat2 <- e18.SLI$DISLI_cat2
#Indigenous language (0:Non-indigenous // 1:Indigenous)
lengua_2018<-read_csv("Bases/LENGUA_2018.csv") %>% mutate("ID_lengua"=paste0(UPM,"_",VIV_SEL,"_",HOGAR,"_",NUMREN)) %>%
  select(ID_lengua,P3_11); ens18_fin <- ens18_fin %>% mutate("ID_lengua"=paste0(UPM.x,"_",VIV_SEL.x,"_",HOGAR.x,"_",NUMREN.x))
ens18_fin<-merge(ens18_fin, lengua_2018, by="ID_lengua", all.x = T); ens18_fin$Indigenous <- with(ens18_fin, case_when(P3_11.y==2~0, P3_11.y==1~1))

setattr(ens18_fin$Sex, "description", "0:Men // 1:Women")
setattr(ens18_fin$Smoking, "description", "0:Never // 1:Former // 2:Current")
setattr(ens18_fin$Area, "description", "0:Rural // 1:Urban")
setattr(ens18_fin$Region, "description", "1:North // 2:Central-West // 3:Central // 4:South-Southeast")
setattr(ens18_fin$Indigenous, "description", "0:Non-indigenous // 1:Indigenous")
setattr(ens18_fin$TX_T2D, "description", "1:Insulin // 2:Pills // 3:Both // 4:Any")


#### ENSANUT 2020 ####
#Join data
social_2020 <- read_csv("Bases/ENSANUT_2020/integrantes_ensanut2020_w.csv") %>%
  mutate("pond_soc"=ponde_f20, "est_final"=est_sel)
salud_2020 <- read_csv("Bases/ENSANUT_2020/adultos_vac_tab_ensanut2020_w.csv") %>%
  mutate("pond_sal"=ponde_g20)
antro_2020 <- read_csv("Bases/ENSANUT_2020/antropometria_ensanut2020_w.csv") %>%
  mutate("pond_ant"=ponde_g20)
sangre_2020 <- read_csv("Bases/ENSANUT_2020/toma_de_sangre_ensanut2020_diab_w.csv") %>%
  mutate("pond_lab"=ponde_g20)
covid_2020 <- read_sav("Bases/ENSANUT_2020/toma_de_sangre_ensanut2020_covid_w.sav") %>%
  mutate("pond_cov"=ponde_g20)

ens20 <- social_2020 %>%
  merge(salud_2020, by="FOLIO_INT", all=T) %>%
  merge(antro_2020, by="FOLIO_INT", all=T) %>%
  left_join(sangre_2020, by="FOLIO_INT", na_matches="never") %>%
  left_join(covid_2020, by="FOLIO_INT", na_matches="never")

#Database cleaning
ens20_fin <- ens20 %>% mutate(
  #SOCIODEMOGRAPHIC
  "Sex"=case_when(H0302.x==1~0, H0302.x==2~1), # 0: Hombre, 1: Mujer
  "Age"=case_when(H0303.x%in%1:110~H0303.x),
  "BMI"=an01_1/(an04_01/100)^2,
  "Waist"=NA,
  "WHtR"=NA,
  "Year"=2020,
  "Smoking"=case_when(ADUL1A04A==3~0, ADUL1A01==3~1, ADUL1A01%in%1:2~2), #0:Never // 1:Former // 2:Current (NA in laboratory subset)
  "Area"=case_when(area_20.x==1~0, area_20.x==2~1), #0:Rural // 1:Urban
  "State"=ENTIDAD.x %>% as.numeric,
  "Region"=case_when(State%in% c(2,3,5,8,19,25,26,28)~1, State%in% c(1,6,10,11,14,16,18,24,32)~2, #1:North // 2:Central-West
                     State%in% c(9,13,15,17,21,22,29)~3, State%in% c(4,7,12,20,23,27,30,31)~4), #3:Central // 4:South-Southeast
  "Indigenous"=case_when(H0311==2~0, H0311==1~1), #0:Non-indigenous // 1:Indigenous,
  #LABORATORY DATA
  "Glucose"=valor.GLU_SUERO,
  "HBA1C"=HB1AC.Valor,
  "Insulin"=valor.INSULINA,
  "TG"=valor.TRIG,
  "CT"=valor.COLEST,
  "HDL"=valor.COL_HDL,
  "LDL"=valor.COL_LDL,
  "Albumin"=valor.ALBUM,
  "Creatinine"=valor.CREAT,
  "Uric_acid"=valor.AC_URICO,
  "ALT"=valor.ALT,
  "AST"=valor.AAT,
  "GGT"=valor.GGT,
  "Fasting"=san04,
  #HEALTH HISTORY 
  "HX_T2D"=ifelse((H0902A.x=="1" | H0902B=="1" | H0902C=="1" | H0902D=="1")==T,1,0), #0:No // 1:Yes
  "TX_T2D"=NA,
  "HX_HBP"=ifelse((H0902A.x=="3" | H0902B=="3" | H0902C=="3" | H0902D=="3")==T,1,0), #0:No // 1:Yes
  "TX_HBP"=NA,
  "HX_CKD"=NA,
  "HX_CVD"=ifelse((H0902A.x=="4" | H0902B=="4" | H0902C=="4" | H0902D=="4")==T,1,0)) #0:No // 1:Yes

#SBP/DBP
ens20_fin$an08_01s[ens20_fin$an08_01s==999]<-NA; ens20_fin$an08_01d[ens20_fin$an08_01d==999]<-NA
ens20_fin$an08_02s[ens20_fin$an08_02s==999]<-NA; ens20_fin$an08_02d[ens20_fin$an08_02d==999]<-NA
ens20_fin$an08_03s[ens20_fin$an08_03s==999]<-NA; ens20_fin$an08_03d[ens20_fin$an08_03d==999]<-NA
ens20_fin$SBP <- ens20_fin %>% select(an08_01s, an08_02s, an08_03s) %>% apply(1, mean, na.rm=T)
ens20_fin$DBP <- ens20_fin %>% select(an08_01d, an08_02d, an08_03d) %>% apply(1, mean, na.rm=T)
#Diabetes
ens20_fin$diabetes_previous <- with(ens20_fin, ifelse(HX_T2D==1, 1, 0))
ens20_fin$diabetes_biochem <- with(ens20_fin, ifelse((Glucose>=126&Fasting>=8)|HBA1C>=6.5, 1, 0))
ens20_fin$diabetes_undx <- with(ens20_fin, ifelse(diabetes_previous==0&diabetes_biochem==1, 1, 0))
ens20_fin$diabetes_fin <- with(ens20_fin, ifelse(diabetes_previous==1|diabetes_biochem==1, 1, 0))
## HOMA-IR
#ens20_fin$FOLIO_INT <- ens20_fin$FOLIO_INT; homa_2020<-read.csv("Bases/ensanut.2020.HOMA.csv")
#homa_2020[homa_2020=="Unable to calculate"]<-NA; homa_2020$HOMA2IR[homa_2020$HOMA2IR==0]<-NA
#homa_2020$HOMA2IR<-as.numeric(homa_2020$HOMA2IR); homa_2020$HOMA2B<-as.numeric(homa_2020$HOMA2B)
#homa_2020$homa_cat<-ifelse(homa_2020$HOMA2IR>=2.5,1,0); ens20_fin<-ens20_fin%>% left_join(homa_2020, by="FOLIO_INT")
## SLI and DISLI
ens20_fin$entidad <- ens20_fin$State; e20.SLI <- merge(ens20_fin, SLI, by="entidad", all.x = T)
ens20_fin$DISLI <- e20.SLI$DISLI; ens20_fin$DISLI_cat <- e20.SLI$DISLI_cat; ens20_fin$DISLI_cat2 <- e20.SLI$DISLI_cat2

setattr(ens20_fin$Sex, "description", "0:Men // 1:Women")
setattr(ens20_fin$Smoking, "description", "0:Never // 1:Former // 2:Current")
setattr(ens20_fin$Area, "description", "0:Rural // 1:Urban")
setattr(ens20_fin$Region, "description", "1:North // 2:Central-West // 3:Central // 4:South-Southeast")
setattr(ens20_fin$Indigenous, "description", "0:Non-indigenous // 1:Indigenous")
setattr(ens20_fin$TX_T2D, "description", "1:Insulin // 2:Pills // 3:Both // 4:Any")


#### ENSANUT 2021 ####
#Join data
social_2021 <- read_csv(
  "Bases/ENSANUT_2021/integrantes_ensanut2021_w_12_01_2022.csv") %>%
  mutate("pond_soc"=ponde_f, "est_final"=est_sel)
salud_2021 <- read_csv(
  "Bases/ENSANUT_2021/ensadul2021_entrega_w_15_12_2021.csv") %>%
  mutate("pond_sal"=ponde_f)
antro_2021 <- read_csv(
  "Bases/ENSANUT_2021/ensaantro21_entrega_w_17_12_2021.csv") %>%
  mutate("pond_ant"=ponde_f)
sangre_2021 <- read_spss(
  "Bases/ENSANUT_2021/ensasangre21_entrega_w_integrada_ok.sav") %>%
  mutate("pond_lab"=ponde_vv)

ens21 <- social_2021 %>%
  merge(salud_2021, by="FOLIO_INT", all.y=T) %>%
  merge(antro_2021, by="FOLIO_INT", all=T) %>%
  left_join(sangre_2021, by="FOLIO_INT")
remove(social_2021,salud_2021,antro_2021,sangre_2021)

#Data cleaning
ens21_fin <- ens21 %>% mutate(
  #SOCIODEMOGRAPHIC
  "Sex"=case_when(sexo==1~0, sexo==2~1),
  "Age"=case_when(edad%in%1:110~edad),
  "Year"=2021,
  "Smoking"=case_when(a1305==3~0, a1301==3~1, a1301%in%1:2~2), #0:Never // 1:Former // 2:Current
  "Area"=case_when(estrato.x==1~0, estrato.x==2~1, estrato.x==3~1), #0:Rural // 1:Urban
  "State"=entidad.x,
  "Region"=case_when(State%in% c(2,3,5,8,19,25,26,28)~1, State%in% c(1,6,10,11,14,16,18,24,32)~2, #1:North // 2:Central-West
                     State%in% c(9,13,15,17,21,22,29)~3, State%in% c(4,7,12,20,23,27,30,31)~4), #3:Central // 4:South-Southeast
  "Indigenous"=case_when(h0311==2~0, h0311==1~1), #0:Non-indigenous // 1:Indigenous
  #LABORATORY DATA
  "Glucose"=valor_GLU_SUERO,
  "HBA1C"=valor_HB1AC,
  "Insulin"=valor_INSULINA,
  "TG"=valor_TRIG,
  "CT"=valor_COLEST,
  "HDL"=valor_COL_HDL,
  "LDL"=valor_COL_LDL,
  "Albumin"=valor_ALBUM,
  "Creatinine"=valor_CREAT,
  "Uric_acid"=valor_AC_URICO,
  "CRP.mg_L"=valor_PROTCREAC,
  "Fasting"=san04,
  #HEALTH HISTORY (0:No // 1:Yes)
  "HX_T2D"=case_when(a0301.x%in%2:3~0, a0301.x==1~1), #0:No // 1:Yes
  "TX_T2D"=a0307, #1:Insulin // 2:Pills // 3:Both // 4:None
  "HX_HBP"=case_when(a0401%in%2:3~0, a0401==1~1), #0:No // 1:Yes
  "TX_HBP"=case_when(a0404==2~0, a0404==1~1), #0:No // 1:Yes
  "HX_CKD"=case_when(a0601c==2~0, a0601c==1~1), #0:No // 1:Yes
  "CKD_dialysis"=case_when((a0602c==2&a0602d==2|a0601c==2)~0, (a0602c==1|a0602d==1)~1), #0:No // 1:Yes
  "HX_AMI"=case_when((a0502a==2|a0501==2)~0, a0502a==1~1), #0:No // 1:Yes
  "HX_HFA"=case_when((a0502c==2|a0501==2)~0, a0502c==1~1), #0:No // 1:Yes
  "HX_EVC"=case_when((a0506==2|a0501==2)~0, a0506==1~1)) #0:No // 1:Yes

## HOMA-IR
ens21_fin$FOLIO_INT <- ens21_fin$FOLIO_INT; homa_2021<-read.csv("Bases/ensanut.2021.HOMA.csv")
homa_2021[homa_2021=="Unable to calculate"]<-NA; homa_2021$HOMA2IR[homa_2021$HOMA2IR==0]<-NA
homa_2021$HOMA2IR<-as.numeric(homa_2021$HOMA2IR); homa_2021$HOMA2B<-as.numeric(homa_2021$HOMA2B)
homa_2021$homa_cat<-ifelse(homa_2021$HOMA2IR>=2.5,1,0); ens21_fin<-ens21_fin%>% left_join(homa_2021, by="FOLIO_INT")
##Unified Height
ens21_fin <- ens21_fin %>% filter(edad>=20)
ens21_fin$Height <- ens21_fin$an04_1
ens21_fin$Height[ens21_fin$Age>=60&!is.na(ens21_fin$an15_1)] <- with(ens21_fin, an15_1[!is.na(an15_1)])
##Unified Weight
ens21_fin$Weight <- ens21_fin$an01_1
ens21_fin$Weight[ens21_fin$Age>=60&!is.na(ens21_fin$an12_1)] <- with(ens21_fin, an12_1[!is.na(an12_1)])
##Unified Waist
ens21_fin$Waist <- ens21_fin$an08_1
ens21_fin$Waist[ens21_fin$Age>=60&!is.na(ens21_fin$an21_1)] <- with(ens21_fin, an21_1[!is.na(an21_1)])
#BMI, WHtR, Central Obesity
ens21_fin <- ens21_fin %>% mutate("BMI"=Weight/(Height/100)^2, "WHtR"=Waist/Height,
                                  "Ob_central"=ifelse((Sex==1 & Waist>=80)|(Sex==0 & Waist>=90),1,0))

#SBP/DBP
ens21_fin$an27_01s[ens21_fin$an27_01s==999]<-NA; ens21_fin$an27_01d[ens21_fin$an27_01d==999]<-NA
ens21_fin$an27_02s[ens21_fin$an27_02s==999]<-NA; ens21_fin$an27_02d[ens21_fin$an27_02d==999]<-NA
ens21_fin$an27_03s[ens21_fin$an27_03s==999]<-NA; ens21_fin$an27_03d[ens21_fin$an27_03d==999]<-NA
ens21_fin$SBP <- ens21_fin %>% select(an27_01s, an27_02s, an27_03s) %>% apply(1, mean, na.rm=T)
ens21_fin$DBP <- ens21_fin %>% select(an27_01d, an27_02d, an27_03d) %>% apply(1, mean, na.rm=T)
#CVD
ens21_fin$HX_CVD <- ((ens21_fin %>% select(HX_AMI, HX_HFA, HX_EVC) %>% apply(1, sum, na.rm=T))==0) %>%
  ifelse(0,1); ens21_fin$HX_CVD[with(ens21_fin, is.na(HX_AMI)&is.na(HX_HFA)&is.na(HX_EVC))] <- NA
#Diabetes
ens21_fin$diabetes_previous <- with(ens21_fin, ifelse(HX_T2D==1|TX_T2D%in%1:3, 1, 0))
ens21_fin$diabetes_biochem <- with(ens21_fin, ifelse((Glucose>=126&Fasting>=8)|HBA1C>=6.5, 1, 0))
ens21_fin$diabetes_undx <- with(ens21_fin, ifelse(diabetes_previous==0&diabetes_biochem==1, 1, 0))
ens21_fin$diabetes_fin <- with(ens21_fin, ifelse(diabetes_previous==1|diabetes_biochem==1, 1, 0))
## SLI and DISLI
ens21_fin$entidad <- ens21_fin$State; e21.SLI <- merge(ens21_fin, SLI, by="entidad", all.x = T)
ens21_fin$DISLI <- e21.SLI$DISLI; ens21_fin$DISLI_cat <- e21.SLI$DISLI_cat; ens21_fin$DISLI_cat2 <- e21.SLI$DISLI_cat2

setattr(ens21_fin$Sex, "description", "0:Men // 1:Women")
setattr(ens21_fin$Smoking, "description", "0:Never // 1:Former // 2:Current")
setattr(ens21_fin$Area, "description", "0:Rural // 1:Urban")
setattr(ens21_fin$Region, "description", "1:North // 2:Central-West // 3:Central // 4:South-Southeast")
setattr(ens21_fin$Indigenous, "description", "0:Non-indigenous // 1:Indigenous")
setattr(ens21_fin$TX_T2D, "description", "1:Insulin // 2:Pills // 3:Both // 4:Any")

#### ENSANUT 2022 ####
ens22 <- read_csv("Bases/ensanut2022.csv")

ens22_fin <- ens22 %>% mutate(
  #SOCIODEMOGRAPHIC
  "Sex"=case_when(sexo==1~0, sexo==2~1), #0:Men // 1:Women
  "Age"=case_when(edad%in%1:110~edad),
  "Year"=2022,
  "Smoking"=case_when(a1305==3~0, a1301==3~1, a1301%in%1:2~2), #0:Never // 1:Former // 2:Current
  "Area"=case_when(estrato.x==1~0, estrato.x==2~1, estrato.x==3~1), #0:Rural // 1:Urban
  "State"=entidad.x,
  "Region"=case_when(State%in% c(2,3,5,8,19,25,26,28)~1, State%in% c(1,6,10,11,14,16,18,24,32)~2, #1:North // 2:Central-West
                     State%in% c(9,13,15,17,21,22,29)~3, State%in% c(4,7,12,20,23,27,30,31)~4), #3:Central // 4:South-Southeast
  "Indigenous"=case_when(h0311==2~0, h0311==1~1), #0:Non-indigenous // 1:Indigenous
  #LABORATORY DATA
  "Glucose"=valor_GLU_SUERO,
  "HBA1C"=valor_HB1AC,
  "Insulin"=valor_INSULINA,
  "TG"=valor_TRIG,
  "CT"=valor_COLEST,
  "HDL"=valor_COL_HDL,
  "LDL"=valor_COL_LDL,
  "Albumin"=valor_ALBU,
  "Creatinine"=valor_CREAT,
  "Uric_acid"=valor_AC_URICO,
  "CRP.mg_L"=valor_PCR, #mg/L 
  "Fasting"=san04,
  #HEALTH HISTORY (0:No // 1:Yes)
  "HX_T2D"=case_when(a0301%in%2:3~0, a0301==1~1), #0:No // 1:Yes
  "TX_T2D"=a0307, #1:Insulin // 2:Pills // 3:Both // 4:None
  "HX_HBP"=case_when(a0401%in%2:3~0, a0401==1~1), #0:No // 1:Yes
  "TX_HBP"=case_when(a0404==2~0, a0404==1~1), #0:No // 1:Yes
  "HX_CKD"=case_when(a0601c==2~0, a0601c==1~1), #0:No // 1:Yes
  "CKD_dialysis"=NA,
  "HX_AMI"=case_when((a0502a==2)~0, a0502a==1~1), #0:No // 1:Yes
  "HX_HFA"=case_when((a0502c==2)~0, a0502c==1~1), #0:No // 1:Yes
  "HX_EVC"=case_when((a0502d==2)~0, a0502d==1~1)) #0:No // 1:Yes

## HOMA-IR
ens22_fin$FOLIO_INT <- ens22_fin$FOLIO_INT; homa_2022<-read.csv("Bases/ensanut.2022.HOMA.csv")
homa_2022[homa_2022=="Unable to calculate"]<-NA; homa_2022$HOMA2IR[homa_2022$HOMA2IR==0]<-NA
homa_2022$HOMA2IR<-as.numeric(homa_2022$HOMA2IR); homa_2022$HOMA2B<-as.numeric(homa_2022$HOMA2B)
homa_2022$homa_cat<-ifelse(homa_2022$HOMA2IR>=2.5,1,0); ens22_fin<-ens22_fin%>% left_join(homa_2022, by="FOLIO_INT")
##Unified Height
ens22_fin <- ens22_fin %>% filter(edad>=20)
ens22_fin$Height <- ens22_fin$an04_1 
ens22_fin$Height[ens22_fin$Age>=60&!is.na(ens22_fin$an15_1 )] <- with(
  ens22_fin, an15_1 [ens22_fin$Age>=60&!is.na(ens22_fin$an15_1 )])
##Unified Weight
ens22_fin$Weight <- ens22_fin$an01_1 
ens22_fin$Weight[ens22_fin$Age>=60&!is.na(ens22_fin$an12_1 )] <- with(
  ens22_fin, an12_1 [ens22_fin$Age>=60&!is.na(ens22_fin$an12_1 )])
##Unified Waist
ens22_fin$waist.A <- ens22_fin %>% select(an08_1 , an08_2 ) %>% apply(1, mean, na.rm=T)
ens22_fin$waist.B <- ens22_fin %>% select(an21_1 , an21_2 ) %>% apply(1, mean, na.rm=T)
ens22_fin$Waist <- ens22_fin$waist.A
ens22_fin$Waist[ens22_fin$Age>=60&!is.na(ens22_fin$waist.B)] <- with(
  ens22_fin, waist.B[ens22_fin$Age>=60&!is.na(ens22_fin$waist.B)])

#BMI, WHtR, Central Obesity
ens22_fin <- ens22_fin %>% mutate("BMI"=Weight/(Height/100)^2, "WHtR"=Waist/Height,
                                  "Ob_central"=ifelse((Sex==1 & Waist>=80)|(Sex==0 & Waist>=90),1,0))
#SBP/DBP
ens22_fin$an27_01s[ens22_fin$an27_01s==999]<-NA; ens22_fin$an27_01d[ens22_fin$an27_01d==999]<-NA
ens22_fin$an27_02s[ens22_fin$an27_02s==999]<-NA; ens22_fin$an27_02d[ens22_fin$an27_02d==999]<-NA
ens22_fin$an27_03s[ens22_fin$an27_03s==999]<-NA; ens22_fin$an27_03d[ens22_fin$an27_03d==999]<-NA
ens22_fin$SBP <- ens22_fin %>% select(an27_01s, an27_02s, an27_03s) %>% apply(1, mean, na.rm=T)
ens22_fin$DBP <- ens22_fin %>% select(an27_01d, an27_02d, an27_03d) %>% apply(1, mean, na.rm=T)
#CVD
ens22_fin$HX_CVD <- ((ens22_fin %>% select(HX_AMI, HX_HFA, HX_EVC) %>% apply(1, sum, na.rm=T))==0) %>%
  ifelse(0,1); ens22_fin$HX_CVD[with(ens22_fin, is.na(HX_AMI)&is.na(HX_HFA)&is.na(HX_EVC))] <- NA
#Diabetes
ens22_fin$diabetes_previous <- with(ens22_fin, ifelse(HX_T2D==1|TX_T2D%in%1:3, 1, 0))
ens22_fin$diabetes_biochem <- with(ens22_fin, ifelse((Glucose>=126&Fasting>=8)|HBA1C>=6.5, 1, 0))
ens22_fin$diabetes_undx <- with(ens22_fin, ifelse(diabetes_previous==0&diabetes_biochem==1, 1, 0))
ens22_fin$diabetes_fin <- with(ens22_fin, ifelse(diabetes_previous==1|diabetes_biochem==1, 1, 0))
## SLI and DISLI
ens22_fin$entidad <- ens22_fin$State; e21.SLI <- merge(ens22_fin, SLI, by="entidad", all.x = T)
ens22_fin$DISLI <- e21.SLI$DISLI; ens22_fin$DISLI_cat <- e21.SLI$DISLI_cat; ens22_fin$DISLI_cat2 <- e21.SLI$DISLI_cat2

setattr(ens22_fin$Sex, "description", "0:Men // 1:Women")
setattr(ens22_fin$Smoking, "description", "0:Never // 1:Former // 2:Current")
setattr(ens22_fin$Area, "description", "0:Rural // 1:Urban")
setattr(ens22_fin$Region, "description", "1:North // 2:Central-West // 3:Central // 4:South-Southeast")
setattr(ens22_fin$Indigenous, "description", "0:Non-indigenous // 1:Indigenous")
setattr(ens22_fin$TX_T2D, "description", "1:Insulin // 2:Pills // 3:Both // 4:Any")




#### Save as .csv ####
write_csv(ens00_fin, "Bases/ENSANUT_final/ensanut2000_fin.csv")
write_csv(ens06_fin, "Bases/ENSANUT_final/ensanut2006_fin.csv")
write_csv(ens12_fin, "Bases/ENSANUT_final/ensanut2012_fin.csv")
write_csv(ens16_fin, "Bases/ENSANUT_final/ensanut2016_fin.csv")
write_csv(ens18_fin, "Bases/ENSANUT_final/ensanut2018_fin.csv")
write_csv(ens20_fin, "Bases/ENSANUT_final/ensanut2020_fin.csv")
write_csv(ens21_fin, "Bases/ENSANUT_final/ensanut2021_fin.csv")
write_csv(ens22_fin, "Bases/ENSANUT_final/ensanut2022_fin.csv")

#### Variable names ####
#2016
ens16$sexo.x #Sex
ens16$edad.x #Age
ens16$imc #BMI
ens16$cintura #Waist
ens16$sistol3.y #SBP
ens16$diastol3.y #DBP
ens16$a1301a; ens16$a1301 #Smoking
ens16$rural #Area
ens16$entidad.x #Region
ens16$h212 #Ind_lan
ens16$valor.GLU_SUERO #Glucose
ens16$valor.HB1AC #HbA1c
ens16$valor.INSULINA #Insulin
ens16$valor.TRIG #Triglycerides
ens16$valor.COLEST #Cholesterol
ens16$valor.COL_HDL #HDLc
ens16$valor.COL_LDL #LDLc
ens16$valor.ALBU #Albumin
ens16$valor.CREAT #Creatinine
ens16$a301.x==1 #Prev_T2D
ens16$a307 %in% c(1,2,3) #T2D_medications
ens16$a401.x==1 #Prev_HBP
ens16$a405.y==1 #HBP_medications
ens16$a605c #Prev_CKD
ens16 %>% select(a502a, a502c, a611) #Prev_CDV


#2018
ens18$SEXO.x #Sex
ens18$EDAD.x #Age
ens18 %>% select(TALLA4_1, TALLA15_1) #Height
ens18 %>% select(PESO1_1, PESO12_1) #Weight
ens18 %>% select(CIRCUNFERENCIA8_1, CINTURA21_1) #Waist
ens18 %>% select(P27_1_1, P27_2_1) #SBP
ens18 %>% select(P27_1_2, P27_2_2) #DBP
ens18 %>% select(P13_2, P13_4) #Smoking
ens18$DOMINIO.x #Area
ens18$ENT.x #Region
ens18$P3_11 #Ind_lan
ens18$VALOR_GLU_SUERO #Glucose
ens18$VALOR_HB1AC #HbA1c
ens18$VALOR_INSULINA #Insulin
ens18$VALOR_TRIG #Triglycerides
ens18$VALOR_COLEST #Cholesterol
ens18$VALOR_COL_HDL #HDLc
ens18$VALOR_COL_LDL #LDLc
ens18$VALOR_ALBUM #Albumin
ens18$VALOR_CREAT #Creatinine
ens18$P3_1==1 #Prev_T2D:  3=No
ens18$P3_8 %in% c(1,2,3) #T2D_medications: 4=No
ens18$P4_1.x==1 #Prev_HBP
ens18$P4_4==1 #HBP_medications
ens18$P6_1_3.x #Prev_CKD: 1=Yes // 2=No
ens18 %>% select(P5_2_1, P5_2_3, P5_6) #Prev_CDV: 1=Yes // 2=No


#2020
ens20$Sexo.y #Sex
ens20$H0303.x #Age
ens20$an01_1 #Weight
ens20$an04_01 #Height
NA #Waist
ens20 %>% select(an08_01s, an08_02s, an08_03s) #SBP
ens20 %>% select(an08_01d, an08_02d, an08_03d) #DBP
ens20 %>% select(ADUL1A01,ADUL1A04A) #Smoking (NA in laboratory subset)
ens20$area_20.x #Area
ens20$ENTIDAD.x #Region
ens20$H0311 #Ind_lan
ens20$valor.GLU_SUERO #Glucose
ens20$HB1AC.Valor #HbA1c
ens20$valor.INSULINA #Insulin
ens20$valor.TRIG #Triglycerides
ens20$valor.COLEST #Cholesterol
ens20$valor.COL_HDL #HDLc
ens20$valor.COL_LDL #LDLc
ens20$valor.ALBUM #Albumin
ens20$valor.CREAT #Creatinine
ens20$valor.AC_URICO
ens20$valor.ALT #ALT
ens20$valor.AAT #AST??
ens20$valor.GGT #GGT
ens20 %>% select(H0902A.x,H0902B,H0902C,H0902D)
#==1: Prev_T2D
#==3: Prev_HBP
#==4: Prev_CDV
NA #T2D_medications
NA #HBP_medications
NA #Prev_CKD


#2021
ens21$sexo #Sex
ens21$edad #Age
ens21 %>% select(an04_1, an15_1)#Height
ens21 %>% select(an01_1, an12_1) #Weight
ens21 %>% select(an08_1, an21_1) #Waist
ens21 %>% select(an27_01s, an27_02s, an27_03s) #SBP
ens21 %>% select(an27_01d, an27_02d, an27_03d) #DBP
ens21 %>% select(a1301, a1305) #Smoking
ens21$estrato.x #Area
ens21$entidad.x #Region
ens21$h0311 #Ind_lan
ens21$valor_GLU_SUERO #Glucose
ens21$valor_HB1AC #HbA1c
ens21$valor_INSULINA #Insulin
ens21$valor_TRIG #Triglycerides
ens21$valor_COLEST #Cholesterol
ens21$valor_COL_HDL #HDLc
ens21$valor_COL_LDL #LDLc
ens21$valor_ALBUM #Albumin
ens21$valor_CREAT #Creatinine
ens21$valor_AC_URICO #Uric acid
ens21$valor_PROTCREAC #CRP
ens21$a0301.x==1 #Prev_T2D:  %in%2:3 = No
ens21$a0307 %in% c(1,2,3) #T2D_medications: 4=No
ens21$a0401==1 #Prev_HBP %in%2:3 = No
ens21$a0404==1 #HBP_medications 2=No
ens21$a0601c==1 #Prev_CKD: 1=Yes // 2=No
with(ens21,(a0602c==1|a0602d==1)) #CKD_dialysis
ens21 %>% select(a0502a, a0502c, a0506) #Prev_CDV: 1=Yes // (VAR==2|a0501==2)


#2022