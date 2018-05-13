#ICD for identifying CVD subgroups in FRS 
#Coronary Heart disease -> IHD ( I20-I25 )
#Stroke -> Cerebrovascular disease (I60-I69)
#Congestive heart failure -> Heart failure (I50)
#Intermitten Claudiation -> Atherosclerotic Arteries disease (I70)

#Code of concern (according to Lagerweij, 2017): 
# MI: I21, I22 ; 
# OCHD: I20, I23, I24, I25; 
# Cardiac arrest: I46, R96; 
# Hemorrhagic stroke: I60, I61, I62; 
# Ischemic stroke: I63, I65; 
# Other Stroke: I64, I66; 
# Other CVD: G45, I67, I69, I70-I74, I50; 
# To add: TIA (G45), R96

pacman::p_load(dplyr, lubridate, tidyr, stringr,Hmisc) #icd)
#The icd package concerns the comorbidity score
setwd("~/Desktop/Validation FHS/Code")

#Hos data -------------------------
hos_r <- spss.get("../Data/hos_CVD.sav")
hos_raw <- hos_r %>% 
  mutate(date_h = ymd(INDATUMA)) %>% 
  subset(select=c(LopNr, date_h, HDIA3)) %>%
  arrange(LopNr) %>%
  distinct(LopNr, .keep_all = TRUE) #remove duplicated entries

#First event ; Grouped based on ICD 9 and 10
#NOI = Not Outcome of Interest
#Note: Pulmonary (I26-28) are excluded from group17
#Missing codes: G45, 798.1, 798.2, 798.9 

hos_icd <- hos_raw %>% 
  #filter(year_h < 1997) %>%  #ICD 9 was used till 1996 
  mutate(icd = str_trim(as.character(droplevels.factor(HDIA3)), "right")) %>% 
  subset(select = c(LopNr, date_h, icd)) %>%
  mutate(group_h = 
    #if_else(str_detect(icd, "39"), "Chr_Rheu_HD",
    if_else(str_detect(icd, "401|402|403|404|405"), "HyperTn",
    if_else(str_detect(icd, paste("4", 10:14, sep = "", collapse = "|")), "IHD",
    if_else(str_detect(icd, paste("4", 15:17, sep = "", collapse = "|")), "Pul_HD",
    if_else(str_detect(icd, paste("4", 20:29, sep = "", collapse = "|")), "Other_HD",
    if_else(str_detect(icd, paste("4", 30:38, sep = "", collapse = "|")), "CerebroVD",
    if_else(str_detect(icd, paste("4", 40:49, sep = "", collapse = "|")), "Arteries",
    if_else(str_detect(icd, paste("4", 51:59, sep = "", collapse = "|")), "Veins",
    if_else(str_detect(icd, paste("I", 10:15, sep = "", collapse = "|")), "HyperTn",
    if_else(str_detect(icd, paste("I", 20:25, sep = "", collapse = "|")), "IHD",
    if_else(str_detect(icd, paste("I", 26:28, sep = "", collapse = "|")), "Pul_HD",
    if_else(str_detect(icd, paste("I", 30:52, sep = "", collapse = "|")), "Other_HD",
    if_else(str_detect(icd, paste("I", 60:69, sep = "", collapse = "|")), "CerebroVD",
    if_else(str_detect(icd, paste("I", 70:79, sep = "", collapse = "|")), "Arteries",
    if_else(str_detect(icd, paste("I", 80:89, sep = "", collapse = "|")), "Veins",
    if_else(str_detect(icd, paste("I", 95:99, sep = "", collapse = "|")), "Unspecified",
    "NOI"))))))))))))))),
    group_h17 = 
    if_else(str_detect(icd, "I21|I22|410"), "MI",
    if_else(str_detect(icd, "I20|I23|I24|I25|411|412|413|414"), "OCHD",
    if_else(str_detect(icd, "I46|427|428"), "Cardiac arrest",
    if_else(str_detect(icd, "I60|I61|I62|430|431|432"), "Hemorrhagic Stroke",
    if_else(str_detect(icd, "I63|I65|433|434|435"), "Ischemic Stroke", #435=G45
    if_else(str_detect(icd, "I64|I66"), "Other Stroke",
    if_else(str_detect(icd, "I67|I69|I70|I71|I72|I73|I74|I50|436|437|438|440|441|442|443|444"), 
            "Other CVD", "NOI")))))))) #I50 = 798.1, 798.2, 798.9

table(hos_icd$group_h17)
#table(hos_icd$group_h17,hos_icd$icd)
#table(hos_icd$group_h,hos_icd$icd)

hos_icd1 <- hos_icd %>% filter(group_h17 != "NOI") %>% 
  mutate(type = "hos")

#filter(hos_icd1, LopNr == 205003)

#Death registy --------------
dod_raw <- spss.get("../Data/dod_CVD.sav", datevars = "dodsdatn")
dod_raw$date_d <- ymd(dod_raw$dodsdatn)

dod_icd <- dod_raw %>% 
  mutate(icd = str_trim(as.character(droplevels.factor(ULORS3)), "right")) %>% 
  subset(select = c(LopNr, date_d, icd)) %>%
  mutate(group_d = 
    #if_else(str_detect(icd, "39"), "Chr_Rheu_HD",
    if_else(str_detect(icd, "401|402|403|404|405"), "HyperTn",
    if_else(str_detect(icd, paste("4", 10:14, sep = "", collapse = "|")), "IHD",
    if_else(str_detect(icd, paste("4", 15:17, sep = "", collapse = "|")), "Pul_HD",
    if_else(str_detect(icd, paste("4", 20:29, sep = "", collapse = "|")), "Other_HD",
    if_else(str_detect(icd, paste("4", 30:38, sep = "", collapse = "|")), "CerebroVD",
    if_else(str_detect(icd, paste("4", 40:49, sep = "", collapse = "|")), "Arteries",
    if_else(str_detect(icd, paste("4", 51:59, sep = "", collapse = "|")), "Veins",
    if_else(str_detect(icd, paste("I", 10:15, sep = "", collapse = "|")), "HyperTn",
    if_else(str_detect(icd, paste("I", 20:25, sep = "", collapse = "|")), "IHD",
    if_else(str_detect(icd, paste("I", 26:28, sep = "", collapse = "|")), "Pul_HD",
    if_else(str_detect(icd, paste("I", 30:52, sep = "", collapse = "|")), "Other_HD",
    if_else(str_detect(icd, paste("I", 60:69, sep = "", collapse = "|")), "CerebroVD",
    if_else(str_detect(icd, paste("I", 70:79, sep = "", collapse = "|")), "Arteries",
    if_else(str_detect(icd, paste("I", 80:89, sep = "", collapse = "|")), "Veins",
    if_else(str_detect(icd, paste("I", 95:99, sep = "", collapse = "|")), "Unspecified",
    "NOI"))))))))))))))),
  group_d17 = 
    if_else(str_detect(icd, "I21|I22|410"), "MI",
    if_else(str_detect(icd, "I20|I23|I24|I25|411|412|413|414"), "OCHD",
    if_else(str_detect(icd, "I46|427|428"), "Cardiac arrest",
    if_else(str_detect(icd, "I60|I61|I62|430|431|432"), "Hemorrhagic Stroke",
    if_else(str_detect(icd, "I63|I65|433|434|435"), "Ischemic Stroke", #435=G45
    if_else(str_detect(icd, "I64|I66"), "Other Stroke",
    if_else(str_detect(icd, "I67|I69|I70|I71|I72|I73|I74|I50|436|437|438|440|441|442|443|444"), 
    "Other CVD", "NOI"))))))), #I50 = 798.1, 798.2, 798.9
  group_score = 
    if_else(str_detect(icd, "401|402|403|404|405"), "CVD",
    if_else(str_detect(icd, paste("4", 10:14, sep = "", collapse = "|")), "CVD",
    if_else(str_detect(icd, paste("4", 26:28, sep = "", collapse = "|")), "CVD",
    if_else(str_detect(icd, paste("4", 30:31, sep = "", collapse = "|")), "CVD",
    if_else(str_detect(icd, paste("4", 33:36, sep = "", collapse = "|")), "CVD",
    if_else(str_detect(icd, paste("4", 38:43, sep = "", collapse = "|")), "CVD",
    if_else(str_detect(icd, paste("I", 10:15, sep = "", collapse = "|")), "CVD",
    if_else(str_detect(icd, paste("I", 20:25, sep = "", collapse = "|")), "CVD",
    if_else(str_detect(icd, paste("I", 44:51, sep = "", collapse = "|")), "CVD",
    if_else(str_detect(icd, paste("I", 61:69, sep = "", collapse = "|")), "CVD",
    if_else(str_detect(icd, paste("I", 70:73, sep = "", collapse = "|")), "CVD",
            "NOI"))))))))))))

table(dod_icd$group_d17, dod_icd$group_score) 
#some cases are considered by one, but not another one

dod_icd1 <- dod_icd %>% 
  mutate(type = "dod" )

#First event data (for export) ---------------------
dh <- full_join(hos_icd1, dod_icd1, by = "LopNr") %>% 
  arrange(LopNr) %>% 
  mutate(
  #group_ev = if_else(is.na(group_h), group_d, group_h),
  group_ev17 = if_else(is.na(group_h17), group_d17, group_h17),
  date_ev = if_else(is.na(date_h), date_d, date_h),
  type_ev = if_else(is.na(group_h), "dod", "hos")) %>% 
  mutate(type_ev = factor(type_ev, levels = c("hos", "dod"))) %>% 
  subset(select=c(LopNr, date_ev, date_d, type_ev, group_ev17,group_score)) %>% 
  distinct(LopNr,date_ev, .keep_all = TRUE)

head(dh)
#table(dh$type_ev, dh$group_ev17)
#table(dh$type_ev, dh$group_ev)
write.csv(dh, file = "dh.csv" ,row.names = FALSE)

##Repeated measure---------------------------
hos_j <- hos_icd %>% rename(group = group_h,date = date_h)
dod_j <- dod_icd %>% rename(group = group_d,date = date_d)

dh1 <- bind_rows(hos_j, dod_j) %>% 
  arrange(LopNr) %>% 
  mutate(type_ev = factor(type_ev, levels = c("hos", "dod"))) %>%
  distinct(LopNr, date, type_ev, .keep_all = TRUE)
head(dh1)
str(dh1)
year(dh1$date)
round(prop.table(table(dh1$type_ev, dh1$group),2),3)
dh1[1:100,]

table(dh$group_ev, dh$type_ev)

#dh1 <- dh %>% 
#  subset(select=c(LopNr, date_ev, group_ev, type_ev))
#write.csv(dh1, file = "dh.csv" )