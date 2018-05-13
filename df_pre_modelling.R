# load packages and function ----------------------
pacman::p_load(dplyr, lubridate, tidyr, tidyverse, labelled, 
               Hmisc, RColorBrewer, gridExtra, gtable, grid,mice)

setwd("~/Desktop/Validation FHS/Code")
vip_r <- spss.get("../Data/VIP_raw.sav", datevars = "provdat")

#One repeated case in this vip dataset
vip_raw <- vip_r %>% 
  distinct(LopNr, provdat, .keep_all = TRUE) %>% 
  mutate(provdat = ymd(provdat),
         year_ex = year(provdat)) %>%
  arrange(LopNr, provdat) 

#Use first blods0 to identify diabetic obs. 
vip_raw$dm <- if_else(vip_raw$blods0 > 7.0, 1, 0)
#Some of them are not captured by blods0, use questionnaire
vip_raw$dm <- if_else(is.na(vip_raw$dm), ifelse(vip_raw$c6 == 1, 1, 0),  vip_raw$dm)
#Take result from blods2
vip_raw$dm <- if_else(is.na(vip_raw$dm), ifelse(vip_raw$blods2 > 12.2, 1, 0),  vip_raw$dm)

#Statin treatment
vip_raw$Tx_statin <- ifelse(vip_raw$c10.5==1|vip_raw$c5e==1, 1, 0)
vip_raw$Tx_statin <- ifelse(is.na(vip_raw$Tx_statin), 0, vip_raw$Tx_statin)
#Anti-HT Tx
vip_raw$Tx_ht <- ifelse(vip_raw$c10.1==1|vip_raw$c5a==1, 1, 0)
vip_raw$Tx_ht <- ifelse(is.na(vip_raw$Tx_ht), 0, vip_raw$Tx_ht)
#Angina Tx
vip_raw$Tx_angina <- ifelse(vip_raw$c10.2==1|vip_raw$c5b==1, 1, 0)
vip_raw$Tx_angina <- ifelse(is.na(vip_raw$Tx_angina), 0, vip_raw$Tx_angina)

#The changes of measurement in skol needs to be adjusted
#method as in Ng et al., 2012
skol_dat <- ymd(20090901)
vip_raw$skol <- ifelse(vip_raw$provdat < skol_dat, (vip_raw$skol- 0.170)/0.939, vip_raw$skol)

#SRH
vip_raw$c1 <- ifelse(vip_raw$c1==9,NA, vip_raw$c1)
table(vip_raw$c1) #impute NA into median value, 2 
vip_raw$c1 <- ifelse(is.na(vip_raw$c1), 2, vip_raw$c1) 

#Family History of CVD
vip_raw$famhis <- if_else(vip_raw$c2 == 1, 1, 0)
vip_raw$famhis <- if_else(is.na(vip_raw$famhis), 0, vip_raw$famhis)
#Yes= 1; No/Don't know/Missing = 0

table(is.na(vip_raw$h8b), vip_raw$year_ex)

smoke_v <- vip_raw %>% 
  reshape2::melt(id.vars = c("LopNr", "provdat"), 
                 measure.vars = c("h1a", "h1b", "h1c", "h1d", "h1e", "h1f", "h1g"), 
                 na.rm = TRUE) %>% 
  mutate(smoke = ifelse(
    variable == "h1a" | variable == "h1f" | variable == "h1g", 
    0, 1)) %>% 
  subset(select = c(LopNr, provdat, smoke)) %>%
  distinct(LopNr, provdat, .keep_all = TRUE) %>% arrange(LopNr) 

vip_j  <- vip_raw %>% subset(
  select=c(LopNr, provdat, year_ex,gender, alder, 
           sbt, skol, hdl, bmi, 
           SjvReg, Tx_angina, Tx_ht, Tx_statin, dm, c8, c1, famhis)) %>%
  left_join(smoke_v) %>% 
  mutate(smoke = ifelse(is.na(smoke), 0, smoke), #set all NA in smoke to non-smoker 
         dm = ifelse(is.na(dm), 0, dm)) #set all NA in dm to non-diabetic 

censor.date <- ymd(20141231)
dh <- read_csv("../Data/dh.csv")
dh <- dh %>%
  mutate(date_ev = ymd(date_ev),
         date_d = ymd(date_d))
str(dh)

#use only the first examination 
vip_tid <- vip_j %>% 
  left_join(dh, by = "LopNr") %>%
  filter(Tx_angina==0 & c8 ==2 & year_ex < 2015) %>%
  distinct(LopNr, .keep_all = TRUE)  %>% #take only the first participation 
  mutate(
    #dup = if_else(duplicated(LopNr) | duplicated(LopNr, fromLast = TRUE), 1,0),
    type_ev = ifelse(is.na(type_ev), "Censored", type_ev),
    grp_d = ifelse(is.na(group_score), "NOI", group_score),
    grp_dh = ifelse(is.na(group_ev17), "NOI", group_ev17),
    status_d = ifelse(grp_d=="NOI", 0, 1),
    status_dh = ifelse(grp_dh=="NOI", 0, 1),
    end_d = if_else(grp_d=="NOI", censor.date, date_d),
    end_dh = if_else(grp_dh=="NOI", censor.date, date_ev),
    time_d = time_length(interval(provdat, end_d), "year"),
    time_dh = time_length(interval(provdat, end_dh), "year")) %>% 
  mutate(
    Tx_angina = NULL,
    c8 = NULL,
    date_ev = NULL,
    date_d = NULL,
    SjvReg = factor(SjvReg, labels = c("Ume\u{00E5} region", "Skellefte\u{00E5} region", "S\u{00F6}dra Lappland")),
    gender = factor(gender, labels = c("Male", "Female"))) %>% 
  filter(alder>35, time_dh > 0) %>% 
  mutate(time_dh = time_dh + 0.1,
         time_d = time_d + 0.1)

#write.csv(vip_tid, file = "vip_first_visit_dh.csv", row.names = FALSE)


#Event in each time period for repeated observations -------- 

#number of repeated participation
table(vip_j %>% group_by(LopNr) %>% group_size())

vip_id_dat <- vip_j %>% 
  filter(Tx_angina==0 & c8 ==2 & year_ex < 2015 & alder>35) %>%
  subset(select=c(LopNr, provdat)) %>% 
  group_by(LopNr) %>% 
  mutate(time_visit = seq(n())) %>%  #create the times of visits 
  ungroup() %>% 
  filter(time_visit<4) %>% #only one observation repeated for 4 times 
  spread(time_visit,provdat) %>% 
  left_join(dh) %>% 
  mutate(date_d = if_else(is.na(date_d) | group_score =="NOI", censor.date, date_d),
         date_dh = if_else(is.na(date_ev) | group_ev17 =="NOI", censor.date, date_ev),
         date_dh = if_else(date_dh > date_d, date_d, date_dh),
         date_ev = NULL,
         group_score = ifelse(group_score == "NOI", NA, group_score),
         group_ev17 = ifelse(group_ev17 =="NOI", NA, group_ev17)) %>% 
  filter(date_dh > `1`) %>% 
  mutate(dh1 = if_else(is.na(`2`), time_length(interval(`1`, date_dh),"year"), 
                       if_else(date_dh < `2`, time_length(interval(`1`, date_dh),"year"),
                               time_length(interval(`1`, `2`),"year"))),
         dh2 = ifelse(is.na(`2`) | date_dh < `2`, NA, 
                       if_else(is.na(`3`), time_length(interval(`2`, date_dh),"year"),
                               if_else(date_dh < `3`, time_length(interval(`2`, date_dh),"year"),
                                       time_length(interval(`2`, `3`),"year")))),
         dh3 = ifelse(is.na(`3`) | date_dh < `3`, NA, time_length(interval(`3`, date_dh),"year")),
         status_dh1 = if_else(is.na(group_ev17), 0, 
                              if_else(is.na(`2`), 1, 
                                      if_else(date_dh < `2`, 1, 0))),
         status_dh2 = if_else(is.na(group_ev17), 0, 
                              if_else(is.na(`2`), 0,
                              if_else(is.na(`3`), 1, 
                                      if_else(date_dh < `3`, 1, 0)))),
         status_dh3 = if_else(is.na(group_ev17), 0,
                              if_else(is.na(`3`), 0, 1)),
         d1 = ifelse(is.na(`2`), time_length(interval(`1`, date_d),"year"), 
                     time_length(interval(`1`, `2`),"year")),
         d2 = ifelse(is.na(`2`), NA, 
                     ifelse(is.na(`3`),time_length(interval(`2`, date_d),"year"), 
                            time_length(interval(`2`, `3`),"year"))),
         d3 = ifelse(is.na(`3`), NA, time_length(interval(`3`, date_d), "year")),
         status_d1 = ifelse(is.na(group_score), 0, 
                            ifelse(is.na(`2`), 1, 0)),
         status_d2 = ifelse(is.na(group_score), 0,
                            ifelse(is.na(`2`), 0, 
                                   ifelse(is.na(`3`), 1, 0))),
         status_d3 = ifelse(is.na(group_score), 0, 
                            ifelse(is.na(`3`), 0,1)))

vip_id_dat1 <- vip_id_dat %>% 
  gather(time_visit, provdat, `1`:`3`) %>%
  drop_na(provdat) %>% 
  mutate(time_dh = if_else(time_visit == 3, dh3, 
                           if_else(time_visit ==2, dh2, dh1)),
         time_d = if_else(time_visit == 3, d3, 
                          if_else(time_visit == 2, d2, d1)),
         status_dh = if_else(time_visit == 3, status_dh3, 
                             if_else(time_visit ==2, status_dh2, status_dh1)),
         status_d = if_else(time_visit == 3, status_d3, 
                             if_else(time_visit ==2, status_d2, status_d1))) %>%
  drop_na(time_dh) %>% 
  subset(select=c(LopNr, provdat, group_ev17, group_score, status_d, time_d, status_dh, time_dh, type_ev))

#Check how the data looks like 
table(vip_id_dat1$time_dh > vip_id_dat1$time_d)
table(vip_j$SjvReg)

#Join data
vip_tid1 <- vip_id_dat1 %>% 
  left_join(vip_j) %>% 
  subset(select=-c(c8,Tx_angina)) %>% 
  mutate(SjvReg = factor(SjvReg, labels = c("S\u{00F6}dra Lappland","Skellefte\u{00E5} region", "Ume\u{00E5} region")),
         gender = factor(gender, labels = c("Male", "Female")))

table(vip_tid1$status_dh, vip_tid1$time_dh>10)
table(vip_tid1$status_d, vip_tid1$time_d>10)
prop.table(table(is.na(vip_tid1$sbt)))*100 #missing: 0.65%
prop.table(table(is.na(vip_tid1$hdl)))*100 #missing: 57.0%
prop.table(table(is.na(vip_tid1$skol)))*100 #missing: 0.53%
prop.table(table(is.na(vip_tid1$bmi)))*100 #missing: 0.49% 
write.csv(vip_tid1, file = "vip_dh.csv", row.names = FALSE)

#imputation for hdl-----------
hdl_df <- vip_tid1 %>% 
  subset(select=c(status_d, status_dh, time_d, time_dh, alder, gender, skol, hdl, sbt, Tx_ht, Tx_statin, dm, smoke,bmi, SjvReg)) %>% 
  drop_na(hdl) #remove missings in HDL

pred <- quickpred(hdl_df, minpuc = 0.25, include = c("gender","status_d","status_dh")) ; pred
pred[, "time_d"] <- 0
pred[, "time_dh"] <- 0
pred[, "SjvReg"] <- 0 ; pred

hdl_df_imp <- hdl_df %>% mice(m=5, maxit = 10, method = "pmm", seed = 1234, print = FALSE, pred = pred)
saveRDS(hdl_df_imp, file = "hdl_imp.rds")
