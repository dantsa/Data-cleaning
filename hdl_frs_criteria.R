# load packages and function ----------------------
#use the imputed dataset to perform modelling
#imputation method: nearest neighbour 
pacman::p_load(dplyr, survival, tidyr, tidyverse, labelled, Hmisc, 
               RColorBrewer,gridExtra,broom, pROC, pec, riskRegression,mice, PredictABEL,boot)

setwd("~/Desktop/Validation FHS/Code")
source(file = "gnd_cal.R")
source(file = "silly_func.R")

#read data ----------------- 

hdl_imp <- readRDS("hdl_imp.rds")

hdl_imp1 <- complete(hdl_imp,2)
hdl_imp1 <- hdl_imp1 %>% 
  mutate(Tx_ht_cat = factor(Tx_ht, labels = c("no Tx", "Tx")),
       dm_cat = factor(dm, labels = c("non-dm", "dm")),
       smoke_cat = factor(smoke, labels = c("non-s", "smoker")),
       alder = as.numeric(alder),
       skol_mgdl = skol/0.0259,
       hdl_mgdl = hdl/0.0259) 

#Preparation for calculating SCORE risk function ------
#Use ALL available sample... the idea is to extract the average mortality and risk factors in the local population 
##1.Expected 10-year cummulative CVD mortality rates
#km <- survfit(Surv(time_d, status_d)~gender+as.factor(alder), data=hdl_imp1)
km <- survfit(Surv(alder, status_d)~gender, data=hdl_imp1)

#km_1 <- survfit(Surv(time_dh, status_dh)~gender, data=hdl_imp1)
summary(km) #extract survival probability around 10 years 

#summary(km, times=min(10))$surv

km_sur_prob <- data.frame(gender = c(rep("Female",3), rep("Male",3)),
                          alder = c(rep(c(40,50,60),2)),
                          sur.prob = c(.998, .988, .936, .994, .969, .867), #incidence endpoint
                          #sur.prob = c(1, .999, .990, .999, .996, .978), #death endpoint -> skip it, totally wrong calibration 
                          hdl_coef = c(rep(log(0.63),3), rep(log(0.78),3)),
                          tc_coef = c(rep(log(1.26),3), rep(log(1.33),3)),
                          sbt_coef = c(rep(log(1.28),3), rep(log(1.37),3)),
                          smk_coef = c(rep(log(1.91),3), rep(log(1.80),3)) )

## 2. Expected SBP, TC, Smoking
#simplify it by using KM estimates and sample mean (proportion for smokers)
#it is the same as using regressionssss (codes of them at the end)
score_fact <- hdl_imp1 %>% 
  group_by(gender,alder) %>% 
  summarise(smk_prev = mean(as.integer(smoke)) - 1,
            mean_tc_mmol = mean(skol),
            mean_sbt = mean(sbt),
            mean_hdl = mean(hdl),
            mean_sum = 0.0185592*mean_sbt+0.004595*mean_tc_mmol+0.7213618*smk_prev) %>%
  ungroup()
score_fact <- as.data.frame(score_fact)

score_df <- km_sur_prob %>% left_join(score_fact) %>% 
  mutate(mean_sum_hdl = hdl_coef*mean_hdl + tc_coef*mean_tc_mmol + sbt_coef*mean_sbt + smk_coef*smk_prev) %>% 
  subset(select=c(gender, alder, sur.prob, mean_sum, mean_sum_hdl))

mean_df <- hdl_imp1 %>% 
  group_by(gender) %>% 
  summarise(smk_prev = mean(as.integer(smoke)) - 1,
            dm_prev = mean(as.integer(dm)) - 1,
            Tx_ht_prev = mean(as.integer(Tx_ht)) - 1,
            mean_sbt = mean(sbt),
            mean_tc = mean(skol_mgdl),
            mean_bmi = mean(bmi), 
            mean_hdl = mean(hdl_mgdl),
            mean_alder = mean(alder)) %>%
  ungroup() 
mean_df <- as.data.frame(mean_df)

str(hdl_imp1$gender)

#Use only the first region for training model ---------
hdl_models <- hdl_imp1 %>% 
  filter(SjvReg == "S\u00f6dra Lappland") %>% 
  group_by_(.dots="gender") %>%
  do(frs = coxph(Surv(time_dh, status_dh)~log(alder)+log(skol_mgdl)+log(hdl_mgdl)+log(sbt)+Tx_ht_cat+dm_cat+smoke_cat, 
                 data = ., model = TRUE),
     score_dh = coxph(Surv(alder, status_dh)~skol+sbt+smoke_cat, 
                   data = filter(.,dm==0), model = TRUE),
     frs_office = coxph(Surv(time_dh, status_dh)~log(alder)+log(bmi)+log(sbt)+Tx_ht_cat+dm_cat+smoke_cat, 
                 data = ., model = TRUE),
     score_dh_hdl = coxph(Surv(alder, status_dh)~skol+sbt+smoke_cat+hdl, 
                       data = filter(.,dm==0), model = TRUE),
     score_d = coxph(Surv(alder, status_d)~skol+sbt+smoke_cat, 
                     data = filter(.,dm==0), model = TRUE),
     score_d_hdl = coxph(Surv(alder, status_d)~skol+sbt+smoke_cat+hdl, 
           data = filter(.,dm==0), model = TRUE)) %>% 
  ungroup()

hdl_models <- as.data.frame(hdl_models) #avoid the error commandssss
summary(hdl_models$frs_office[[2]])
summary(hdl_models$score_dh_hdl[[1]])

basehaz(hdl_models$frs_office[[2]])
#sum_tab_exp(hdl_models$frs_office[[1]], file = "frs_f.txt")
#sum_tab_exp(hdl_models$frs_office[[2]], file = "frs_m.txt")

#Compute absolute risk - female ---  
ts_f <- hdl_imp1 %>% 
  filter(gender =="Female" & SjvReg != "Ume<U+00E5> region" ) %>% 
  mutate(status_dh_10 = ifelse(time_dh>10, 0, status_dh)) %>% 
  left_join(score_df) %>% 
  left_join(mean_df) %>% 
  mutate(
    pred.score.dh = 1- sur.prob^predict(hdl_models$score_dh[[1]], newdata = ., type = "risk"),
    pred.score.hdl.dh = 1- sur.prob^predict(hdl_models$score_dh_hdl[[1]], newdata = ., type = "risk"),
    pred.frs = risk.bh.predict(coxph.object = hdl_models$frs[[1]], newdata = ., time.p = 10),
    pred.frs.of = risk.bh.predict(coxph.object = hdl_models$frs_office[[1]], newdata = ., time.p = 10),
    a = 2.32888*log(alder) + 1.20904*log(skol_mgdl) - 0.70833*log(hdl_mgdl) + 2.76157*Tx_ht + 2.82263*log(sbt) + 0.52873*smoke + 0.69154*dm, 
    b = 2.32888*3.8686 + 1.20904*5.3504 - 0.70833*4.0176 + 2.76157*4.2400 + 2.82263*0.5826 + 0.52873*0.3423 + 0.69154*0.0376,  
    frs_ori = 1 - 0.95012^exp(a - b),
    a1 =  2.72107*log(alder) + 0.51125*log(bmi) + 2.88267*Tx_ht + 2.81291*log(sbt) + 0.61868*smoke + 0.77763*dm, 
    b1 =  3.11296*3.8560 + 0.79277*log(26.49642) + 1.92672*4.3544 + 1.85508*0.5019 + 0.70953*0.3522 + 0.53160*0.0650, 
    frsbmi_ori = 1 - 0.94833^exp(a1 - b1),
    a_s = 0.0185592*sbt+0.004595*skol+0.7213618*smoke,
    score_local = 1 - sur.prob^exp(a_s-mean_sum),
    a_s_hdl = log(0.63)*hdl+log(1.26)*skol+log(1.28)*sbt+log(1.91)*smoke,
    score_local_hdl = 1 - sur.prob^exp(a_s_hdl-mean_sum_hdl))

#show differences in distribution ------ 
plot_df1 <- data.frame(score = ts_f$pred.score.dh,
                       model = "SCORE",
                       outcome = ts_f$status_dh)
plot_df2 <- data.frame(score = ts_f$pred.score.hdl.dh,
                       model = "SCORE HDL",
                       outcome = ts_f$status_dh)
plot_df3 <- data.frame(score = ts_f$pred.frs,
                       model = "FRS",
                       outcome = ts_f$status_dh)
table(plot_df$outcome)
plot_df <- rbind(plot_df1,plot_df2,plot_df3)
plot_df$outcome <- factor(plot_df$outcome, labels = c("Survived", "Event"))
plot_df$mod_out <- interaction(plot_df$model, plot_df$outcome)
ggplot(data = filter(plot_df,model!="FRS"), aes(score, color = mod_out))+geom_density()

# compute AUC - female ------ 
#Rees-Score 
#dh-dh, dh-d, d-d, d-dh
roc.score.f <- roc(status_dh_10~pred.score.dh, ts_f, algorithm=2, direction = "<")
plot(roc.score.f, print.thres = TRUE, print.auc = TRUE) #0.731, 0.013, 0.536, 0.815
#(with km(death endpoint): 0.724, .001, .536, .815 )

roc.score.hdl.f <- roc(status_dh~pred.score.hdl.dh, ts_f, algorithm=2, direction = "<")
plot(roc.score.hdl.f, print.thres = TRUE, print.auc = TRUE) #0.767, 0.056, 0.789, 0.611
#(with km(death endpoint): 0.756, .0009, .789, .606 )

#Rees-FRS
roc.frs.f <- roc(status_dh~pred.frs, ts_f, algorithm=2, direction = "<")
auc(roc.frs.f) #0.7757
plot(roc.frs.f, print.thres = TRUE) #0.045, 0.679, 0.764

roc.frs.of.f <- roc(status_dh~pred.frs.of, ts_f, algorithm=2, direction = "<")
auc(roc.frs.of.f) #0.7595
plot(roc.frs.of.f, print.thres = TRUE) #0.041, 0.616, 0.787

#Original ------- 
#frs - dh - d
roc.frs_ori_f <- roc(status_dh_10~frs_ori, ts_f, algorithm=2, direction = "<")
auc(roc.frs_ori_f) #0.742
plot(roc.frs_ori_f, print.thres= TRUE) #0.084, 0.580, 0.845

roc.frsbmi_ori_f <- roc(status_dh_10~frsbmi_ori, ts_f, algorithm=2, direction = "<")
auc(roc.frsbmi_ori_f) #0.739
plot(roc.frsbmi_ori_f, print.thres= TRUE) #0.35, 0.595, 0.780

#score - dh - d
roc.score_local_f <- roc(status_dh_10~score_local, ts_f, algorithm=2, direction = "<")
auc(roc.score_local_f) #0.7241 (with death endpoint: 0.7125)
plot(roc.score_local_f, print.thres= TRUE) #0.027, 0.526, .826 (0.002, 0.526, 0.826)

roc.score_local_hdl_f <- roc(status_dh_10~score_local_hdl, ts_f, algorithm=2, direction = "<")
auc(roc.score_local_hdl_f) #0.6992 (with death endpoint: 0.7127)
plot(roc.score_local_hdl_f, print.thres= TRUE) #0.117, 0.722, .595 (0.002, 0.622, 0.727)



# male ---------
ts_m <- hdl_imp1 %>% 
  filter(gender =="Male" & SjvReg != "Ume<U+00E5> region" ) %>% 
  mutate(status_dh_10 = ifelse(time_dh > 10, 0, status_dh))%>%
  left_join(score_df) %>% 
  left_join(mean_df) %>% 
  mutate(
    pred.score.dh = 1- sur.prob^predict(hdl_models$score_dh[[2]], newdata = ., type = "risk"),
    pred.score.hdl.dh = 1- sur.prob^predict(hdl_models$score_dh_hdl[[2]], newdata = ., type = "risk"),
    pred.frs = risk.bh.predict(coxph.object = hdl_models$frs[[2]], newdata = ., time.p = 10),
    pred.frs.of = risk.bh.predict(coxph.object = hdl_models$frs_office[[2]], newdata = ., time.p = 10),
    a = 3.06117*log(alder) + 1.12370*log(skol_mgdl) - 0.93263*log(hdl_mgdl) + 1.93303*Tx_ht + 1.99881*log(sbt) + 0.65451*smoke + 0.57367*dm,
    b = 3.06117*3.8560 + 1.12370*5.3420 - 0.93263*3.7686 + 1.93303*4.3544 + 1.99881*0.5019 + 0.65451*0.3522 + 0.57367*0.0650,  
    frs_ori = 1 - 0.88936^exp(a - b),
    a1 =  3.11296*log(alder) + 0.79277*log(bmi) + 1.92672*Tx_ht + 1.85508*log(sbt) + 0.70953*smoke + 0.53160*dm, 
    b1 =  3.11296*3.8560 + 0.79277*26.49642 + 1.92672*4.3544 + 1.85508*0.5019 + 0.70953*0.3522 + 0.53160*0.0650, 
    frsbmi_ori = 1 - 0.88431^exp(a1 - b1),
    a_s = 0.0185592*sbt+0.004595*skol+0.7213618*smoke,
    score_local = 1 - sur.prob^exp(a_s-mean_sum),
    a_s_hdl = log(0.78)*hdl+log(1.33)*skol+log(1.37)*sbt+log(1.80)*smoke,
    score_local_hdl = 1 - sur.prob^exp(a_s_hdl-mean_sum_hdl)) 

roc.score.m <- roc(status_dh_10~pred.score.dh, ts_m, algorithm=2, direction = "<")
plot(roc.score.m, print.thres = TRUE, print.auc=TRUE) #0.726, 0.029,0.527,0.820
#(with death endpoint: 0.725, .004, .526, .821 )

roc.score.hdl.m <- roc(status_dh_10~pred.score.hdl.dh, ts_m, algorithm=2, direction = "<")
plot(roc.score.hdl.m, print.thres = TRUE, print.auc=TRUE) #0.732, 0.029, 0.517, 0.836
#(with death endpoint: 0.730, .004, .516, .836 )

roc.frs.m <- roc(status_dh_10~pred.frs, ts_m, algorithm=2, direction = "<")
plot(roc.frs.m, print.thres = TRUE, print.auc=TRUE) #0.747, .074,.537,.830
roc.frs.of.m <- roc(status_dh~pred.frs.of, ts_m, algorithm=2, direction = "<")
plot(roc.frs.of.m, print.thres = TRUE, print.auc=TRUE) #0.721,.086,.577,.763

#Original
roc.frs_ori_m <- roc(status_dh~frs_ori, ts_m, algorithm=2, direction = "<")
plot(roc.frs_ori_m, print.thres = TRUE, print.auc=TRUE) # .722, .173, .573, .796

roc.frsbmi_ori_m <- roc(status_dh~frsbmi_ori, ts_m, algorithm=2, direction = "<")
plot(roc.frsbmi_ori_m, print.thres = TRUE, print.auc=TRUE) #.698, 0.0, .530, .813

roc.score_local_m <- roc(status_dh~score_local, ts_m, algorithm=2, direction = "<")
plot(roc.score_local_m, print.thres = TRUE, print.auc=TRUE) #.722, .073, .594, .755
#(with death endpoint: 0.721, .010, .594, .755 )

roc.score_local_hdl_m <- roc(status_dh~score_local_hdl, ts_m, algorithm=2, direction = "<")
plot(roc.score_local_hdl_m, print.thres = TRUE, print.auc=TRUE) #.629, .073, .594, .755
#(with death endpoint: 0.628, .044, .717, .483 )

#Calibration (GND chi sqaure)-------------
#The GND for SCORE is not so trustworthy, since it only capture ONE baseline hazard
library(ResourceSelection)
hoslem.test(ts_f$status_dh_10,ts_f$pred.frs, 5 )
hoslem.test(ts_f$status_dh_10,ts_f$pred.frs.of, 5 )
hoslem.test(ts_f$status_dh_10,ts_f$pred.score.dh, 5 )
hoslem.test(ts_f$status_dh_10,ts_f$pred.score.hdl.dh, 5 )
hoslem.test(ts_f$status_dh_10,ts_f$frs_ori, 5 )
hoslem.test(ts_f$status_dh_10,ts_f$frsbmi_ori, 5 )
hoslem.test(ts_f$status_dh_10,ts_f$score_local, 5 )
hoslem.test(ts_f$status_dh_10,ts_f$score_local_hdl, 5 )

hoslem.test(ts_m$status_dh_10,ts_m$pred.frs, 5 )
hoslem.test(ts_m$status_dh_10,ts_m$pred.frs.of, 5 )
hoslem.test(ts_m$status_dh_10,ts_m$pred.score.dh, 5 )
hoslem.test(ts_m$status_dh_10,ts_m$pred.score.hdl.dh, 5 )
hoslem.test(ts_m$status_dh_10,ts_m$frs_ori, 5 )
hoslem.test(ts_m$status_dh_10,ts_m$frsbmi_ori, 5 )
hoslem.test(ts_m$status_dh_10,ts_m$score_local, 5 )
hoslem.test(ts_m$status_dh_10,ts_m$score_local_hdl, 5 )

ts_f$de <- as.numeric(cut2(ts_f$pred.frs, g=5))
#ts_f$de2 <- ifelse(round(ts_f$pred.frs,2)<=0.04, 1,2)
table(ts_f$de, ts_f$status_dh)
gnd.frs.f <- GND.calib(pred=ts_f$pred.frs, tvar=ts_f$time_dh, out=ts_f$status_dh_10, 
                     cens.t=10, groups = ts_f$de, adm.cens=10) ; gnd.frs.f



ts_f$de <- as.numeric(cut2(ts_f$pred.frs.of, g=5))
#ts_f$de2 <- ifelse(round(ts_f$pred.frs,2)<=0.04, 1,2)
table(ts_f$de, ts_f$status_dh)
gnd.frs.f <- GND.calib(pred=ts_f$pred.frs.of, tvar=ts_f$time_dh, out=ts_f$status_dh_10, 
                       cens.t=10, groups = ts_f$de, adm.cens=10) ; gnd.frs.f

ts_f$de <- as.numeric(cut2(ts_f$frs_ori, g=5))
#ts_f$de2 <- ifelse(round(ts_f$pred.frs,2)<=0.04, 1,2)
table(ts_f$de, ts_f$status_dh)
gnd.frs.f <- GND.calib(pred=ts_f$frs_ori, tvar=ts_f$time_dh, out=ts_f$status_dh_10, 
                       cens.t=10, groups = ts_f$de, adm.cens=10) ; gnd.frs.f

ts_f$de <- as.numeric(cut2(ts_f$frsbmi_ori, g=5))
#ts_f$de2 <- ifelse(round(ts_f$pred.frs,2)<=0.04, 1,2)
table(ts_f$de, ts_f$status_dh)
gnd.frs.f <- GND.calib(pred=ts_f$frsbmi_ori, tvar=ts_f$time_dh, out=ts_f$status_dh_10, 
                       cens.t=10, groups = ts_f$de, adm.cens=10) ; gnd.frs.f

ts_f_4 <- ts_f %>% 
  filter(alder == 60)

ts_f_4$de <- as.numeric(cut2(ts_f_4$pred.score.dh, g=5))
#ts_f$de2 <- ifelse(round(ts_f$pred.frs,2)<=0.04, 1,2)
table(ts_f_4$de, ts_f_4$status_dh_10)
ts_f_4$de <- ifelse(ts_f_4$de==5,4,ts_f_4$de)
ts_f_4$de <- ifelse(ts_f_4$de==1,2,ts_f_4$de)
gnd.frs.f <- GND.calib(pred=ts_f_4$pred.score.dh, tvar=ts_f_4$alder, out=ts_f_4$status_dh_10, 
                       cens.t=10, groups = ts_f_4$de, adm.cens=10) ; gnd.frs.f
#Age 40 = 10.6,  0.01389747 
#Age 50 - 50.4,  2.928905e-10
#age 60 - 1931.914,    0.000 


ts_f$de <- as.numeric(cut2(ts_f$pred.score.dh, g=5))
#ts_f$de2 <- ifelse(round(ts_f$pred.frs,2)<=0.04, 1,2)
table(ts_f$de, ts_f$status_dh)
gnd.frs.f <- GND.calib(pred=ts_f$pred.score.dh, tvar=ts_f$alder, out=ts_f$status_dh, 
                       cens.t=60, groups = ts_f$de, adm.cens=60) ; gnd.frs.f


ts_f$de <- as.numeric(cut2(ts_f$pred.score.hdl.dh, g=5))
#ts_f$de2 <- ifelse(round(ts_f$pred.frs,2)<=0.04, 1,2)
table(ts_f$de, ts_f$status_dh)
gnd.frs.f <- GND.calib(pred=ts_f$pred.score.hdl.dh, tvar=ts_f$alder, out=ts_f$status_dh, 
                       cens.t=60, groups = ts_f$de, adm.cens=60) ; gnd.frs.f

ts_f$de <- as.numeric(cut2(ts_f$score_local, g=5))
#ts_f$de2 <- ifelse(round(ts_f$pred.frs,2)<=0.04, 1,2)
table(ts_f$de, ts_f$status_dh)
gnd.frs.f <- GND.calib(pred=ts_f$score_local, tvar=ts_f$alder, out=ts_f$status_dh, 
                       cens.t=60, groups = ts_f$de, adm.cens=60) ; gnd.frs.f

ts_f$de <- as.numeric(cut2(ts_f$score_local_hdl, g=5))
#ts_f$de2 <- ifelse(round(ts_f$pred.frs,2)<=0.04, 1,2)
table(ts_f$de, ts_f$status_dh)
gnd.frs.f <- GND.calib(pred=ts_f$score_local_hdl, tvar=ts_f$alder, out=ts_f$status_dh, 
                       cens.t=60, groups = ts_f$de, adm.cens=60) ; gnd.frs.f


#Male
ts_m$de <- as.numeric(cut2(ts_m$pred.frs, g=5))
table(ts_m$de, ts_m$status_dh)
gnd.m <- GND.calib(pred=ts_m$pred.frs, tvar=ts_m$time_dh, out=ts_m$status_dh, 
                   cens.t=10, groups = ts_m$de, adm.cens=10) ; gnd.m

ts_m$de <- as.numeric(cut2(ts_m$pred.frs.of, g=5))
table(ts_m$de, ts_m$status_dh)
gnd.m <- GND.calib(pred=ts_m$pred.frs.of, tvar=ts_m$time_dh, out=ts_m$status_dh, 
                   cens.t=10, groups = ts_m$de, adm.cens=10) ; gnd.m

ts_m$de <- as.numeric(cut2(ts_m$frs_ori, g=5))
table(ts_m$de, ts_m$status_dh)
gnd.m <- GND.calib(pred=ts_m$frs_ori, tvar=ts_m$time_dh, out=ts_m$status_dh_10, 
                   cens.t=10, groups = ts_m$de, adm.cens=10) ; gnd.m

ts_m$de <- as.numeric(cut2(ts_m$frsbmi_ori, g=5))
table(ts_m$de, ts_m$status_dh)
gnd.m <- GND.calib(pred=ts_m$frsbmi_ori, tvar=ts_m$time_dh, out=ts_m$status_dh_10, 
                   cens.t=10, groups = ts_m$de, adm.cens=10) ; gnd.m

ts_m$de <- as.numeric(cut2(ts_m$pred.score.dh, g=5))
table(ts_m$de, ts_m$status_dh)
gnd.m <- GND.calib(pred=ts_m$pred.score.dh, tvar=ts_m$alder, out=ts_m$status_dh_10, 
                   cens.t=60, groups = ts_m$de, adm.cens=60) ; gnd.m

ts_m$de <- as.numeric(cut2(ts_m$pred.score.hdl.dh, g=5))
table(ts_m$de, ts_m$status_dh)
gnd.m <- GND.calib(pred=ts_m$pred.score.hdl.dh, tvar=ts_m$alder, out=ts_m$status_dh_10, 
                   cens.t=60, groups = ts_m$de, adm.cens=60) ; gnd.m

ts_m$de <- as.numeric(cut2(ts_m$score_local, g=5))
table(ts_m$de, ts_m$status_dh)
gnd.m <- GND.calib(pred=ts_m$score_local, tvar=ts_m$alder, out=ts_m$status_dh, 
                   cens.t=60, groups = ts_m$de, adm.cens=60) ; gnd.m

ts_m$de <- as.numeric(cut2(ts_m$score_local_hdl, g=5))
table(ts_m$de, ts_m$status_dh)
gnd.m <- GND.calib(pred=ts_m$score_local_hdl, tvar=ts_m$alder, out=ts_m$status_dh, 
                   cens.t=60, groups = ts_m$de, adm.cens=60) ; gnd.m

#Calibration plots ------------- 
#riskRegression + pec... so slow... 
x <- Score(list("FRS"=hdl_models$frs_dh[[1]],
                "SCORE" = hdl_models$score_dh[[1]]),data=ts_f,formula=Surv(time_dh, status_dh)~1,plots="cali",times=10)
plotCalibration(x,ylim=c(0,0.3),xlim=c(0,0.5))

#AUC plots -------------- 
#male
par(pty="s")
plot.roc(roc.17.h.m, lty = 5, main="Receiver operator curves for males FRS origin definition", col="royalblue")
plot.roc(roc.17.h.m.frs, add = TRUE, col="red" )
testobj1 <- roc.test(roc.17.h.m, roc.17.h.m.frs, method = "delong")
text(0.42, 0.55, labels=paste("Recalibrated AUC:", format.pval(auc(roc.17.h.m), digits = 3)), adj=c(0, .5), cex = 0.8)
text(0.42, 0.45, labels=paste("Original FRS AUC:", format.pval(auc(roc.17.h.m.frs), digits = 3)), adj=c(0, .5), cex = 0.8)
text(0.42, 0.35, labels=paste("DeLong's test:", format.pval(testobj1$p.value, digits = 3)), adj=c(0, .5), cex = 0.8)
legend("bottomright", legend=c("Recalibrated", "Origin FRS"), col=c("royalblue", "red"), lty = c(5,1),
       lwd=2, cex=0.8, bty = "n", title.adj=0.15, text.font = 4)

#female
par(pty="s")
plot.roc(roc.17.h.f, lty = 5, main="Receiver operator curves for females FRS origin definition", col="royalblue")
plot.roc(roc.17.h.f.frs, add = TRUE, col="red" )
testobj1 <- roc.test(roc.17.h.f, roc.17.h.f.frs, method = "delong")
text(0.42, 0.55, labels=paste("Recalibrated AUC:", format.pval(auc(roc.17.h.f), digits = 3)), adj=c(0, .5), cex = 0.8)
text(0.42, 0.45, labels=paste("Original FRS AUC:", format.pval(auc(roc.17.h.f.frs), digits = 3)), adj=c(0, .5), cex = 0.8)
text(0.42, 0.35, labels=paste("DeLong's test:", format.pval(testobj1$p.value, digits = 3)), adj=c(0, .5), cex = 0.8)
legend("bottomright", legend=c("Recalibrated", "Origin FRS"), col=c("royalblue", "red"), lty = c(5,1),
       lwd=2, cex=0.8, bty = "n", title.adj=0.15, text.font = 4)


#examine PH assumptions----------
ph.17.h.f <- cox.zph(hdl_models$score_dh[[1]])
ph.17.h.f #All follows PH assumptions

ph.17.h.m <- cox.zph(hdl_17_models$cox[[2]])
ph.17.h.m #log(skol) violates the assumption, plot to examine
plot(ph.17.h.m, var = "log(skol)")


#regressions for population mean -----------
#exp_score_factors <- hdl_imp1 %>% 
#  group_by(gender) %>% 
#  do(sbt_mod = lm(sbt~poly(alder,2),data = .),
#     tc_mod = lm(skol_mmol~poly(alder,2),data=.),
#     smk_mod = glm(smoke~poly(alder,2), family = "binomial", data = .))


#score_fact_f <- km_sur_prob %>% 
#  filter(gender=="Female") %>%
#  mutate(
#    mean_sbp = predict(exp_score_factors$sbt_mod[[1]], newdata = .),
#    mean_tc = predict(exp_score_factors$tc_mod[[1]], newdata = .),
#    smk_prev = predict(exp_score_factors$smk_mod[[1]], newdata = ., type = "response"),
#    mean_sum = 0.0185592*mean_sbp+0.004595*mean_tc+0.7213618*smk_prev)

# plot km curves by gender and age ------------
library(ggfortify)

mod_f <- survfit(Surv(time_dh, status_dh) ~ alder, data = filter(hdl_imp1,gender=="Female"))
mod_m <- survfit(Surv(time_dh, status_dh) ~ alder, data = filter(hdl_imp1,gender=="Male"))
# extract results to a data.frame
res_f <- mod_f %>% 
  fortify() %>% 
  mutate(gender = "Female")
res_m <- mod_m %>% 
  fortify() %>% 
  mutate(gender = "Male")

res <- rbind(res_f, res_m)
res <- res %>% 
  mutate(gen_alder = interaction(as.factor(gender), strata))
str(res)

ggplot(data = filter(res, time <= 10), aes(x = time, y = surv, color = gen_alder)) +
  geom_line() + 
  # plot censor marks
  #  geom_point(aes(shape = factor(ifelse(n.censor >= 1, 1, NA)))) + 
  # format censor shape as "+"
  scale_shape_manual(values = 3) + 
  # hide censor legend 
  guides(shape = "none") 
#summary stat---------
sum_df <- hdl_imp1 %>% 
  group_by(SjvReg, gender) %>% 
  summarize_at( 
    c("status_d", "status_dh","time_dh","alder", "sbt", "skol", "hdl", "smoke", "dm", "Tx_ht", "Tx_statin", "bmi"), mean) %>% 
  mutate(
    status_d = status_d*100,
    status_dh = status_dh*100,
    Tx_ht = Tx_ht*100,
    Tx_statin = Tx_statin*100,
    dm = dm*100,
    smoke = smoke*100) %>% 
  ungroup() %>%
  mutate_if(is.numeric, round,1) %>% 
  t() 

table(hdl_imp1$gender, hdl_imp1$Tx_statin)

sum_df_sd <- hdl_imp1 %>% 
  group_by(SjvReg, gender) %>% 
  summarize_at( 
    c("time_dh","alder", "sbt", "skol", "hdl","bmi"), sd) %>% 
  ungroup() %>%
  mutate_if(is.numeric, round,1) %>% 
  t() 
#write.table(sum_df, file = "sum_df.txt", sep = "&", quote = FALSE, row.names=TRUE, col.names = FALSE)
#write.table(sum_df_sd, file = "sum_df_sd.txt", sep = "&", quote = FALSE, row.names=TRUE, col.names = FALSE)
