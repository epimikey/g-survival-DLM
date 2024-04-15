rm(list=ls())
setwd("C:/Users/Michael/Dropbox/BIDMC/Preterm/Projects/PM/")

# install required packages
requiredPackages = c('tidyverse','mgcv','lubridate','dlnm','data.table','survival','survminer','splitstackshape','boot','svMisc','splines','ggplot2','ggsci')
for (p in requiredPackages) {
  if (!require(p, character.only = TRUE))
    install.packages(p)
  library(p, character.only = TRUE)
}

df <- readRDS("C:/Users/Michael/Dropbox/BIDMC/Preterm/Data/df_fetGrow.rds") %>%
  select(c(baby_mrn, mom_mrn, ga_weeks, doc, doy, conceptyear, season, mage, meduc, mrace, ins, 
           parity_first, sex, smoking, adi_natrank)) %>%
  mutate(preterm = ifelse(ga_weeks < 37, 1, 0),
         harm = sin((2*pi*doy)/(365)) + cos((2*pi*doy)/(365)),
         doy = yday(doc),
         ga_weeks = ceiling(ga_weeks)) %>%
  filter(ga_weeks >= 28) %>%
  distinct() %>%
  na.omit()

pm.week <- readRDS("C:/Users/Michael/Dropbox/BIDMC/Fetal Growth/Data/exposures/pm_week.rds") %>% rename(time=week) %>% select(-pmsd.week)
nd.week <- readRDS("C:/Users/Michael/Dropbox/BIDMC/Fetal Growth/Data/exposures/nd_week.rds") %>% rename(time=week) %>% select(-ndsd.week)
temp.week <- readRDS("C:/Users/Michael/Dropbox/BIDMC/Fetal Growth/Data/exposures/tempPrism_week.rds") %>% rename(time=week) %>% select(-tsd)

# reformat exposure
pm.week <- pm.week %>% filter(time <= 36) %>% 
  mutate(week = str_pad(1:n(), 2, pad="0"), lag = paste("pmlag",week,sep="")) %>%
  select(-c(week, time)) %>% spread(lag, pm.week) %>%
  filter(baby_mrn %in% df$baby_mrn)

nd.week <- nd.week %>% filter(time <= 36) %>% 
  mutate(week = str_pad(1:n(), 2, pad="0"), lag = paste("ndlag", week,sep="")) %>%
  select(-c(week, time)) %>% spread(lag, nd.week) %>%
  filter(baby_mrn %in% df$baby_mrn)

temp.week <- temp.week %>% filter(time <= 36) %>% 
  mutate(week = str_pad(1:n(), 2, pad="0"), lag = paste("templag",week,sep="")) %>%
  select(-c(week, time)) %>% spread(lag, tmean) %>%
  filter(baby_mrn %in% df$baby_mrn)


###****************************
### 0: Read in Data & Setup ###
###****************************

preterm.boot <- function(data, indices) {
  
  d <- data[indices,]
  p$tick()$print()  # update progress bar
  
  # store results
  pred_obs <- data.frame(baby_mrn=character(), time=numeric(), p=numeric(), week_intervention=numeric(), intervention=numeric())
  pred_intervention <- data.frame(baby_mrn=character(), time=numeric(), p=numeric(), week_intervention=numeric(), intervention=numeric())
  
  for(i in 28:36){
    pm_week <- pm.week %>% select(1:i+1)
    nd_week <- nd.week %>% select(1:i+1)
    temp_week <- temp.week %>% select(1:i+1)
    
    # STEP 0: creation of person-week data
    df.surv <- d %>% mutate(survtime = ifelse(ga_weeks >= 37, 37, ga_weeks))
    df.surv <- expandRows(df.surv, "survtime", drop=F)
    df.surv$time <- sequence(rle(df.surv$baby_mrn)$lengths)-1
    df.surv$event <- ifelse(df.surv$time==df.surv$survtime-1 &
                              df.surv$preterm==1, 1, 0)
    df.surv <- df.surv %>% filter(time == i) %>% left_join(pm_week, by="baby_mrn") %>% left_join(nd_week, by="baby_mrn") %>% left_join(temp_week, by="baby_mrn") 
    
    meteor <- df.surv %>% select(-c(baby_mrn:event))
    pm.hist <- meteor[,1:(i+1)]
    nd.hist <- meteor[,(i+1):(2*i)]
    temp.hist <- meteor[,(2*i+1):(3*i)]
    
    # set AIC
    aic <- 3
    
    # create crossbasis
    cb.pm.obs <- crossbasis(pm.hist, argvar = list("lin"), arglag = list(df=aic)) %>% data.frame()
    cb.nd.obs <- crossbasis(nd.hist, argvar = list("lin"), arglag = list(df=aic)) %>% data.frame()
    cb.temp.obs <- crossbasis(temp.hist, argvar = list("lin"), arglag = list(df=aic)) %>% data.frame()
    
    names(cb.pm.obs) <- paste0('cb.', "pm", '.', names(cb.pm.obs))
    df.surv <- df.surv %>% bind_cols(cb.pm.obs)
    
    names(cb.nd.obs) <- paste0('cb.', "nd", '.', names(cb.nd.obs))
    df.surv <- df.surv %>% bind_cols(cb.nd.obs)
    
    names(cb.temp.obs) <- paste0('cb.', "temp", '.', names(cb.temp.obs))
    df.surv <- df.surv %>% bind_cols(cb.temp.obs)
    
    gf.model <- glm(event ~ as.factor(conceptyear) + as.factor(season) +
                      poly(mage,2) + as.factor(meduc) + as.factor(mrace) + as.factor(parity_first) +
                      as.factor(sex) + as.factor(ins) + poly(adi_natrank,2) + 
                      cb.pm.v1.l1 + cb.pm.v1.l2 + cb.pm.v1.l3 +  
                      cb.nd.v1.l1 + cb.nd.v1.l2 + cb.nd.v1.l3 +  
                      cb.temp.v1.l1 + cb.temp.v1.l2 + cb.temp.v1.l3 + , 
                    family = binomial(), 
                    data = df.surv)
    
    # predict observed
    pred0 <- expandRows(d, count=37, count.is.col=F)
    pred0$time <- rep(seq(0, 36), nrow(d))
    pred0 <- pred0 %>% filter(time == i)
    pred0 <- pred0 %>% left_join(pm_week, by="baby_mrn") %>% left_join(nd_week, by="baby_mrn") %>% left_join(temp_week, by="baby_mrn") 
    
    meteor.cf0 <- pred0 %>% select(-c(baby_mrn:time))
    pm_week.cf0 <- meteor.cf0[,1:(i)]
    nd_week.cf0 <- meteor.cf0[,(i+1):(2*i)]
    temp_week.cf0 <- meteor.cf0[,(2*i+1):(3*i)]
    
    cb.pm.obs0 <- crossbasis(pm_week.cf0, argvar = list("lin"), arglag = list(df=aic)) %>% data.frame()
    names(cb.pm.obs0) <- paste0('cb.', "pm", '.', names(cb.pm.obs0))
    pred0 <- pred0 %>% bind_cols(cb.pm.obs0)
    
    cb.nd.obs0 <- crossbasis(nd_week.cf0, argvar = list("lin"), arglag = list(df=aic)) %>% data.frame()
    names(cb.nd.obs0) <- paste0('cb.', "nd", '.', names(cb.nd.obs0))
    pred0 <- pred0 %>% bind_cols(cb.nd.obs0)
    
    cb.temp.obs0 <- crossbasis(temp_week.cf0, argvar = list("lin"), arglag = list(df=4)) %>% data.frame()
    names(cb.temp.obs0) <- paste0('cb.', "temp", '.', names(cb.temp.obs0))
    pred0 <- pred0 %>% bind_cols(cb.temp.obs0)
    
    pred0$p_noevent <- 1 - predict(gf.model, pred0, type="response", se.fit=F)
    pred0 <- pred0 %>% select(baby_mrn, time, p_noevent) %>% mutate(week_intervention = 0, intervention = 0)
    
    pred_obs <- rbind(pred_obs, pred0) %>% arrange(baby_mrn, time)
    
    
    # simulate interventions
    for(k in 1:i){
      
      # create prediction dataset (needs to include those lost in the previous risk set as well)
      pred <- expandRows(d, count=37, count.is.col=F)
      pred$time <- rep(seq(0, 36), nrow(d))
      pred <- pred %>% filter(time == i)
      pred <- pred %>% left_join(pm_week, by="baby_mrn") %>% left_join(nd_week, by="baby_mrn") %>% left_join(temp_week, by="baby_mrn") 
      
      pm_week.cf1 <- pm_week.cf0
      pm_week.cf2 <- pm_week.cf0
      
      nd_week.cf <- nd_week.cf0
      temp_week.cf <- temp_week.cf0
      
      # intervention 1: 20% reduction
      pm_week.cf1[,k] <- pm_week.cf1[,k]*0.8 
      
      # intervention 2: set to week-specific median
      med_week <- rep(quantile(pm_week.cf2[,k], probs=0.5, na.rm=T), dim(df)[1])
      pm_week.cf2 <- cbind(pm_week.cf2, med_week)
      pm_week.cf2[,k] <- ifelse(pm_week.cf2[,k] > pm_week.cf2$med_week, pm_week.cf2$med_week, pm_week.cf2[,k])
      pm_week.cf2 <- pm_week.cf2 %>% select(-med_week)
      
      # predict with crossbasis
      cb.pm.obs1 <- crossbasis(pm_week.cf1, argvar = list("lin"), arglag = list(df=aic)) %>% data.frame()
      names(cb.pm.obs1) <- paste0('cb.', "pm", '.', names(cb.pm.obs1))
      pred1 <- pred %>% bind_cols(cb.pm.obs1)
      
      cb.pm.obs2 <- crossbasis(pm_week.cf2, argvar = list("lin"), arglag = list(df=aic)) %>% data.frame()
      names(cb.pm.obs2) <- paste0('cb.', "pm", '.', names(cb.pm.obs2))
      pred2 <- pred %>% bind_cols(cb.pm.obs2)
      
      cb.nd.obs1 <- crossbasis(nd_week.cf, argvar = list("lin"), arglag = list(df=aic)) %>% data.frame()
      names(cb.nd.obs1) <- paste0('cb.', "nd", '.', names(cb.nd.obs1))
      pred1 <- pred1 %>% bind_cols(cb.nd.obs1)
      pred2 <- pred2 %>% bind_cols(cb.nd.obs1)
      
      cb.temp.obs1 <- crossbasis(temp_week.cf, argvar = list("lin"), arglag = list(df=aic)) %>% data.frame()
      names(cb.temp.obs1) <- paste0('cb.', "temp", '.', names(cb.temp.obs1))
      pred1 <- pred1 %>% bind_cols(cb.temp.obs1)
      pred2 <- pred2 %>% bind_cols(cb.temp.obs1)
      
      # predict under observed
      pred1$p_noevent <- 1 - predict(gf.model, pred1, type="response", se.fit=F)
      pred2$p_noevent <- 1 - predict(gf.model, pred2, type="response", se.fit=F)
      
      pred1 <- pred1 %>% select(baby_mrn, time, p_noevent) %>% mutate(week_intervention = k, intervention = 1)
      pred2 <- pred2 %>% select(baby_mrn, time, p_noevent) %>% mutate(week_intervention = k, intervention = 2)
      
      pred_intervention <- rbind(pred_intervention, pred1, pred2) %>% arrange(baby_mrn, week_intervention, time)
    }
    progress(i,36)
  }
  
  pred_all <- rbind(pred_obs, pred_intervention)
  remove(pred_intervention)
  
  pred1 <- pred_all %>% filter(intervention != 2)
  pred2 <- pred_all %>% filter(intervention != 1)
  
  ### These are for when the week intervention are before study entry
  pred_average1 <- pred1 %>% filter(week_intervention < 29) %>% group_by(baby_mrn, week_intervention) %>% mutate(surv = cumprod(p_noevent)) %>% na.omit() %>%
    group_by(week_intervention, time) %>% summarise(surv1 = mean(surv, na.rm=T)) %>% mutate(risk = 1 - surv1)
  
  pred_average2 <- pred2 %>% filter(week_intervention < 29) %>% group_by(baby_mrn, week_intervention) %>% mutate(surv = cumprod(p_noevent)) %>% na.omit() %>%
    group_by(week_intervention, time) %>% summarise(surv1 = mean(surv, na.rm=T)) %>% mutate(risk = 1 - surv1)
  
  ### These are for intervention that occur during the follow-up
  pred29_36_all1 <- data.frame(baby_mrn=character(), time=numeric(), p=numeric(), week_intervention=numeric())
  for(m in 29:36) {
    pred29_36_1 <- rbind((pred_obs %>% filter(time < m)), (pred1 %>% filter(week_intervention == m))) %>% arrange(baby_mrn, time) %>% mutate(week_intervention = m)
    pred29_36_all1 <- rbind(pred29_36_all1, pred29_36_1)
  }
  
  pred29_36_all2 <- data.frame(baby_mrn=character(), time=numeric(), p=numeric(), week_intervention=numeric())
  for(m in 29:36) {
    pred29_36_2 <- rbind((pred_obs %>% filter(time < m)), (pred2 %>% filter(week_intervention == m))) %>% arrange(baby_mrn, time) %>% mutate(week_intervention = m)
    pred29_36_all2 <- rbind(pred29_36_all2, pred29_36_2)
  }
  
  pred_average1_29_36 <- pred29_36_all1 %>% group_by(baby_mrn, week_intervention) %>% mutate(surv = cumprod(p_noevent)) %>% na.omit() %>% 
    group_by(week_intervention, time) %>% summarise(surv1 = mean(surv)) %>% mutate(risk = 1 - surv1)
  
  pred_average2_29_36 <- pred29_36_all2 %>% group_by(baby_mrn, week_intervention) %>% mutate(surv = cumprod(p_noevent)) %>% na.omit() %>% 
    group_by(week_intervention, time) %>% summarise(surv1 = mean(surv)) %>% mutate(risk = 1 - surv1)
  
  ### Combine the two (i.e., interventions prior to and then during)
  pred_average_all1 <- rbind(pred_average1, pred_average1_29_36) %>% filter(time >= week_intervention)
  pred_average_all2 <- rbind(pred_average2, pred_average2_29_36) %>% filter(time >= week_intervention)
  
  pred_average_obs1 <- pred_average_all1 %>% filter(week_intervention == 0) %>% ungroup() %>% select(-week_intervention) %>% rename(risk0 = risk)
  pred_average_obs2 <- pred_average_all2 %>% filter(week_intervention == 0) %>% ungroup() %>% select(-week_intervention) %>% rename(risk0 = risk)
  
  pred_average_all1 <- pred_average_all1 %>% left_join(pred_average_obs1, by = "time") %>% mutate(rd = risk - risk0)
  pred_average_all2 <- pred_average_all2 %>% left_join(pred_average_obs2, by = "time") %>% mutate(rd = risk - risk0)
  
  return(c(sum((pred_average_all1 %>% filter(time == 36))$rd), sum((pred_average_all2 %>% filter(time == 36))$rd)))
}

tot_rep <- 30
p <- progress_estimated(tot_rep+1)
set.seed(2345)
survival.std.results <- boot(data = df, statistic = preterm.boot, R = tot_rep)

boot.original <- data.frame(survival.std.results$t0)
boot.std <- data.frame(survival.std.results$t)

boot.iteration <- data.frame(cbind(original=data.frame(survival.std.results$t0), t(survival.std.results$t)))

# saveRDS(boot.iteration, "C:/Users/Michael/Dropbox/BIDMC/Preterm/Projects/PM/Revision/boot_seed0123_20.rds")
saveRDS(boot.iteration, "C:/Users/Michael/Dropbox/BIDMC/Preterm/Projects/PM/Revision/boot_seed02345_30.rds")

