###*******************************************************************************************
###
### R code for analysis in:
###
### "Using parametric g-computation for time-to-event data and distributed lag models 
###  to identify critical exposure windows for preterm birth: An illustrative example 
###  using PM2.5 in a retrospective birth cohort based in Eastern Massachusetts (2011-2016)"
###
###*******************************************************************************************

###*************
### N: Notes ###
###*************

#' @param id individual id
#' @param preterm preterm status indicator (1: preterm, 0: not)
#' @param ga_Week gestational week of delivery
#' @param year year of conception
#' @param season  season of conception
#' @param pmlag average PM2.5 in a specific gestational week (e.g., pmlag01 is PM2.5 in gestational week 1, pmlag02 is PM2.5 in gestational week 2 etc.)
#' @param ndlag average NO2 in a specific gestational week (e.g., ndlag01 is NO2 in gestational week 1, ndlag02 is NO2 in gestational week 2 etc.)
#' @param templag average temperature in a specific gestational week (e.g., templag01 is temperature in gestational week 1, templag02 is temperature in gestational week 2 etc.)


###******************************
### Required Packages & Setup ###
###******************************

# install required packages
requiredPackages = c('tidyverse','dlnm','splitstackshape','boot','svMisc','ggplot2')
for (p in requiredPackages) {
  if (!require(p, character.only = TRUE))
    install.packages(p)
  library(p, character.only = TRUE)
}

# Create dataframes to store results
pred_obs <- data.frame(id=character(), time=numeric(), p=numeric(), week_intervention=numeric(), intervention=numeric()) # for conditional probabilities under no intervention
pred_intervention <- data.frame(id=character(), time=numeric(), p=numeric(), week_intervention=numeric(), intervention=numeric()) # for conditional probabilities under intervention

###****************************
### G-Survival-DLM Analysis ###
###****************************

### STEP 0: READ IN DATASETS

df <- readRDS("C:/Users/Michael/Dropbox/BIDMC/Preterm/Projects/PM/Revision/To Submit/Code/df_mock.rds")
pm.week <- readRDS("C:/Users/Michael/Dropbox/BIDMC/Preterm/Projects/PM/Revision/To Submit/Code/pm_week.rds")
nd.week <- readRDS("C:/Users/Michael/Dropbox/BIDMC/Preterm/Projects/PM/Revision/To Submit/Code/nd_week.rds")
temp.week <- readRDS("C:/Users/Michael/Dropbox/BIDMC/Preterm/Projects/PM/Revision/To Submit/Code/temp_week.rds")

for(i in 28:36){ # loop through all the risk sets
  
  ### STEP 1: FORMAT DATA
  
  # 1a. transform dataset into long format (person-week for each row)
  df.surv <- df %>% mutate(survtime = ifelse(ga_weeks >= 37, 37, ga_weeks))
  df.surv <- expandRows(df.surv, "survtime", drop=F)
  df.surv <- df.surv %>%
    mutate(time = sequence(rle(id)$lengths)-1,
           event = ifelse(time == survtime-1 & preterm==1, 1, 0))

  # 1b. split long dataset into separate datasets for each risk set
  df.surv <- df.surv %>% filter(time == i) %>% left_join(pm.week, by="id") %>% left_join(nd.week, by="id") %>% left_join(temp.week, by="id")
  
  # select weekly exposures up until the risk set
  pm.hist <- df.surv %>% select(pmlag01:paste("pmlag",i,sep="")) 
  nd.hist <- df.surv %>% select(ndlag01:paste("ndlag",i,sep=""))
  temp.hist <- df.surv %>% select(templag01:paste("templag",i,sep=""))
  
  ### STEP 2: NO INTERVENTION
  
  # determine degrees of freedom for the lag response (can be determined by model fit like AIC or subject matter knowledge)
  n_dfs <- 4
  
  # create crossbasis (need to manually do this to allow back-predict absolute probabilities)
  cb.pm.obs <- crossbasis(pm.hist, argvar = list("lin"), arglag = list(df=n_dfs)) %>% data.frame()
  cb.nd.obs <- crossbasis(nd.hist, argvar = list("lin"), arglag = list(df=n_dfs)) %>% data.frame()
  cb.temp.obs <- crossbasis(temp.hist, argvar = list("lin"), arglag = list(df=n_dfs)) %>% data.frame()
  
  names(cb.pm.obs) <- paste0('cb.', "pm", '.', names(cb.pm.obs))
  names(cb.nd.obs) <- paste0('cb.', "nd", '.', names(cb.nd.obs))
  names(cb.temp.obs) <- paste0('cb.', "temp", '.', names(cb.temp.obs))
  
  # rejoin manually created crossbases to original dataset
  df.surv <- df.surv %>% bind_cols(cb.pm.obs) %>% bind_cols(cb.nd.obs) %>% bind_cols(cb.temp.obs)
  
  # Step 2a. fit a distributed lag model to the observed data 
  gf.model <- glm(event ~ 
                    cb.pm.v1.l1 + cb.pm.v1.l2 + cb.pm.v1.l3 + cb.pm.v1.l4 +   # manually specified crossbasis for PM2.5
                    cb.nd.v1.l1 + cb.nd.v1.l2 + cb.nd.v1.l3 + cb.nd.v1.l4 +   # manually specified crossbasis for NO2
                    cb.temp.v1.l1 + cb.temp.v1.l2 + cb.temp.v1.l3 + cb.temp.v1.l4, # manually specified crossbasis for temperature
                  family = binomial(), 
                  data = df.surv)

  # Step 2b. predict conditional probability of not being born 
  pred0 <- expandRows(df, count=37, count.is.col=F) # need to create new dataset for all risk sets regardless of when individual was born
  pred0$time <- rep(seq(0, 36), nrow(df))
  pred0 <- pred0 %>% filter(time == i)
  pred0 <- pred0 %>% left_join(pm.week, by="id") %>% left_join(nd.week, by="id") %>% left_join(temp.week, by="id") 
  
  # create new crossbasis based on complete exposure history of everyone
  pm.week.cf0 <- pred0 %>% select(pmlag01:paste("pmlag",i,sep=""))
  nd.week.cf0 <- pred0 %>% select(ndlag01:paste("ndlag",i,sep=""))
  temp.week.cf0 <- pred0 %>% select(templag01:paste("templag",i,sep=""))
  
  cb.pm.obs0 <- crossbasis(pm.week.cf0, argvar = list("lin"), arglag = list(df=n_dfs)) %>% data.frame()
  cb.nd.obs0 <- crossbasis(nd.week.cf0, argvar = list("lin"), arglag = list(df=n_dfs)) %>% data.frame()
  cb.temp.obs0 <- crossbasis(temp.week.cf0, argvar = list("lin"), arglag = list(df=n_dfs)) %>% data.frame()
  
  names(cb.pm.obs0) <- paste0('cb.', "pm", '.', names(cb.pm.obs0))
  names(cb.nd.obs0) <- paste0('cb.', "nd", '.', names(cb.nd.obs0))
  names(cb.temp.obs0) <- paste0('cb.', "temp", '.', names(cb.temp.obs0))
  
  # predict the conditional probabilities
  pred0 <- pred0 %>% bind_cols(cb.pm.obs0) %>% bind_cols(cb.nd.obs0) %>% bind_cols(cb.temp.obs0)
  pred0$p_noevent <- 1 - predict(gf.model, pred0, type="response", se.fit=F)
  pred0 <- pred0 %>% select(id, time, p_noevent) %>% mutate(week_intervention = 0, intervention = 0)
  pred_obs <- rbind(pred_obs, pred0) %>% arrange(id, time)
  
  ### STEP 3: SIMULATED INTERVENTION
  
  for(k in 1:i){ # for each risk set loop through all the exposure weeks up until the risk set
    
    # create new prediction dataset for interventions
    pred <- expandRows(df, count=37, count.is.col=F)  # again, need to create new dataset for all risk sets regardless of when individual was born
    pred$time <- rep(seq(0, 36), nrow(df))
    pred <- pred %>% filter(time == i)
    pred <- pred %>% left_join(pm.week, by="id") %>% left_join(nd.week, by="id") %>% left_join(temp.week, by="id") 
    
    pm.week.cf1 <- pm.week.cf0
    nd.week.cf <- nd.week.cf0
    temp.week.cf <- temp.week.cf0
    
    # intervention: 20% reduction in each gestational week
    pm.week.cf1[,k] <- pm.week.cf1[,k]*0.8 
    
    # create new crossbasis based on counterfactual exposure history under intervention
    cb.pm.obs1 <- crossbasis(pm.week.cf1, argvar = list("lin"), arglag = list(df=n_dfs)) %>% data.frame()
    cb.nd.obs1 <- crossbasis(nd.week.cf, argvar = list("lin"), arglag = list(df=n_dfs)) %>% data.frame()
    cb.temp.obs1 <- crossbasis(temp.week.cf, argvar = list("lin"), arglag = list(df=n_dfs)) %>% data.frame()
    
    names(cb.pm.obs1) <- paste0('cb.', "pm", '.', names(cb.pm.obs1))
    names(cb.nd.obs1) <- paste0('cb.', "nd", '.', names(cb.nd.obs1))
    names(cb.temp.obs1) <- paste0('cb.', "temp", '.', names(cb.temp.obs1))
    
    pred1 <- pred %>% bind_cols(cb.pm.obs1)
    pred1 <- pred1 %>% bind_cols(cb.nd.obs1)
    pred1 <- pred1 %>% bind_cols(cb.temp.obs1)

    # predict conditional probability under intervention
    pred1$p_noevent <- 1 - predict(gf.model, pred1, type="response", se.fit=F)
    pred1 <- pred1 %>% select(id, time, p_noevent) %>% mutate(week_intervention = k, intervention = 1)
    pred_intervention <- rbind(pred_intervention, pred1) %>% arrange(id, week_intervention, time)
  }
  progress(i,36)
}

### STEPS 2C-D & 3C-E: POOLING and STANDARDIZATION

# combine conditional probabilities under no intervention and under intervention
pred1 <- rbind(pred_obs, pred_intervention)

# These are for when the intervention week does not occur during follow-up
pred_average1 <- pred1 %>% 
  filter(week_intervention < 29) %>% 
  group_by(id, week_intervention) %>% 
  mutate(surv = cumprod(p_noevent)) %>% # pool conditional probability of no event
  na.omit() 

pred_average1 <- pred_average1 %>%
  group_by(week_intervention, time) %>% 
  summarise(surv1 = mean(surv, na.rm=T)) %>% # standardize by week of intervention and risk set
  mutate(risk = 1 - surv1)

# These are for intervention weeks that occur during the follow-up

# need to combine observed conditional probabilities up until the week of the intervention, then use conditional probabilites under intervention
pred29_36_all1 <- data.frame(id=character(), time=numeric(), p=numeric(), week_intervention=numeric())
for(m in 29:36) {
  pred29_36_1 <- rbind((pred_obs %>% filter(time < m)), (pred1 %>% filter(week_intervention == m))) %>% 
    arrange(id, time) %>% 
    mutate(week_intervention = m)
  pred29_36_all1 <- rbind(pred29_36_all1, pred29_36_1)
}

# pool
pred_average1_29_36 <- pred29_36_all1 %>% 
  group_by(id, week_intervention) %>% 
  mutate(surv = cumprod(p_noevent)) %>% 
  na.omit()  

# standardize
pred_average1_29_36 <- pred_average1_29_36 %>%
  group_by(week_intervention, time) %>% 
  summarise(surv1 = mean(surv)) %>% 
  mutate(risk = 1 - surv1)

# Combine all the intervention weeks together (i.e., interventions prior to and then during)
pred_average_all1 <- rbind(pred_average1, pred_average1_29_36) %>% filter(time >= week_intervention)
pred_average_obs1 <- pred_average_all1 %>% filter(week_intervention == 0) %>% ungroup() %>% select(-week_intervention) %>% rename(risk0 = risk)

### STEP 4: CALCULATE RISK DIFFERENCE

# calculate the risk differnce for each risk set and week of intervention
pred_average_all1 <- pred_average_all1 %>% left_join(pred_average_obs1, by = "time") %>% mutate(rd = risk - risk0)

# Calculate the cumulative risk difference (i.e., what is the risk of preterm overall if we intervened in every week)
sum((pred_average_all1 %>% filter(time == 36))$rd)

### STEP 5 (ADDITIONAL STEPS)

# Creating Figure 3 (Risk Difference for each risk set and week of intervention)
ggplot(pred_average_all1, aes(x=week_intervention, y=(time), fill=rd)) +
  geom_raster() +
  scale_fill_gradient2(high="blue", mid="white",low="red",midpoint=0) +
  labs(x="Gestational Week of Intervention", y="Risk Set (Gestational Week of Follow-Up)", fill="Risk Difference") +
  scale_x_continuous(breaks=seq(0, 36, 5)) +
  scale_y_reverse(breaks=seq(28, 36, 1)) +
  theme(panel.spacing.y=unit(1.5, "lines"),
        panel.spacing.x=unit(1.25, "lines"),
        panel.border=element_rect(colour="black", linewidth=0.2, fill=NA),
        panel.background=element_blank(),
        strip.text=element_text(size=16),
        axis.title=element_text(size=16))

# Creating Figure S3 (Risks under No Intervention)
pred_average_all1 %>% filter(week_intervention == 0) %>% ungroup() %>% select(time, risk) %>% mutate(risk = round(risk,4)) %>% 
  ggplot() +
  geom_line(aes(x=time, y=risk)) +
  labs(x="Gestational Weeks", y="Risk of Birth") +
  scale_x_continuous(limits = c(28, 36), breaks=seq(28,36,1)) +
  scale_y_continuous(limits = c(0, 0.1), breaks=seq(0,0.1,0.01)) +
  coord_cartesian(xlim=c(28,36)) +
  theme_bw() +
  theme(panel.spacing.y=unit(1.5,"lines"),
        panel.spacing.x=unit(1.25,"lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position=c(0.15, 0.87))
