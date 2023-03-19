# this file provides the codes for GAM main analyses(O3-adjusted model)
# For single-factor model, just delete the O3 covariant 

library(tidyverse); library(dlnm);
library(lubridate); library(tsModel);library(ggsci)
library(mgcv); library(scales); library(splines)

################################## Load data ###################################

load("hko_23ap_ha_0516.rda")
df <- hko_ap_ha_0516 %>% select(date,DOW,Holiday,Tmean,Humidity,co,fsp,no2,o3h8max,copd)
names(df) <- c("date","dow","holiday","temp","rh","co","fsp","no2","o3h8max","copd")

dfA <- df

dfA$t <- 1:length(dfA$date)

str(dfA)
summary(dfA)
dim(dfA)
names(dfA)

data=dfA

################################## Multi plt ###################################

FTjust <- function(dis, plt, coplt1, coplt2,df_t, coplt3){
  
  cb.plt = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt1 = crossbasis(data[[coplt1]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt2 = crossbasis(data[[coplt2]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  
  plist <- ifelse(plt == coplt3, "cb.plt+cb.coplt1+cb.coplt2", "cb.plt+cb.coplt1+cb.coplt2+cb.coplt3")
  
  eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
        s(t,k=",df_t,"*12,fx=T,bs='cr')+

        as.factor(dow)+as.factor(holiday),
 
        family=quasipoisson(link = 'log'), scale=-1,
        control=gam.control(epsilon=0.0000001,maxit=200),
        na.action= na.omit, data=data)", sep='')))
  
  
  iqr<-IQR(data[[plt]],na.rm=TRUE)
  eval(parse(text=paste("pred <- crosspred(cb.plt, model,
                        at=iqr, cumul=TRUE)",sep='')))
  matRRfit <- data.frame(fit = "matfit", dis = dis, plt = plt, iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "mathigh", dis = dis, plt = plt, iqr=iqr,
                          coplt = coplt3,df_t=df_t,
                          pred$matfit+3.29*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matlow", dis = dis, plt = plt,iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit-3.29*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}

plist <- list(dis = c('copd'),
              plt = c("co","fsp","no2","o3h8max"),
              coplt1 = c('temp'),
              coplt2 = c("rh"),
              df_t=c(7),
              coplt3 = c("o3h8max"))

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)


bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

output <- bbA %>% mutate(
  sig=as.factor(ifelse(mathigh<=0,1,ifelse(matlow>=0,1,0))))

names(output) <- c("dis","plt","iqr","coplt","df_t","lagA","beta","h","l","lag","Significance")
output

################################## temp###################################

FTjust <- function(dis, plt, coplt1, coplt2,df_t, coplt3){
  
  # cb.coplt1 = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.plt  = crossbasis(data[[coplt1]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt2 = crossbasis(data[[coplt2]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  
  plist <- ifelse(plt==coplt3,"cb.plt+cb.coplt2","cb.plt+cb.coplt2+cb.coplt3")
    
  eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
        s(t,k=",df_t,"*12,fx=T,bs='cr')+

        as.factor(dow)+as.factor(holiday),
 
        family=quasipoisson(link = 'log'), scale=-1,
        control=gam.control(epsilon=0.0000001,maxit=200),
        na.action= na.omit, data=data)", sep='')))
  
  
  iqr<-IQR(data[[plt]],na.rm=TRUE)
  eval(parse(text=paste("pred <- crosspred(cb.plt, model,
                        at=iqr, cumul=TRUE)",sep='')))
  matRRfit <- data.frame(fit = "matfit", dis = dis, plt = plt, iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "mathigh", dis = dis, plt = plt, iqr=iqr,
                          coplt = coplt3,df_t=df_t,
                          pred$matfit+3.29*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matlow", dis = dis, plt = plt,iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit-3.29*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}

plist <- list(dis = c('copd'),
              plt = c('temp'),
              coplt1 = c('temp'),
              coplt2 = c("rh"),
              df_t=c(7),
              coplt3 = c("o3h8max"))

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)


bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

output <- bbA %>% mutate(
  sig=as.factor(ifelse(mathigh<=0,1,ifelse(matlow>=0,1,0))))
output 

################################## rh###################################

FTjust <- function(dis, plt, coplt1, coplt2,df_t, coplt3){
  
  # cb.coplt1 = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt2 = crossbasis(data[[coplt1]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.plt = crossbasis(data[[coplt2]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  
  plist <- ifelse(plt==coplt3,"cb.plt+cb.coplt2","cb.plt+cb.coplt2+cb.coplt3")
  
  eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
        s(t,k=",df_t,"*12,fx=T,bs='cr')+

        as.factor(dow)+as.factor(holiday),
 
        family=quasipoisson(link = 'log'), scale=-1,
        control=gam.control(epsilon=0.0000001,maxit=200),
        na.action= na.omit, data=data)", sep='')))
  
  
  iqr<-IQR(data[[plt]],na.rm=TRUE)
  eval(parse(text=paste("pred <- crosspred(cb.plt, model,
                        at=iqr, cumul=TRUE)",sep='')))
  matRRfit <- data.frame(fit = "matfit", dis = dis, plt = plt, iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "mathigh", dis = dis, plt = plt, iqr=iqr,
                          coplt = coplt3,df_t=df_t,
                          pred$matfit+3.29*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matlow", dis = dis, plt = plt,iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit-3.29*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}

plist <- list(dis = c('copd'),
              plt = c('rh'),
              coplt1 = c('temp'),
              coplt2 = c("rh"),
              df_t=c(7),
              coplt3 = c("o3h8max"))

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)


bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

output <- bbA %>% mutate(
  sig=as.factor(ifelse(mathigh<=0,1,ifelse(matlow>=0,1,0))))

names(output) <- c("dis","plt","iqr","coplt","df_t","lagA","beta","h","l","lag","Significance")

output
