library(tidyverse); library(dlnm);
library(lubridate); library(tsModel);library(ggsci)
library(mgcv); library(scales); library(splines)

################################## Load data ###################################

load("hko_23ap_ha_0516.rda")
df <- hko_ap_ha_0516 %>% select(date,DOW,Holiday,Tmean,Humidity,co,fsp,no,no2,nox,o3,o3h8max,rsp,so2,copd)
names(df) <- c("date","dow","holiday","temp","rh","co","fsp","no","no2","nox","o3","o3h8max","rsp","so2","copd")

dfA <- df

dfA$t <- 1:length(dfA$date)

str(dfA)
summary(dfA)
dim(dfA)
names(dfA)

data=dfA

################################## Single + Multi plt ###################################

FTjust <- function(dis, plt, coplt1, coplt2,df_t, coplt3){
  ### coplt1 = temp
  ### coplt2 = rh
  ### coplt3 = real coplt
  ### df_t = df of time variable t
  # cb.plt = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
  # cb.coplt1 = crossbasis(data[[coplt1]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt2 = crossbasis(data[[coplt2]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  
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
  matRRfit <- data.frame(fit = "matRRfit", dis = dis, plt = plt, iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "matRRhigh", dis = dis, plt = plt, iqr=iqr,
                          coplt = coplt3,df_t=df_t,
                          pred$matfit+1.96*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matRRlow", dis = dis, plt = plt,iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit-1.96*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}

plist <- list(dis = c('copd'),
              plt = c("co","fsp","no2","o3h8max"),
              coplt1 = c('temp'),
              coplt2 = c("rh"),
              df_t=c(7,12),
              coplt3 = c("co","fsp","no2","o3h8max"))

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)


bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

# print(as_tibble(bbA), n=Inf)
dat1 <- bbA

bbA <- dat1
output <- bbA %>% mutate(
  sig=as.factor(ifelse(matRRhigh<=0,1,ifelse(matRRlow>=0,1,0))))

names(output) <- c("dis","plt","iqr","coplt","df_t","lagA","beta","h","l","lag","Significance")

output$lag <- factor(output$lag)
output$Significance <- factor(output$Significance,levels = c(0,1),labels = c("Non-sig","Sig"))


output$Coplt <- output$coplt
output$plt <- factor(output$plt,levels = c("o3h8max","co","fsp","no2"),labels = c("O3","CO","PM2.5","NO2"))
output$Coplt <- factor(output$Coplt,levels = c("o3h8max","co","fsp","no2"),labels = c("O3","CO","PM2.5","NO2"))

plt.labs <- c("O3","CO","PM2.5","NO2")
names(plt.labs) <- c("O3","CO","PM2.5","NO2")
Coplt.labs <- paste("coplt_",c("O3","CO","PM2.5","NO2"),sep = "")
names(Coplt.labs) <- c("O3","CO","PM2.5","NO2")
P2_2 <- ggplot(output, aes(lag,beta,ymin = l, ymax = h)) +
  geom_hline(yintercept = 0,linetype='dashed') +
  geom_errorbar(aes(group = plt, col = plt, width=0), 
                size=1, position = position_dodge((width=0.5))) +
  scale_color_jama()+
  geom_point(aes(group = plt, col = plt,shape=Significance),size = 2, 
             position = position_dodge(width=0.5),,fontface = "bold") +
  # scale_color_manual(values = c('#0272B7',"red")) +
  scale_shape_manual(values = c(16,8)) +
  # geom_ribbon(aes(x=lag,ymin =l, ymax = h),
  #             color=NA,fill = "gray10",alpha=0.2) +
  # geom_linerange(aes(x=lag,ymin =l, ymax = h,col=plt2),
  #                size=0.8,color='gray10',alpha=0.5)+
  # facet_grid(Coplt~plt,scales='free',labeller = "label_parsed") +
  facet_grid(Coplt~plt,scales='free',labeller = labeller(Coplt = Coplt.labs, plt = plt.labs)) +
  ylab(expression(paste(beta,"(Risk estimates air pollutants → COPD)", sep = "")))+
  xlab("Lag")+
  theme_bw() +
  # mytheme1
  theme(axis.text.x=element_text(size = 10, color= "black"),
        axis.text.y=element_text(size = 10, color= "black"),
        axis.title.x=element_text(size=11),
        axis.title.y=element_text(size=11),
        axis.ticks.length=unit(0.2,'cm'),
        axis.line = element_line(colour = "black"),
        legend.position=c("bottom"),
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.box="vertical",
        legend.margin=margin(-15,unit="pt"),
        # legend.box.margin=margin(10,10,10,10),
        legend.box.spacing = margin(15.5),
        legend.background = element_rect(fill="transparent"),
        strip.background = element_rect(
          color = "white", fill = "white"),
        panel.border = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        panel.grid = element_blank(),
        # panel.border = element_rect(fill=NA,color="black", size=0.9, linetype="solid"),
        plot.title = element_text(size = 18, hjust=0.5))
P2_2
# ggsave("COPD_gam_coplt_sensitive.jpg", width = 10, height = 5.5, dpi = 300)


################################## o3_adjusted t df sensitive ###################################

FTjust <- function(dis, plt, coplt1, coplt2,df_t, coplt3){
  ### coplt1 = temp
  ### coplt2 = rh
  ### coplt3 = real coplt
  ### df_t = df of time variable t
  # cb.plt = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
  # cb.coplt1 = crossbasis(data[[coplt1]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt2 = crossbasis(data[[coplt2]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  
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
  matRRfit <- data.frame(fit = "matRRfit", dis = dis, plt = plt, iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "matRRhigh", dis = dis, plt = plt, iqr=iqr,
                          coplt = coplt3,df_t=df_t,
                          pred$matfit+3.29*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matRRlow", dis = dis, plt = plt,iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit-3.29*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}

plist <- list(dis = c('copd'),
              plt = c("co","fsp","no2","o3h8max"),
              coplt1 = c('temp'),
              coplt2 = c("rh"),
              df_t=c(5:12),
              coplt3 = c("o3h8max"))

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)


bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

# print(as_tibble(bbA), n=Inf)
dat1 <- bbA

bbA <- dat1
output <- bbA %>% mutate(
  sig=as.factor(ifelse(matRRhigh<=0,1,ifelse(matRRlow>=0,1,0))))

names(output) <- c("dis","plt","iqr","coplt","df_t","lagA","beta","h","l","lag","Significance")

output$lag <- factor(output$lag)
output$Significance <- factor(output$Significance,levels = c(0,1),labels = c("Non-sig","Sig"))


output$Coplt <- output$coplt
output$plt <- factor(output$plt,levels = c("o3h8max","co","fsp","no2"),labels = c("O3","CO","PM2.5","NO2"))
output$df_t <- factor(output$df_t)

# plt.labs <- c("O3","CO","PM2.5","NO2")
# names(plt.labs) <- c("O3","CO","PM2.5","NO2")
# Coplt.labs <- paste("coplt_",c("O3","CO","PM2.5","NO2"),sep = "")
# names(Coplt.labs) <- c("O3","CO","PM2.5","NO2")
P2_2 <- ggplot(output, aes(lag,beta,ymin = l, ymax = h)) +
  geom_hline(yintercept = 0,linetype='dashed') +
  geom_errorbar(aes(group = plt, col = plt, width=0), 
                size=1, position = position_dodge((width=0.5))) +
  scale_color_jama()+
  geom_point(aes(group = plt, col = plt,shape=Significance),size = 2, 
             position = position_dodge(width=0.5),,fontface = "bold") +
  # scale_color_manual(values = c('#0272B7',"red")) +
  scale_shape_manual(values = c(16,8)) +
  # geom_ribbon(aes(x=lag,ymin =l, ymax = h),
  #             color=NA,fill = "gray10",alpha=0.2) +
  # geom_linerange(aes(x=lag,ymin =l, ymax = h,col=plt2),
  #                size=0.8,color='gray10',alpha=0.5)+
  # facet_grid(Coplt~plt,scales='free',labeller = "label_parsed") +
  facet_grid(df_t~plt,scales='free',labeller = "label_parsed") +
  ylab(expression(paste(beta,"(Risk estimates air pollutants → COPD)", sep = "")))+
  xlab("Lag")+
  theme_bw() +
  # mytheme1
  theme(axis.text.x=element_text(size = 10, color= "black"),
        axis.text.y=element_text(size = 10, color= "black"),
        axis.title.x=element_text(size=11),
        axis.title.y=element_text(size=11),
        axis.ticks.length=unit(0.2,'cm'),
        axis.line = element_line(colour = "black"),
        legend.position=c("bottom"),
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.box="vertical",
        legend.margin=margin(-15,unit="pt"),
        # legend.box.margin=margin(10,10,10,10),
        legend.box.spacing = margin(15.5),
        legend.background = element_rect(fill="transparent"),
        strip.background = element_rect(
          color = "white", fill = "white"),
        panel.border = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        panel.grid = element_blank(),
        # panel.border = element_rect(fill=NA,color="black", size=0.9, linetype="solid"),
        plot.title = element_text(size = 18, hjust=0.5))
P2_2
ggsave("COPD_gam_dft_sensitive.jpg", width = 12, height = 8, dpi = 300)







################################## temp###################################

FTjust <- function(dis, plt, coplt1, coplt2,df_t, coplt3){
  ### coplt1 = temp
  ### coplt2 = rh
  ### coplt3 = real coplt
  ### df_t = df of time variable t
  # cb.plt = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
  # cb.coplt1 = crossbasis(data[[coplt1]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt2 = crossbasis(data[[coplt2]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  
  cb.coplt1 = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
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
  matRRfit <- data.frame(fit = "matRRfit", dis = dis, plt = plt, iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "matRRhigh", dis = dis, plt = plt, iqr=iqr,
                          coplt = coplt3,df_t=df_t,
                          pred$matfit+1.96*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matRRlow", dis = dis, plt = plt,iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit-1.96*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}
#"co","fsp","no2","o3h8max"
plist <- list(dis = c('copd'),
              plt = c('temp'),
              coplt1 = c('temp'),
              coplt2 = c("rh"),
              df_t=c(7,12),
              coplt3 = c('temp',"co","fsp","no2","o3h8max"))

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)


bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

# print(as_tibble(bbA), n=Inf)
dat2 <- bbA

bbA <- dat2
output <- bbA %>% mutate(
  sig=as.factor(ifelse(matRRhigh<=0,1,ifelse(matRRlow>=0,1,0))))

names(output) <- c("dis","plt","iqr","coplt","df_t","lagA","beta","h","l","lag","Significance")

output$lag <- factor(output$lag)
output$Significance <- factor(output$Significance,levels = c(0,1),labels = c("Non-sig","Sig"))


output$Coplt <- output$coplt
# output$plt <- factor(output$plt,levels = c("o3h8max","co","fsp","no2"),labels = c("O3","CO","PM2.5","NO2"))
output$df_t <- factor(output$df_t)

# plt.labs <- c("O3","CO","PM2.5","NO2")
# names(plt.labs) <- c("O3","CO","PM2.5","NO2")
# Coplt.labs <- paste("coplt_",c("O3","CO","PM2.5","NO2"),sep = "")
# names(Coplt.labs) <- c("O3","CO","PM2.5","NO2")
P2_2 <- ggplot(output, aes(lag,beta,ymin = l, ymax = h)) +
  geom_hline(yintercept = 0,linetype='dashed') +
  geom_errorbar(aes(group = plt, col = plt, width=0), 
                size=1, position = position_dodge((width=0.5))) +
  scale_color_jama()+
  geom_point(aes(group = plt, col = plt,shape=Significance),size = 2, 
             position = position_dodge(width=0.5),,fontface = "bold") +
  # scale_color_manual(values = c('#0272B7',"red")) +
  scale_shape_manual(values = c(16,8)) +
  # geom_ribbon(aes(x=lag,ymin =l, ymax = h),
  #             color=NA,fill = "gray10",alpha=0.2) +
  # geom_linerange(aes(x=lag,ymin =l, ymax = h,col=plt2),
  #                size=0.8,color='gray10',alpha=0.5)+
  # facet_grid(Coplt~plt,scales='free',labeller = "label_parsed") +
  facet_grid(~coplt,scales='free',labeller = "label_parsed") +
  ylab(expression(paste(beta,"(Risk estimates air pollutants → COPD)", sep = "")))+
  xlab("Lag")+
  theme_bw() +
  # mytheme1
  theme(axis.text.x=element_text(size = 10, color= "black"),
        axis.text.y=element_text(size = 10, color= "black"),
        axis.title.x=element_text(size=11),
        axis.title.y=element_text(size=11),
        axis.ticks.length=unit(0.2,'cm'),
        axis.line = element_line(colour = "black"),
        legend.position=c("bottom"),
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.box="vertical",
        legend.margin=margin(-15,unit="pt"),
        # legend.box.margin=margin(10,10,10,10),
        legend.box.spacing = margin(15.5),
        legend.background = element_rect(fill="transparent"),
        strip.background = element_rect(
          color = "white", fill = "white"),
        panel.border = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        panel.grid = element_blank(),
        # panel.border = element_rect(fill=NA,color="black", size=0.9, linetype="solid"),
        plot.title = element_text(size = 18, hjust=0.5))
P2_2
# ggsave("COPD_gam_dft_sensitive.jpg", width = 12, height = 8, dpi = 300)

################################## rh###################################

FTjust <- function(dis, plt, coplt1, coplt2,df_t, coplt3){
  ### coplt1 = temp
  ### coplt2 = rh
  ### coplt3 = real coplt
  ### df_t = df of time variable t
  # cb.plt = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
  # cb.coplt1 = crossbasis(data[[coplt1]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt2 = crossbasis(data[[coplt2]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  
  cb.coplt1 = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
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
  matRRfit <- data.frame(fit = "matRRfit", dis = dis, plt = plt, iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "matRRhigh", dis = dis, plt = plt, iqr=iqr,
                          coplt = coplt3,df_t=df_t,
                          pred$matfit+1.96*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matRRlow", dis = dis, plt = plt,iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit-1.96*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}
#"co","fsp","no2","o3h8max"
plist <- list(dis = c('copd'),
              plt = c('rh'),
              coplt1 = c('temp'),
              coplt2 = c("rh"),
              df_t=c(7,12),
              coplt3 = c("rh","co","fsp","no2","o3h8max"))

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)


bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

# print(as_tibble(bbA), n=Inf)
dat3 <- bbA

bbA <- dat3
output <- bbA %>% mutate(
  sig=as.factor(ifelse(matRRhigh<=0,1,ifelse(matRRlow>=0,1,0))))

names(output) <- c("dis","plt","iqr","coplt","df_t","lagA","beta","h","l","lag","Significance")

output$lag <- factor(output$lag)
output$Significance <- factor(output$Significance,levels = c(0,1),labels = c("Non-sig","Sig"))


output$Coplt <- output$coplt
# output$plt <- factor(output$plt,levels = c("o3h8max","co","fsp","no2"),labels = c("O3","CO","PM2.5","NO2"))
output$df_t <- factor(output$df_t)

# plt.labs <- c("O3","CO","PM2.5","NO2")
# names(plt.labs) <- c("O3","CO","PM2.5","NO2")
# Coplt.labs <- paste("coplt_",c("O3","CO","PM2.5","NO2"),sep = "")
# names(Coplt.labs) <- c("O3","CO","PM2.5","NO2")
P2_2 <- ggplot(output, aes(lag,beta,ymin = l, ymax = h)) +
  geom_hline(yintercept = 0,linetype='dashed') +
  geom_errorbar(aes(group = plt, col = plt, width=0), 
                size=1, position = position_dodge((width=0.5))) +
  scale_color_jama()+
  geom_point(aes(group = plt, col = plt,shape=Significance),size = 2, 
             position = position_dodge(width=0.5),,fontface = "bold") +
  # scale_color_manual(values = c('#0272B7',"red")) +
  scale_shape_manual(values = c(16,8)) +
  # geom_ribbon(aes(x=lag,ymin =l, ymax = h),
  #             color=NA,fill = "gray10",alpha=0.2) +
  # geom_linerange(aes(x=lag,ymin =l, ymax = h,col=plt2),
  #                size=0.8,color='gray10',alpha=0.5)+
  # facet_grid(Coplt~plt,scales='free',labeller = "label_parsed") +
  facet_grid(~coplt,scales='free',labeller = "label_parsed") +
  ylab(expression(paste(beta,"(Risk estimates air pollutants → COPD)", sep = "")))+
  xlab("Lag")+
  theme_bw() +
  # mytheme1
  theme(axis.text.x=element_text(size = 10, color= "black"),
        axis.text.y=element_text(size = 10, color= "black"),
        axis.title.x=element_text(size=11),
        axis.title.y=element_text(size=11),
        axis.ticks.length=unit(0.2,'cm'),
        axis.line = element_line(colour = "black"),
        legend.position=c("bottom"),
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.box="vertical",
        legend.margin=margin(-15,unit="pt"),
        # legend.box.margin=margin(10,10,10,10),
        legend.box.spacing = margin(15.5),
        legend.background = element_rect(fill="transparent"),
        strip.background = element_rect(
          color = "white", fill = "white"),
        panel.border = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        panel.grid = element_blank(),
        # panel.border = element_rect(fill=NA,color="black", size=0.9, linetype="solid"),
        plot.title = element_text(size = 18, hjust=0.5))
P2_2
# ggsave("COPD_gam_dft_sensitive.jpg", width = 12, height = 8, dpi = 300)
gam_res <- rbind(dat1,dat2,dat3)
save(gam_res, file = "0313gam_res95.rda")



























# 
# 
# 
# FTjust <- function(dis, plt, coplt1, coplt2, l1 = 0,l2 = 3){
#   cb.plt = crossbasis(data[[plt]], lag=3, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
#   cb.coplt1 = crossbasis(data[["temp"]], lag=3, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
#   cb.coplt2 = crossbasis(data[["rh"]], lag=3, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
#   
#   # plist <- ifelse(plt == coplt1, "cb.plt", "cb.plt+cb.coplt1")
#   plist <-  "cb.plt+cb.coplt1+cb.coplt2"
#   # s(rh,k=6+1,fx=T,bs='cr')+
#   
#   eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
#         s(t,k=7*12,fx=T,bs='cr')+
# 
#         as.factor(dow)+as.factor(holiday),
#  
#         family=quasipoisson(link = 'log'), scale=-1,
#         control=gam.control(epsilon=0.0000001,maxit=200),
#         na.action= na.omit, data=data)", sep='')))
#   
#   attr(cb.plt,"lag") <- c(l1,l2)
#   iqr<-IQR(data[[plt]],na.rm=TRUE)
#   eval(parse(text=paste("pred <- crosspred(cb.plt, model,
#                         at=iqr, cumul=TRUE)",sep='')))
#   matRRfit <- data.frame(fit = "matRRfit", dis = dis, plt = plt, iqr=iqr,
#                          coplt = coplt1,
#                          pred$matfit, stringsAsFactors = FALSE)
#   matRRhigh <- data.frame(fit = "matRRhigh", dis = dis, plt = plt, iqr=iqr,
#                           coplt = coplt1,
#                           pred$matfit+2.58*pred$matse, stringsAsFactors = FALSE)
#   matRRlow <- data.frame(fit = "matRRlow", dis = dis, plt = plt,iqr=iqr,
#                          coplt = coplt1,
#                          pred$matfit-2.58*pred$matse, stringsAsFactors = FALSE)
#   out <- bind_rows(matRRfit, matRRhigh, matRRlow)
#   out
# }





























# 
# 
# 
# FTjust <- function(dis, plt, coplt1, coplt2, l1 = 0,l2 = 3){
#   cb.plt = crossbasis(data[[plt]], lag=3, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
#   cb.coplt1 = crossbasis(data[["temp"]], lag=3, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
#   cb.coplt2 = crossbasis(data[["rh"]], lag=3, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
#   
#   # plist <- ifelse(plt == coplt1, "cb.plt", "cb.plt+cb.coplt1")
#   plist <-  "cb.plt+cb.coplt1+cb.coplt2"
#   # s(rh,k=6+1,fx=T,bs='cr')+
#   
#   eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
#         s(t,k=7*12,fx=T,bs='cr')+
# 
#         as.factor(dow)+as.factor(holiday),
#  
#         family=quasipoisson(link = 'log'), scale=-1,
#         control=gam.control(epsilon=0.0000001,maxit=200),
#         na.action= na.omit, data=data)", sep='')))
#   
#   attr(cb.plt,"lag") <- c(l1,l2)
#   iqr<-IQR(data[[plt]],na.rm=TRUE)
#   eval(parse(text=paste("pred <- crosspred(cb.plt, model,
#                         at=iqr, cumul=TRUE)",sep='')))
#   matRRfit <- data.frame(fit = "matRRfit", dis = dis, plt = plt, iqr=iqr,
#                          coplt = coplt1,
#                          pred$matfit, stringsAsFactors = FALSE)
#   matRRhigh <- data.frame(fit = "matRRhigh", dis = dis, plt = plt, iqr=iqr,
#                           coplt = coplt1,
#                           pred$matfit+2.58*pred$matse, stringsAsFactors = FALSE)
#   matRRlow <- data.frame(fit = "matRRlow", dis = dis, plt = plt,iqr=iqr,
#                          coplt = coplt1,
#                          pred$matfit-2.58*pred$matse, stringsAsFactors = FALSE)
#   out <- bind_rows(matRRfit, matRRhigh, matRRlow)
#   out
# }


















################################## temp and rh ###################################

FTjust <- function(dis, plt, coplt1, coplt2,df_t, coplt3){
  ### coplt1 = temp
  ### coplt2 = rh
  ### coplt3 = real coplt
  ### df_t = df of time variable t
  # cb.plt = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
  # cb.coplt1 = crossbasis(data[[coplt1]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt2 = crossbasis(data[[coplt2]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  
  cb.coplt1 = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.plt= crossbasis(data[[coplt1]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt2 = crossbasis(data[[coplt2]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  
  plist <- "cb.plt+cb.coplt1+cb.coplt2"
  # plist <- "cb.plt+cb.coplt2"
  
  eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
        s(t,k=",df_t,"*12,fx=T,bs='cr')+

        as.factor(dow)+as.factor(holiday),
 
        family=quasipoisson(link = 'log'), scale=-1,
        control=gam.control(epsilon=0.0000001,maxit=200),
        na.action= na.omit, data=data)", sep='')))
  
  
  iqr<-IQR(data[[coplt1]],na.rm=TRUE)
  eval(parse(text=paste("pred <- crosspred(cb.plt, model,
                        at=iqr, cumul=TRUE)",sep='')))
  matRRfit <- data.frame(fit = "matRRfit", dis = dis, plt = plt, iqr=iqr,
                         coplt = coplt1,df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "matRRhigh", dis = dis, plt = plt, iqr=iqr,
                          coplt = coplt1,df_t=df_t,
                          pred$matfit+3.29*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matRRlow", dis = dis, plt = plt,iqr=iqr,
                         coplt = coplt1,df_t=df_t,
                         pred$matfit-3.29*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}

plist <- list(dis = c('copd'),
              plt = c("o3h8max"),
              coplt1 = c('temp'),
              coplt2 = c("rh"),
              df_t=7,
              coplt3 = c("o3h8max"))

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)


bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

# print(as_tibble(bbA), n=Inf)
dat3 <- bbA 

FTjust <- function(dis, plt, coplt1, coplt2,df_t, coplt3){
  ### coplt1 = temp
  ### coplt2 = rh
  ### coplt3 = real coplt
  ### df_t = df of time variable t
  # cb.plt = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
  # cb.coplt1 = crossbasis(data[[coplt1]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt2 = crossbasis(data[[coplt2]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  
  # cb.coplt1 = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt2= crossbasis(data[[coplt1]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.plt = crossbasis(data[[coplt2]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  
  plist <- "cb.plt+cb.coplt3+cb.coplt2"
  # plist <- "cb.plt+cb.coplt2"
  
  eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
        s(t,k=",df_t,"*12,fx=T,bs='cr')+

        as.factor(dow)+as.factor(holiday),
 
        family=quasipoisson(link = 'log'), scale=-1,
        control=gam.control(epsilon=0.0000001,maxit=200),
        na.action= na.omit, data=data)", sep='')))
  
  
  iqr<-IQR(data[[coplt2]],na.rm=TRUE)
  eval(parse(text=paste("pred <- crosspred(cb.plt, model,
                        at=iqr, cumul=TRUE)",sep='')))
  matRRfit <- data.frame(fit = "matRRfit", dis = dis, plt = plt, iqr=iqr,
                         coplt = coplt2,df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "matRRhigh", dis = dis, plt = plt, iqr=iqr,
                          coplt = coplt2,df_t=df_t,
                          pred$matfit+3.29*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matRRlow", dis = dis, plt = plt,iqr=iqr,
                         coplt = coplt2,df_t=df_t,
                         pred$matfit-3.29*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}

plist <- list(dis = c('copd'),
              plt = c("o3h8max"),
              coplt1 = c('temp'),
              coplt2 = c("rh"),
              df_t=7,
              coplt3 = c("o3h8max"))

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)


bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

# print(as_tibble(bbA), n=Inf)
dat5 <- bbA 
res1 <- dat1 %>% filter(plt==coplt|coplt=="o3h8max") %>% mutate(Model=ifelse(plt==coplt,"Single model","Ozone-adjusted model"))
res2 <- dat2 %>% mutate(Model="Single model",plt="temp")
res3 <- dat3 %>% mutate(Model="Ozone-adjusted model",plt="temp")  
res4 <- dat4 %>% mutate(Model="Single model",plt="rh")
res5 <- dat5 %>% mutate(Model="Ozone-adjusted model",plt="rh")    
res <- rbind(res1,res2,res3,res4,res5)
save(res,file = "GAM.rda")  
  
  
  
bbA <- rbind(dat1,dat2)
output <- bbA %>% mutate(
  sig=as.factor(ifelse(matRRhigh<=0,1,ifelse(matRRlow>=0,1,0))))

names(output) <- c("dis","plt","iqr","coplt","df_t","lagA","beta","h","l","lag","Significance")

output$lag <- factor(output$lag)
output$Significance <- factor(output$Significance,levels = c(0,1),labels = c("Non-sig","Sig"))


output$Coplt <- output$coplt
output1 <- output %>% mutate(group="CO")

P2_2 <- ggplot(output, aes(lag,beta,ymin = l, ymax = h)) +
  geom_hline(yintercept = 0,linetype='dashed') +
  geom_errorbar(aes(group = plt, col = plt, width=0), 
                size=1, position = position_dodge((width=0.5))) +
  scale_color_jco()+
  geom_point(aes(group = plt, col = plt,shape=Significance),size = 2, 
             position = position_dodge(width=0.5),,fontface = "bold") +
  # scale_color_manual(values = c('#0272B7',"red")) +
  scale_shape_manual(values = c(16,8)) +
  # geom_ribbon(aes(x=lag,ymin =l, ymax = h),
  #             color=NA,fill = "gray10",alpha=0.2) +
  # geom_linerange(aes(x=lag,ymin =l, ymax = h,col=plt2),
  #                size=0.8,color='gray10',alpha=0.5)+
  # facet_grid(Coplt~plt,scales='free',labeller = "label_parsed") +
  facet_grid(~Coplt,scales='free') +
  ylab(expression(paste(beta,"(Risk estimates environment factors → COPD)", sep = "")))+
  xlab("Lag")+
  theme_bw() +
  # mytheme1
  theme(axis.text.x=element_text(size = 10, color= "black"),
        axis.text.y=element_text(size = 10, color= "black"),
        axis.title.x=element_text(size=11),
        axis.title.y=element_text(size=11),
        axis.ticks.length=unit(0.2,'cm'),
        axis.line = element_line(colour = "black"),
        legend.position=c("bottom"),
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.box="vertical",
        legend.margin=margin(-15,unit="pt"),
        # legend.box.margin=margin(10,10,10,10),
        legend.box.spacing = margin(15.5),
        legend.background = element_rect(fill="transparent"),
        strip.background = element_rect(
          color = "white", fill = "white"),
        panel.border = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        panel.grid = element_blank(),
        # panel.border = element_rect(fill=NA,color="black", size=0.9, linetype="solid"),
        plot.title = element_text(size = 18, hjust=0.5))
P2_2
# ggsave("GAM_co.jpg", width = 8, height = 3.5, dpi = 300)

################################## temp###################################

FTjust <- function(dis, plt, coplt1, coplt2,df_t, coplt3){
  ### coplt1 = temp
  ### coplt2 = rh
  ### coplt3 = real coplt
  ### df_t = df of time variable t
  # cb.plt = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
  # cb.coplt1 = crossbasis(data[[coplt1]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt2 = crossbasis(data[[coplt2]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  
  cb.coplt1 = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.plt  = crossbasis(data[[coplt1]], lag=12, argvar=list(fun="lin"), arglag=list(fun="integer"))
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
                        at=24.5, cumul=TRUE,cen=23.5)",sep='')))
  matRRfit <- data.frame(fit = "matRRfit", dis = dis, plt = plt, iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "matRRhigh", dis = dis, plt = plt, iqr=iqr,
                          coplt = coplt3,df_t=df_t,
                          pred$matfit+3.29*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matRRlow", dis = dis, plt = plt,iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit-3.29*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}
#"co","fsp","no2","o3h8max"
plist <- list(dis = c('copd'),
              plt = c('temp'),
              coplt1 = c('temp'),
              coplt2 = c("rh"),
              df_t=c(7),
              coplt3 = c("co","fsp","no2","o3h8max"))

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)


bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

# print(as_tibble(bbA), n=Inf)
dat2 <- bbA

bbA <- dat2
output <- bbA %>% mutate(
  sig=as.factor(ifelse(matRRhigh<=0,1,ifelse(matRRlow>=0,1,0))))

names(output) <- c("dis","plt","iqr","coplt","df_t","lagA","beta","h","l","lag","Significance")

output$lag <- factor(output$lag)
output$Significance <- factor(output$Significance,levels = c(0,1),labels = c("Non-sig","Sig"))


output$Coplt <- output$coplt
# output$plt <- factor(output$plt,levels = c("o3h8max","co","fsp","no2"),labels = c("O3","CO","PM2.5","NO2"))
output$df_t <- factor(output$df_t)

# plt.labs <- c("O3","CO","PM2.5","NO2")
# names(plt.labs) <- c("O3","CO","PM2.5","NO2")
# Coplt.labs <- paste("coplt_",c("O3","CO","PM2.5","NO2"),sep = "")
# names(Coplt.labs) <- c("O3","CO","PM2.5","NO2")
P2_2 <- ggplot(output, aes(lag,beta,ymin = l, ymax = h)) +
  geom_hline(yintercept = 0,linetype='dashed') +
  geom_errorbar(aes(group = plt, col = plt, width=0), 
                size=1, position = position_dodge((width=0.5))) +
  scale_color_jama()+
  geom_point(aes(group = plt, col = plt,shape=Significance),size = 2, 
             position = position_dodge(width=0.5),,fontface = "bold") +
  # scale_color_manual(values = c('#0272B7',"red")) +
  scale_shape_manual(values = c(16,8)) +
  # geom_ribbon(aes(x=lag,ymin =l, ymax = h),
  #             color=NA,fill = "gray10",alpha=0.2) +
  # geom_linerange(aes(x=lag,ymin =l, ymax = h,col=plt2),
  #                size=0.8,color='gray10',alpha=0.5)+
  # facet_grid(Coplt~plt,scales='free',labeller = "label_parsed") +
  facet_grid(~coplt,scales='free',labeller = "label_parsed") +
  ylab(expression(paste(beta,"(Risk estimates air pollutants → COPD)", sep = "")))+
  xlab("Lag")+
  theme_bw() +
  # mytheme1
  theme(axis.text.x=element_text(size = 10, color= "black"),
        axis.text.y=element_text(size = 10, color= "black"),
        axis.title.x=element_text(size=11),
        axis.title.y=element_text(size=11),
        axis.ticks.length=unit(0.2,'cm'),
        axis.line = element_line(colour = "black"),
        legend.position=c("bottom"),
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.box="vertical",
        legend.margin=margin(-15,unit="pt"),
        # legend.box.margin=margin(10,10,10,10),
        legend.box.spacing = margin(15.5),
        legend.background = element_rect(fill="transparent"),
        strip.background = element_rect(
          color = "white", fill = "white"),
        panel.border = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        panel.grid = element_blank(),
        # panel.border = element_rect(fill=NA,color="black", size=0.9, linetype="solid"),
        plot.title = element_text(size = 18, hjust=0.5))
P2_2
# ggsave("COPD_gam_dft_sensitive.jpg", width = 12, height = 8, dpi = 300)






























# 
# 
# 
# FTjust <- function(dis, plt, coplt1, coplt2, l1 = 0,l2 = 3){
#   cb.plt = crossbasis(data[[plt]], lag=3, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
#   cb.coplt1 = crossbasis(data[["temp"]], lag=3, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
#   cb.coplt2 = crossbasis(data[["rh"]], lag=3, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
#   
#   # plist <- ifelse(plt == coplt1, "cb.plt", "cb.plt+cb.coplt1")
#   plist <-  "cb.plt+cb.coplt1+cb.coplt2"
#   # s(rh,k=6+1,fx=T,bs='cr')+
#   
#   eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
#         s(t,k=7*12,fx=T,bs='cr')+
# 
#         as.factor(dow)+as.factor(holiday),
#  
#         family=quasipoisson(link = 'log'), scale=-1,
#         control=gam.control(epsilon=0.0000001,maxit=200),
#         na.action= na.omit, data=data)", sep='')))
#   
#   attr(cb.plt,"lag") <- c(l1,l2)
#   iqr<-IQR(data[[plt]],na.rm=TRUE)
#   eval(parse(text=paste("pred <- crosspred(cb.plt, model,
#                         at=iqr, cumul=TRUE)",sep='')))
#   matRRfit <- data.frame(fit = "matRRfit", dis = dis, plt = plt, iqr=iqr,
#                          coplt = coplt1,
#                          pred$matfit, stringsAsFactors = FALSE)
#   matRRhigh <- data.frame(fit = "matRRhigh", dis = dis, plt = plt, iqr=iqr,
#                           coplt = coplt1,
#                           pred$matfit+2.58*pred$matse, stringsAsFactors = FALSE)
#   matRRlow <- data.frame(fit = "matRRlow", dis = dis, plt = plt,iqr=iqr,
#                          coplt = coplt1,
#                          pred$matfit-2.58*pred$matse, stringsAsFactors = FALSE)
#   out <- bind_rows(matRRfit, matRRhigh, matRRlow)
#   out
# }

################################## temp and rh ###################################

FTjust <- function(dis, plt, coplt1, coplt2,df_t, coplt3){
  ### coplt1 = temp
  ### coplt2 = rh
  ### coplt3 = real coplt
  ### df_t = df of time variable t
  # cb.plt = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
  # cb.coplt1 = crossbasis(data[[coplt1]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt2 = crossbasis(data[[coplt2]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  
  cb.coplt1 = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.plt= crossbasis(data[[coplt1]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt2 = crossbasis(data[[coplt2]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  
  plist <- "cb.plt+cb.coplt1+cb.coplt2"
  # plist <- "cb.plt+cb.coplt2"
  
  eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
        s(t,k=",df_t,"*12,fx=T,bs='cr')+

        as.factor(dow)+as.factor(holiday),
 
        family=quasipoisson(link = 'log'), scale=-1,
        control=gam.control(epsilon=0.0000001,maxit=200),
        na.action= na.omit, data=data)", sep='')))
  
  
  iqr<-IQR(data[[coplt1]],na.rm=TRUE)
  eval(parse(text=paste("pred <- crosspred(cb.plt, model,
                        at=iqr, cumul=TRUE)",sep='')))
  matRRfit <- data.frame(fit = "matRRfit", dis = dis, plt = plt, iqr=iqr,
                         coplt = coplt1,df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "matRRhigh", dis = dis, plt = plt, iqr=iqr,
                          coplt = coplt1,df_t=df_t,
                          pred$matfit+3.29*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matRRlow", dis = dis, plt = plt,iqr=iqr,
                         coplt = coplt1,df_t=df_t,
                         pred$matfit-3.29*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}

plist <- list(dis = c('copd'),
              plt = c("o3h8max"),
              coplt1 = c('temp'),
              coplt2 = c("rh"),
              df_t=7,
              coplt3 = c("o3h8max"))

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)


bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

# print(as_tibble(bbA), n=Inf)
dat3 <- bbA 

FTjust <- function(dis, plt, coplt1, coplt2,df_t, coplt3){
  ### coplt1 = temp
  ### coplt2 = rh
  ### coplt3 = real coplt
  ### df_t = df of time variable t
  # cb.plt = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=2))
  # cb.coplt1 = crossbasis(data[[coplt1]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt2 = crossbasis(data[[coplt2]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  # cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="poly",degree=4))
  
  # cb.coplt1 = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt2= crossbasis(data[[coplt1]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.plt = crossbasis(data[[coplt2]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  
  plist <- "cb.plt+cb.coplt3+cb.coplt2"
  # plist <- "cb.plt+cb.coplt2"
  
  eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
        s(t,k=",df_t,"*12,fx=T,bs='cr')+

        as.factor(dow)+as.factor(holiday),
 
        family=quasipoisson(link = 'log'), scale=-1,
        control=gam.control(epsilon=0.0000001,maxit=200),
        na.action= na.omit, data=data)", sep='')))
  
  
  iqr<-IQR(data[[coplt2]],na.rm=TRUE)
  eval(parse(text=paste("pred <- crosspred(cb.plt, model,
                        at=iqr, cumul=TRUE)",sep='')))
  matRRfit <- data.frame(fit = "matRRfit", dis = dis, plt = plt, iqr=iqr,
                         coplt = coplt2,df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "matRRhigh", dis = dis, plt = plt, iqr=iqr,
                          coplt = coplt2,df_t=df_t,
                          pred$matfit+3.29*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matRRlow", dis = dis, plt = plt,iqr=iqr,
                         coplt = coplt2,df_t=df_t,
                         pred$matfit-3.29*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}

plist <- list(dis = c('copd'),
              plt = c("o3h8max"),
              coplt1 = c('temp'),
              coplt2 = c("rh"),
              df_t=7,
              coplt3 = c("o3h8max"))

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)


bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

# print(as_tibble(bbA), n=Inf)
dat5 <- bbA 
res1 <- dat1 %>% filter(plt==coplt|coplt=="o3h8max") %>% mutate(Model=ifelse(plt==coplt,"Single model","Ozone-adjusted model"))
res2 <- dat2 %>% mutate(Model="Single model",plt="temp")
res3 <- dat3 %>% mutate(Model="Ozone-adjusted model",plt="temp")  
res4 <- dat4 %>% mutate(Model="Single model",plt="rh")
res5 <- dat5 %>% mutate(Model="Ozone-adjusted model",plt="rh")    
res <- rbind(res1,res2,res3,res4,res5)
save(res,file = "GAM.rda")  



bbA <- rbind(dat1,dat2)
output <- bbA %>% mutate(
  sig=as.factor(ifelse(matRRhigh<=0,1,ifelse(matRRlow>=0,1,0))))

names(output) <- c("dis","plt","iqr","coplt","df_t","lagA","beta","h","l","lag","Significance")

output$lag <- factor(output$lag)
output$Significance <- factor(output$Significance,levels = c(0,1),labels = c("Non-sig","Sig"))


output$Coplt <- output$coplt
output1 <- output %>% mutate(group="CO")

P2_2 <- ggplot(output, aes(lag,beta,ymin = l, ymax = h)) +
  geom_hline(yintercept = 0,linetype='dashed') +
  geom_errorbar(aes(group = plt, col = plt, width=0), 
                size=1, position = position_dodge((width=0.5))) +
  scale_color_jco()+
  geom_point(aes(group = plt, col = plt,shape=Significance),size = 2, 
             position = position_dodge(width=0.5),,fontface = "bold") +
  # scale_color_manual(values = c('#0272B7',"red")) +
  scale_shape_manual(values = c(16,8)) +
  # geom_ribbon(aes(x=lag,ymin =l, ymax = h),
  #             color=NA,fill = "gray10",alpha=0.2) +
  # geom_linerange(aes(x=lag,ymin =l, ymax = h,col=plt2),
  #                size=0.8,color='gray10',alpha=0.5)+
  # facet_grid(Coplt~plt,scales='free',labeller = "label_parsed") +
  facet_grid(~Coplt,scales='free') +
  ylab(expression(paste(beta,"(Risk estimates environment factors → COPD)", sep = "")))+
  xlab("Lag")+
  theme_bw() +
  # mytheme1
  theme(axis.text.x=element_text(size = 10, color= "black"),
        axis.text.y=element_text(size = 10, color= "black"),
        axis.title.x=element_text(size=11),
        axis.title.y=element_text(size=11),
        axis.ticks.length=unit(0.2,'cm'),
        axis.line = element_line(colour = "black"),
        legend.position=c("bottom"),
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.box="vertical",
        legend.margin=margin(-15,unit="pt"),
        # legend.box.margin=margin(10,10,10,10),
        legend.box.spacing = margin(15.5),
        legend.background = element_rect(fill="transparent"),
        strip.background = element_rect(
          color = "white", fill = "white"),
        panel.border = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        panel.grid = element_blank(),
        # panel.border = element_rect(fill=NA,color="black", size=0.9, linetype="solid"),
        plot.title = element_text(size = 18, hjust=0.5))
P2_2
# ggsave("GAM_co.jpg", width = 8, height = 3.5, dpi = 300)



