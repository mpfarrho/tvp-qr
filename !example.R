library(lubridate)
library(GIGrvg)
library(lemon)

load("data_raw.rda")
source("!qrdhs.R")

grid.p <- c(seq(0.05,0.95,by=0.05))

nburn <- 3000
nsave <- 9000
thinfac <- 3

prior <- "dhs"
sl.cn <- "US"
quant_ls <- list()

# ---------------------------------------------------------------------------------------------------------
# insample
Y <- matrix(data_raw[[sl.cn]],ncol=1)
T <- nrow(Y)-1
dates <- as.character(format(date_decimal(as.numeric(time(data_raw[[sl.cn]]))),"%Y-%m-%d")[-1])

# UC-QR
quant_store <- array(NA,dim=c(2,T,length(grid.p)))
dimnames(quant_store) <- list(c("UC-QR","UC-QR-SV"),dates,paste0("p",grid.p*100))

for(sv in c(FALSE,TRUE)){
  message("Estimating: ",sl.cn,", sv=",sv,".")
  Yt <- Y[-1,,drop=FALSE]
  Xt <- Yt^0
  est <- tvpqr.grid(Y=Yt,X=Xt,Xout=NULL,p=grid.p,cpu=6,tvp=prior,sv=sv,fhorz=0,
                    nburn=nburn,nsave=nsave,thinfac=thinfac,out="mcmc")
  bt_store <- est$bt
  sig2_store <- est$sig2
  
  if(sv){
    quant_store["UC-QR-SV",,] <- t(apply(apply(bt_store,c(2,3),remove_outliers),c(2,3),mean,na.rm=TRUE))
  }else{
    quant_store["UC-QR",,] <- t(apply(apply(bt_store,c(2,3),remove_outliers),c(2,3),mean,na.rm=TRUE))
  }
}
quant_ls[[sl.cn]] <- quant_store

# ------------------------------------
# output
library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

shade.col <- "#1974D2"
r80col <- "#2c7bb6"
r90col <- "#d7191c"
shade.col.ser <- c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')
alpha <- 0.9
sl.mods <- c("UCQR-TIS","UCQR-TVS")

Y <- matrix(data_raw[[sl.cn]],ncol=1)[-1,]
quant_tmp <- melt(quant_ls[[sl.cn]]) %>% rename("Model"="Var1","Time"="Var2","Quantile"="Var3") %>% pivot_wider(names_from=Quantile,values_from=value)
levels(quant_tmp$Model) <- c("UCQR-TIS","UCQR-TVS")

quant_long <- melt(quant_tmp) %>% rename("Quantile"="variable")
quant_tmp <- cbind(quant_tmp,Y)

pp <- quant_tmp %>%
  subset(Model %in% sl.mods) %>%
  ggplot(aes(x=as.Date(Time))) +
  
  geom_ribbon(aes(ymin=p5,ymax=p95),alpha=alpha,fill=shade.col.ser[1]) + 
  geom_ribbon(aes(ymin=p10,ymax=p90),alpha=alpha,fill=shade.col.ser[2]) + 
  geom_ribbon(aes(ymin=p15,ymax=p85),alpha=alpha,fill=shade.col.ser[3]) + 
  geom_ribbon(aes(ymin=p20,ymax=p80),alpha=alpha,fill=shade.col.ser[4]) + 
  geom_ribbon(aes(ymin=p25,ymax=p75),alpha=alpha,fill=shade.col.ser[5]) + 
  geom_ribbon(aes(ymin=p30,ymax=p70),alpha=alpha,fill=shade.col.ser[6]) + 
  geom_ribbon(aes(ymin=p35,ymax=p65),alpha=alpha,fill=shade.col.ser[7]) + 
  geom_ribbon(aes(ymin=p40,ymax=p60),alpha=alpha,fill=shade.col.ser[8]) + 
  geom_ribbon(aes(ymin=p45,ymax=p55),alpha=alpha,fill=shade.col.ser[9]) + 
  geom_line(aes(y=value,group=Quantile),color="black",size=0.7,data=subset(subset(quant_long,Model %in% sl.mods),Quantile %in% c("p5","p95"))) +
  
  geom_hline(yintercept=0,size=0.5,color="red") +
  geom_line(aes(y=p50),color="white",size=0.7) +
  facet_rep_wrap(.~Model,scales="fixed",repeat.tick.labels = 'all') + xlab("") + ylab(expression(pi[t])) +
  coord_cartesian(expand=FALSE) +
  theme_cowplot() + 
  theme(axis.line.x = element_blank(), 
        panel.grid.major=element_line(color="grey80",linetype="solid",size=0.5), panel.grid.minor=element_line(color="grey80",linetype="dashed",size=0.5),
        strip.background = element_blank(), strip.text = element_text(face="bold",size=17),
        panel.spacing.x = unit(1,"cm"),
        axis.text = element_text(size=16), axis.title = element_text(size=17),
        plot.margin = unit(c(0,1,0,0), "cm"))

# pdf(paste0("!latex_R1/plots/insample/",sl.cn,"_",prior,".pdf"),width=16,height=3)
print(pp)
# dev.off()
