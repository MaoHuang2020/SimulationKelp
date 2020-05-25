rm(list=ls())
setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/Run_Function/Different_varE_nDH/nDH96_Mean")
N.DH<-96

setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/Run_Function/Different_varE_nDH/nDH96_Sd")
N.DH<-96


filenames <- list.files(full.names=TRUE)  
filenames

All <- lapply(filenames,function(i){
  read.csv(i, header=TRUE,row.names=1)
})

  class(All)
  length(All)
  dim(All[[1]])
  All[[1]][1:5,1:6]
  
All.df<-do.call(rbind.data.frame,All) 
  dim(All.df)
  head(All.df)    
  tail(All.df)
  str(All.df)
library(stringr)  
StringSplit<-str_split_fixed(string=All.df$Shorten,"_",7) ### This text string becomes a 336x7 matrix
  head(StringSplit)
  dim(StringSplit)
All.df$nDH<-StringSplit[,4]
All.df$Ne<-StringSplit[,6]
All.df$varE<-StringSplit[,5]


write.csv(All.df,paste(N.DH,"_All.df.csv",sep=""))

All.df %>%
group_by(Year) %>%                         # Specify group indicator (TestSP/SelectSP/Year)
  summarise_at(vars(Mean),              # Specify column
               list(name = mean)) 



library(ggplot2)
library(dplyr)

dev.set()
tiff(file=paste("nDH",N.DH,".tiff",sep=""),width=1400,height=1000,units="px",pointsize=12,res=150)

#par(mfrow=c(2,3))

plot<-ggplot(data=All.df,mapping=aes(x=Cycles,y=Mean))+
  geom_point(aes(shape=Year))+
  geom_line(aes(group=Shorten,color=interaction(TestSP,Year),linetype=SelectSP))+
  theme_bw()+
  labs(x="Year",y="SP genetic mean")+ 
  facet_grid(rows=vars(Ne),cols=vars(varE))+
  ggtitle(paste("number of GPs:",N.DH,sep=""))

plot+scale_color_manual(breaks = c("1000.1yr", "400.1yr", "1000.2yr","400.2yr"),
                        values=c("magenta2","darkturquoise","orangered", "gray17")) 

dev.off()

All.df<-na.omit(All.df)

lmm<-lm(Mean~SelectSP+as.factor(TestSP)+Year+as.factor(Ne)+as.factor(varE),data=All.df)
anova(lmm)
