rm(list=ls())
#setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/Run_Function/Different_varE_nDH/Final_combine_Run_1,2_year/Output_from_bioHPC/Files_Sum")

setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/201030Update/100Cycles/FileSum_GP/Input")

filenames <- list.files(full.names=TRUE)  
filenames

# values<-"Mean"
# nDH<-24

#### Split plots in nDH/nGP
#### Figure 2 (Mean) and 3 (Sd)

for(values in c("Mean","Sd")){
  
  ylab<-ifelse(values=="Mean","GP Genetic Mean","GP Genetic Variance")
  
  for (nDH in c(24,96)){ #different nGP/nDH
  SumFile<-filenames[grepl(paste0("nDH",nDH),filenames)]
  
  MeanFile<-SumFile[grep(paste0(values,".csv"),SumFile)]
  
  All <- lapply(MeanFile,function(i){
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
  All.df$nGP<-StringSplit[,4]  # Add this nGP level
  
  All.df$Ne<-StringSplit[,6]
  All.df$varE<-StringSplit[,5]
  
  # Change Label the varE to h2=0.2 when varE4 and h2=0.5 when varE1
  All.df$varE[All.df$varE=="varE1"]<-"heritability: 0.5"
  All.df$varE[All.df$varE=="varE4"]<-"heritability: 0.2"
    
  # Change Label the Ne
  All.df$Ne[All.df$Ne=="Ne60"]<-"Ne=60"
  All.df$Ne[All.df$Ne=="Ne600"]<-"Ne=600"
  
  
  dim(All.df)
  head(All.df)
  
  #write.csv(All.df,paste0("nDH",nDH,"_All.df","_",values,".csv"))
  
  library(ggplot2)
  library(dplyr)
  
  All.df$NumPlots<-as.factor(All.df$TestSP) # changed the names TestSP to "NumPlots"
  All.df$CycleTime<-as.factor(All.df$Year)  # changed the names "Year" to "CycleTime"
  #dev.set(i)
    
    tiff(file=paste("Number_of_GPs_nDH",nDH,"_",values,".tiff",sep=""),width=1400,height=1000,units="px",pointsize=12,res=150)
    
    #par(mfrow=c(2,3))
    ##### !!!!! Do y=Mean or Y=Sd
    
    if (values=="Mean"){
      ylim<-c(-0.5,7)
    }else{
      ylim<-c(0.5,1.1)
    }
    
    
    plot<-ggplot(data=All.df,mapping=aes(x=Cycles,y=Mean))+
      geom_point(aes(shape=NumPlots))+
      geom_line(aes(group=Shorten,color=CycleTime,linetype=SelectSP)) +
      theme_bw()+
      theme(panel.grid.major.x=element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black"))+
      labs(x="Year",y=ylab)+ 
      facet_grid(rows=vars(Ne),cols=vars(varE))+
      ylim(ylim)+
      ggtitle(paste("Number of GPs:",nDH,sep=""))+
      scale_color_manual(breaks = c("1yr","2yr"),values=c("orangered", "gray17"))+
      scale_shape_manual(values=c(3,16))
    print(plot)   ### Need to print(plot) otherwise cannot save out as tiff
    dev.off()
    
    #### Add stderr bar  +
    #geom_errorbar(aes(ymin=Mean, ymax=Mean+StdErr), width=.2,position=position_dodge(.9))
    #strip.background = element_rect(colour="white", fill="white"),
    #
  }
}

## how to add 2 newlines
## cat("\n\n", file = "tests.txt", append = TRUE)




######
###### Split Plots in NumPlots !!!!!!!!!!
for(values in c("Mean","Sd")){
  #for (nDH in c(24,96)){ #different nGP/nDH
  #  SumFile<-filenames[grepl(paste0("nDH",nDH),filenames)]
    SumFile<-filenames
    MeanFile<-SumFile[grep(paste0(values,".csv"),SumFile)]
    
    All <- lapply(MeanFile,function(i){
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
    All.df$nGP<-StringSplit[,4]  # Add this nGP level
    
    All.df$Ne<-StringSplit[,6]
    All.df$varE<-StringSplit[,5]
    
    dim(All.df)
    head(All.df)

    #write.csv(All.df,paste0("nDH",nDH,"_All.df","_",values,".csv"))
    
    library(ggplot2)
    library(dplyr)
    
    All.df$NumPlots<-as.factor(All.df$TestSP) # changed the names TestSP to "NumPlots"
    All.df$CycleTime<-as.factor(All.df$Year)  # changed the names "Year" to "CycleTime"
  #dev.set(i)
    All.df0<-All.df
  for (NumPlots in c(400, 1000)){
 
    tiff(file=paste("NumPlots",NumPlots,"_",values,".tiff",sep=""),width=1400,height=1000,units="px",pointsize=12,res=150)
    
    #par(mfrow=c(2,3))
    ##### !!!!! Do y=Mean or Y=Sd
    
      if (values=="Mean"){
  ylim<-c(-0.5,7)
  }else{
  ylim<-c(0.5,1.1)
  }

    All.df<-droplevels(All.df0[which(All.df0$NumPlots==NumPlots),])
    
    plot<-ggplot(data=All.df,mapping=aes(x=Cycles,y=Mean))+
      geom_point(aes(shape=nGP))+
      geom_line(aes(group=Shorten,color=CycleTime,linetype=SelectSP)) +
      theme_bw()+
      labs(x="Year",y="GP genetic mean")+ 
      facet_grid(rows=vars(Ne),cols=vars(varE))+
      ylim(ylim)+
      ggtitle(paste("number of Plots:",NumPlots,sep=""))+
      scale_color_manual(breaks = c("1yr","2yr"),values=c("orangered", "gray17"))+
      scale_shape_manual(values=c(3,16))
    

    print(plot)   ### Need to print(plot) otherwise cannot save out as tiff
    dev.off()
 
 #### Add stderr bar   
   #geom_errorbar(aes(ymin=Mean, ymax=Mean+StdErr), width=.2,position=position_dodge(.9))+ 
  }
}

## how to add 2 newlines
## cat("\n\n", file = "tests.txt", append = TRUE)





#### Conducting ANOVA
#### Table 1 and Supplemental Table 1
rm(list=ls())
setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/201030Update/100Cycles/FileSum_GP/Plots_ANOVA")
#values<-"Mean"
#nDH<-24
library(plyr)

for(values in c("Mean","Sd")){
  
  All_2nDH<-NULL  ## cannot be outside the loop !!
  for (nDH in c(24,96)){

All.df<-read.csv(paste0("nDH",nDH,"_All.df","_",values,".csv"),sep=",",header=TRUE)


All_2nDH<-rbind(All_2nDH,All.df)
  }
 # ANOVA under each Ne and varE combination, split All_2nDH into Ne and varE combinations
  All_2nDH<-na.omit(All_2nDH)
  All_split<-split(x=All_2nDH,f=list(All_2nDH$Ne,All_2nDH$varE))   
  
 for (i in 1:length(All_split)){
   for (col in c("SelectSP","TestSP","Year","nDH")){
     All_split[[i]][,col]<-as.factor(All_split[[i]][,col])
     
     print(t.test(All_split[[i]]$Mean[All_split[[i]]$TestSP==1000],All_split[[i]]$Mean[All_split[[i]]$TestSP==400]),paired=FALSE)
     
   }
   print(unique(paste0(All_split[[i]]$Ne,"_",All_split[[i]]$varE)))
   lmm<-lm(Mean~SelectSP+TestSP+Year+nDH+SelectSP*TestSP+SelectSP*Year+SelectSP*nDH+TestSP*Year+TestSP*nDH+Year*nDH,data=All_split[[i]])
   anova_test<-anova(lmm)
   print(anova_test)
   cat("Anova Test\n", file = paste0("CombinedANOVA_",names(All_split)[i],"_",values,".txt"), append = TRUE)
   ## This add "Anova Test into the file in the first line
   capture.output(anova_test, file = paste0("CombinedANOVA_",names(All_split)[i],"_",values,".txt"), append = TRUE)
  }   
}  



### calculate the % change across senarios
library(dplyr)

### Changed this to average across nDH, but just each scenario (Dec 22nd, 2020)
### 
MeanLS<-NULL
for(values in c("Mean","Sd")){
  All.df<-NULL
  for (nDH in c(24,96)){
  tmp<-read.csv(paste0("nDH",nDH,"_All.df","_",values,".csv"),sep=",",header=TRUE)
  All.df<-rbind(All.df,tmp)  
  }
   Means<-All.df %>% group_by(nDH) %>%                        # Specify group indicator (SelectSP/Year/TestSP)
      summarise_at(vars(Mean),              # Specify column
                list(name = mean))    ###!!! list(meant=mean)
   Means<-as.data.frame(Means)

  # MeanSave<-as.data.frame(Means[,2])
  # rownames(MeanSave)<-Means[,1]  # 1st value is pheno/400, 2nd value is rand/1000
  # colnames(MeanSave)<-paste0(values) #"_nDH",nDH
  # MeanLS<-c(MeanLS,MeanSave)
  #}
  
}  




### OR do
aggregate(All.df[,2],list(All.df$SelectSP),mean)

print((MeanLS$Mean[1]-MeanLS$Mean[2])/MeanLS$Mean[2])
#SelectSP (pheno-rand)/rand    #1.014842
#Year (Means[Means$Year=="1yr",2]-Means[Means$Year=="2yr",2])/Means[Means$Year=="2yr",2]  #0.4540391
#TestSP (Means[Means$TestSP==1000,2]-Means[Means$TestSP==400,2])/Means[Means$TestSP==400,2] #0.1054314


#Mean(nDH24)
# 1.34296
#Mean(nDH96)
# 0.3674465

## SelectSP
#1st value is pheno, 2nd value is rand
print((MeanLS$Mean_nDH24[1]-MeanLS$Mean_nDH24[2])/MeanLS$Mean_nDH24[2])
print((MeanLS$Mean_nDH96[1]-MeanLS$Mean_nDH96[2])/MeanLS$Mean_nDH96[2])

## SelectSP
#1st value is 400, 2nd value is 1000
print((MeanLS$Mean_nDH96[2]-MeanLS$Mean_nDH96[1])/MeanLS$Mean_nDH96[1])
print((MeanLS$Mean_nDH24[2]-MeanLS$Mean_nDH24[1])/MeanLS$Mean_nDH24[1])
## Year
# 1st value is 1yr, 2nd value is 2yr
print((MeanLS$Mean_nDH96[1]-MeanLS$Mean_nDH96[2])/MeanLS$Mean_nDH96[2])
print((MeanLS$Mean_nDH24[1]-MeanLS$Mean_nDH24[2])/MeanLS$Mean_nDH24[2])


## SelectSP, 
#nDH96  0.72
#nDH24  1.55
## TestSP
#nDH96 0.13
#nDH24 0.08
## Year
#nDH96 0.46
#nDH24 0.46

# MeanLS
## 1st value is pheno, 2nd value is rand
#$Mean_nDH24
#[1] 1.9293656 0.7565552

#$Mean_nDH96
#[1] 2.320251 1.352602


#### Adding in plot, x-axis on the 5-generation end Gain from each scheme; y-axis on the 5-generation end Genetic Variance from each scheme
#### Figure 4

values<-"Mean"
tmp<-NULL
names<-NULL
  for (nDH in c(24,96)){
    All.df<-read.csv(paste0("nDH",nDH,"_All.df","_",values,".csv"),sep=",",header=TRUE)
    
    tmp<-c(tmp,All.df$Mean[All.df$Cycles==7])
    names<-c(names,as.character(All.df$Shorten[All.df$Cycles==7]))
  }
X<-tmp
names(X)<-names

values<-"Sd"
tmp<-NULL
names<-NULL
for (nDH in c(24,96)){
  All.df<-read.csv(paste0("nDH",nDH,"_All.df","_",values,".csv"),sep=",",header=TRUE)
  
  tmp<-c(tmp,All.df$Mean[All.df$Cycles==7])
  names<-c(names,as.character(All.df$Shorten[All.df$Cycles==7]))
}
Y<-tmp
names(Y)<-names

library(stringr)
identical(str_split_fixed(names(X),"_",7)[1:6],str_split_fixed(names(Y),"_",7)[1:6]) # TRUE=same treatment order
dataf<-data.frame(X,Y)

tiff(file=paste("Figure4.tiff",sep=""),width=1400,height=1000,units="px",pointsize=12,res=150)

library(ggplot2)
plot<-ggplot(data=dataf,mapping=aes(x=X,y=Y))+
  geom_point()+
  labs(x="Genetic gain",y="Genetic variance")+
  theme_bw()

print(plot) 
dev.off()

head(dataf)
colnames(dataf)[1]<-"Mean"
colnames(dataf)[2]<-"Sd"
dataf$SelectSP<-str_split_fixed(rownames(dataf),"_",7)[,1]
dataf$NumPlots<-str_split_fixed(rownames(dataf),"_",7)[,2]
dataf$Year<-str_split_fixed(rownames(dataf),"_",7)[,3]
dataf$nGP<-str_split_fixed(rownames(dataf),"_",7)[,4]
dataf$varE<-str_split_fixed(rownames(dataf),"_",7)[,5]
dataf$Ne<-str_split_fixed(rownames(dataf),"_",7)[,6]

# install.packages("reshape2")
# library("reshape2")
# dataf_long<-melt(dataf,id.vars=c("varE","Ne"))
write.csv(dataf_long,"dataf_FinalGain_Variance.csv")

#### Figure 4 JL version Dec 22, 2020
library(tidyverse)
gainSD <- read.csv("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/201030Update/100Cycles/FileSum_GP/Plots_ANOVA/dataf_FinalGain_Variance.csv",sep=",",header=T)

colnames(gainSD)[1] <- "Name"
gainSD <- gainSD[order(gainSD$varE, gainSD$Ne),]
plot(gainSD$Mean, gainSD$Sd, pch=16)

pdf("gainSDbyInterventionLines.pdf", height=8, width=8)
op <- par(mfrow=c(2,2))
for (f in 4:7){
  sb <- (4:7)[-(f-3)]
  print(c(sb, f))
  gainSD <- gainSD[order(gainSD[,8], gainSD[,9], gainSD[,sb[1]], gainSD[,sb[2]], gainSD[,sb[3]], gainSD[,f]),]
  print(gainSD)
  curNxt <- c(2, 1, 2, 1)[f-3]
  mainTitle <- c("Selection on SP", "Number of Plots Eval.", "Cycle Time in Years", "nGP per Parental SP")[f-3]
  plot(gainSD$Mean, gainSD$Sd, pch=16, col=rep(curNxt:(3-curNxt), 32), xlab="Final Genetic Mean", ylab="Final Genetic Variance", main=mainTitle)
  for (i in 0:31){
    lines(gainSD$Mean[i*2 + 1:2], gainSD$Sd[i*2 + 1:2], lwd=0.5, col="gray")
  }
    #text(6.6, 0.93, labels=c("A", "B", "C", "D")[f-3], cex=1.5)
}

dev.off()


##########
rm(list=ls())
setwd("/Users/maohuang/Desktop/Kelp/SNP_calling/Blast_SNPs")
GO_list<-read.csv("nearest.gene_FLK_GWAS_Output_withGO.csv",sep=",",header=T)
  dim(GO_list)
  head(GO_list)  
library(MASS)  

GO<-GO_list$GoTerm
GO.freq<-table(GO)
pie(GO.freq)
  