### I,  making plots

rm(list=ls())

WD<-"/results_folder/"
setwd(paste0(WD,"output"))
plotdatafdr<-paste0(WD,"output/plottingdata/")

filenames <- list.files(full.names=TRUE)  
filenames

#### Split plots in nDH/nGP

for(values in c("Mean","Sd","SPselInt","GPselInt","GPcor")){
  
  if (values=="Mean"){
    ylim<-c(-0.5,10)
    ylab<-"GP Genetic Mean"
  }else if (values=="Sd"){
    ylim<-c(0.45,1.1)
    ylab<-"GP Genetic Variance"
  }else if (values=="GPcor"){
    ylim<-c(0.1,0.8)
    ylab<-"GP GS cor"
  }else if (values=="GPselInt"){
    ylim<-c(0,1.8)
    ylab<-"GP Selection Intensity"
  }else if (values=="SPselInt"){
    ylim<-c(-0.5,2.0)
    ylab<-"SP Selection Intensity"
  }
  
  Plotfnc(values=values,ylim=ylim,ylab=ylab)
} 
  

#### Plotting function for each "values"

Plotfnc<-function(values,ylim,ylab){  
  for (nDH in c(24,96)){ #different nGP/nDH
  SumFile<-filenames[grepl(paste0("nDH",nDH),filenames)]
  
  MeanFile<-SumFile[grep(paste0(values,".csv"),SumFile)]
  
  All <- lapply(MeanFile,function(i){
    read.csv(i, header=TRUE,row.names=1)
  })

  All.df<-do.call(rbind.data.frame,All) 

  ### RM cycle1, RM other cycles for 2yr scenario
  All.df2<-All.df
  All.df<-All.df[!All.df$Cycles==1,]
  All.df<-All.df[!(All.df$Year=="2yr"& All.df$Cycles%in%c(3,5,7,9)),]
  
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
  
  write.csv(All.df,paste0(plotdatafdr,"nDH",nDH,"_All.df","_",values,".csv"))
  
  library(ggplot2)
  library(dplyr)
  
  #All.df$NumPlots<-as.factor(All.df$TestSP) # changed the names TestSP to "NumPlots"
  All.df$NumCross<-as.factor(All.df$TestSP) #Changed NumPlots to NumCross per reveiwer request
  All.df$CycleTime<-as.factor(All.df$Year)  # changed the names "Year" to "CycleTime"
  #dev.set(i)
  tiff(file=paste0(plotdatafdr,"Number_of_GPs_nDH",nDH,"_",values,".tiff"),width=1400,height=1000,units="px",pointsize=12,res=150)
    
    #par(mfrow=c(2,3))
    plot<-ggplot(data=All.df,mapping=aes(x=Cycles,y=Mean))+
      geom_point(aes(shape=NumCross))+
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
  }
}



#### II, Conducting ANOVA

rm(list=ls())
setwd("/results_folder/output/plottingdata")


for(values in c("Mean","Sd","GPcor")){
  
All_2nDH<-NULL
  for (nDH in c(24,96)){
All.df<-read.csv(paste0("nDH",nDH,"_All.df","_",values,".csv"),sep=",",header=TRUE)
All_2nDH<-rbind(All_2nDH,All.df)
     }
 # ANOVA under each Ne and varE combination, split All_2nDH into Ne and varE combinations
  All_2nDH<-na.omit(All_2nDH)
All_2nDH$NumCross<-All_2nDH$TestSP    ## Rename NumCross to replace TestSP
All_2nDH$CycleTime<-All_2nDH$Year     ## Rename CycleTime to replace Year

  All_split<-All_2nDH %>%
    split(list(All_2nDH$Ne,All_2nDH$varE))   
  
  # All_split$NumCross<-All_split$TestSP
  
 for (i in 1:length(All_split)){
   for (col in c("SelectSP","NumCross","CycleTime","nGP")){
     All_split[[i]][,col]<-as.factor(All_split[[i]][,col])
   }
   lmm<-lm(Mean~SelectSP+NumCross+CycleTime+nGP+SelectSP*NumCross+SelectSP*CycleTime+SelectSP*nGP+NumCross*CycleTime+NumCross*nGP+CycleTime*nGP,data=All_split[[i]])
   anova_test<-anova(lmm)
   anova_test
   cat("Anova Test\n", file = paste0("CombinedANOVA_",names(All_split)[i],"_",values,".txt"), append = TRUE)
   # This add "Anova Test into the file in the first line
   capture.output(anova_test, file = paste0("CombinedANOVA_",names(All_split)[i],"_",values,".txt"), append = TRUE)
  }   
}  


### III, calculate the % change across senarios
library(dplyr)

MeanLS<-NULL
for(values in c("Mean","Sd")){
  for (nDH in c(24,96)){
  All.df<-read.csv(paste0(plotdatafdr,"nDH",nDH,"_All.df","_",values,".csv"),sep=",",header=TRUE)
    
   Means<-All.df %>% group_by(TestSP) %>%                        # Specify group indicator (TestSP/SelectSP/Year)
      summarise_at(vars(Mean),              # Specify column
                list(name = mean))    ###!!! list(meant=mean)
   Means<-as.data.frame(Means)

   MeanSave<-as.data.frame(Means[,2])
   rownames(MeanSave)<-Means[,1]  # 1st value is pheno/400, 2nd value is rand/1000
   colnames(MeanSave)<-paste0(values,"_nDH",nDH)
   MeanLS<-c(MeanLS,MeanSave)
  }
}  


## SelectSP
#1st value is pheno, 2nd value is rand
print((MeanLS$Mean_nDH96[1]-MeanLS$Mean_nDH96[2])/MeanLS$Mean_nDH96[2])
print((MeanLS$Mean_nDH24[1]-MeanLS$Mean_nDH24[2])/MeanLS$Mean_nDH24[2])

## TestSP/NumCross
#1st value is 400, 2nd value is 1000
print((MeanLS$Mean_nDH96[2]-MeanLS$Mean_nDH96[1])/MeanLS$Mean_nDH96[1])
print((MeanLS$Mean_nDH24[2]-MeanLS$Mean_nDH24[1])/MeanLS$Mean_nDH24[1])
## Year
# 1st value is 1yr, 2nd value is 2yr
print((MeanLS$Mean_nDH96[1]-MeanLS$Mean_nDH96[2])/MeanLS$Mean_nDH96[2])
print((MeanLS$Mean_nDH24[1]-MeanLS$Mean_nDH24[2])/MeanLS$Mean_nDH24[2])



MeanLS<-NULL
for(values in c("GPcor")){
  for (nDH in c(24,96)){
    All.df<-read.csv(paste0(plotdatafdr,"nDH",nDH,"_All.df","_",values,".csv"),sep=",",header=TRUE)
    
    Means<-All.df %>% group_by(Year) %>%                        # Specify group indicator (TestSP/SelectSP/Year)
      summarise_at(vars(Mean),              # Specify column
                   list(name = mean))    ###!!! list(meant=mean)
    Means<-as.data.frame(Means)
    
    MeanSave<-as.data.frame(Means[,2])
    rownames(MeanSave)<-Means[,1]  # 1st value is pheno/400, 2nd value is rand/1000
    colnames(MeanSave)<-paste0(values,"_nDH",nDH)
    MeanLS<-c(MeanLS,MeanSave)
  }
}  


## 1st one is when nDH=96
MeanLS[1]
## 2nd ones is when nDH=24
MeanLS[2]





  
