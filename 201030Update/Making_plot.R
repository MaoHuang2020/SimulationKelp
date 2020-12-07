rm(list=ls())
#setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/Run_Function/Different_varE_nDH/Final_combine_Run_1,2_year/Output_from_bioHPC/Files_Sum")

setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/201030Update/100Cycles/FileSum_GP/Input")

filenames <- list.files(full.names=TRUE)  
filenames

# values<-"Mean"
# nDH<-24

#### Split plots in nDH/nGP

for(values in c("Mean","Sd")){
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
      labs(x="Year",y="GP genetic mean")+ 
      facet_grid(rows=vars(Ne),cols=vars(varE))+
      ylim(ylim)+
      ggtitle(paste("Number of GPs:",nDH,sep=""))+
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
setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/201030Update/100Cycles/FileSum_GP/Plots_ANOVA")
#values<-"Mean"
#nDH<-24

All_2nDH<-NULL
for(values in c("Mean","Sd")){
  for (nDH in c(24,96)){

All.df<-read.csv(paste0("nDH",nDH,"_All.df","_",values,".csv"),sep=",",header=TRUE)


All_2nDH<-rbind(All_2nDH,All.df)
  }
 # ANOVA under each Ne and varE combination, split All_2nDH into Ne and varE combinations
  All_2nDH<-na.omit(All_2nDH)
  All_split<-All_2nDH %>%
    split(list(All_2nDH$Ne,All_2nDH$varE))   
 for (i in 1:length(All_split)){
   for (col in c("SelectSP","TestSP","Year","nDH")){
     All_split[[i]][,col]<-as.factor(All_split[[i]][,col])
   }
   lmm<-lm(Mean~SelectSP+TestSP+Year+nDH+SelectSP*TestSP+SelectSP*Year+SelectSP*nDH+TestSP*Year+TestSP*nDH+Year*nDH,data=All_split[[i]])
   anova_test<-anova(lmm)
   anova_test
   cat("Anova Test\n", file = paste0("CombinedANOVA_",names(All_split)[i],"_",values,".txt"), append = TRUE)
   # This add "Anova Test into the file in the first line
   capture.output(anova_test, file = paste0("CombinedANOVA_",names(All_split)[i],"_",values,".txt"), append = TRUE)
  }   
}  



### calculate the % change across senarios
library(dplyr)

MeanLS<-NULL
for(values in c("Mean","Sd")){
  for (nDH in c(24,96)){
  All.df<-read.csv(paste0("nDH",nDH,"_All.df","_",values,".csv"),sep=",",header=TRUE)
    
   Means<-All.df %>% group_by(SelectSP) %>%                        # Specify group indicator (TestSP/SelectSP/Year)
      summarise_at(vars(Mean),              # Specify column
                list(name = mean))    ###!!! list(meant=mean)
   Means<-as.data.frame(Means)

   MeanSave<-as.data.frame(Means[,2])
   rownames(MeanSave)<-Means[,1]  # 1st value is pheno/400, 2nd value is rand/1000
   colnames(MeanSave)<-paste0(values,"_nDH",nDH)
   MeanLS<-c(MeanLS,MeanSave)
  }
}  

#Mean(nDH24)
# 1.34296
#Mean(nDH96)
# 0.3674465

## SelectSP
#1st value is pheno, 2nd value is rand
print((MeanLS$Mean_nDH96[1]-MeanLS$Mean_nDH96[2])/MeanLS$Mean_nDH96[2])
print((MeanLS$Mean_nDH24[1]-MeanLS$Mean_nDH24[2])/MeanLS$Mean_nDH24[2])

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
  