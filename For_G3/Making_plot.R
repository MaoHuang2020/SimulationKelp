### STEP 2,making plots

rm(list=ls())
#setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/201030Update/100Cycles/FileSum_GP/Input")

WD<-"/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/210602Update/20Cycles/"
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



#### Conducting ANOVA
#setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/201030Update/100Cycles/FileSum_GP/Plots_ANOVA")
#values<-"Mean"
#nDH<-24
rm(list=ls())
setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/210602Update/100Cycles/output/plottingdata")


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



### calculate the % change across senarios
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

### Mean
## SelectSP, % change from random to pheno
#nDH96   0.5277912
#nDH24   1.177964
## TestSP  % change from 400 to 1000
#nDH96  0.1691327
#nDH24  0.07015513
## Year   % change from 2yr to 1yr
#nDH96  0.698914
#nDH24  0.7332221

# MeanLS
## 1st value is pheno, 2nd value is rand
#MeanLS$Mean_nDH24   $Sd_nDH24
#[1] 2.178771        0.8821729

#MeanLS$Mean_nDH96   $Sd_nDH96
#[1] 2.949128        0.790105

#change from nDH24 to nDH96 
# 0.3535741*100%


########## Sd output
#$Sd_nDH24
# 0.8821729
#$Sd_nDH96
# 0.790105

### TestSP  400->1000
# $Sd_nDH24 (1000-400)/400= 0.04034996
# [1] 0.8647270 0.8996187

# $Sd_nDH96 (1000-400)/400= 0.08701718
# [1] 0.7571620 0.8230481

### SelectSP  pheno-> rand
# $Sd_nDH24 (pheno-rand)/rand= -0.06925011
# [1] 0.8505320 0.9138137

# $Sd_nDH96  (pheno-rand)/rand= -0.05844567
# [1] 0.7663209 0.8138892

### Year  1yr -> 2yr
# $Sd_nDH24  (1yr-2yr)/2yr= -0.04712685
# [1] 0.8668611 0.9097340

# $Sd_nDH96 (1yr-2yr)/2yr= -0.108452
# [1] 0.7572085 0.8493188


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

####### GPcor
#nGP
#nDH96   0.454859
#nDH24   0.4466945

#TestSP 400 -> 1000 (increased, not always sig)
#nDH24 0.4283803 0.4650087
#nDH96 0.4297344 0.4799837

#SelectSP  pheno -> rand (decreased, always sig)
#nDH24 0.4162794 0.4771097
#nDH96 0.4292512 0.4804669

#Year  1yr -> 2yrs (decreased, not always sig)
#nDH24 0.4522131 0.4370370
#nDH96 0.4591701 0.4473146






#### Add stderr bar   
#geom_errorbar(aes(ymin=Mean, ymax=Mean+StdErr), width=.2,position=position_dodge(.9))+ 

## how to add 2 newlines
## cat("\n\n", file = "tests.txt", append = TRUE)




# ###### DID not pick this approach in the end
# ###### Split Plots in NumPlots !!!!!!!!!!
# for(values in c("Mean","Sd")){
#   #for (nDH in c(24,96)){ #different nGP/nDH
#   #  SumFile<-filenames[grepl(paste0("nDH",nDH),filenames)]
#     SumFile<-filenames
#     MeanFile<-SumFile[grep(paste0(values,".csv"),SumFile)]
#     
#     All <- lapply(MeanFile,function(i){
#       read.csv(i, header=TRUE,row.names=1)
#     })
#     
#     class(All)
#     length(All)
#     dim(All[[1]])
#     All[[1]][1:5,1:6]
#     
#     All.df<-do.call(rbind.data.frame,All) 
#     dim(All.df)
#     head(All.df)    
#     tail(All.df)
#     str(All.df)
#     
#     library(stringr)  
#     
#     StringSplit<-str_split_fixed(string=All.df$Shorten,"_",7) ### This text string becomes a 336x7 matrix
#     
#     head(StringSplit)
#     dim(StringSplit)
#     All.df$nGP<-StringSplit[,4]  # Add this nGP level
#     
#     All.df$Ne<-StringSplit[,6]
#     All.df$varE<-StringSplit[,5]
#     
#     dim(All.df)
#     head(All.df)
# 
#     #write.csv(All.df,paste0("nDH",nDH,"_All.df","_",values,".csv"))
#     
#     library(ggplot2)
#     library(dplyr)
#     
#     All.df$NumPlots<-as.factor(All.df$TestSP) # changed the names TestSP to "NumPlots"
#     All.df$CycleTime<-as.factor(All.df$Year)  # changed the names "Year" to "CycleTime"
#   #dev.set(i)
#     All.df0<-All.df
#   for (NumPlots in c(400, 1000)){
#  
#     tiff(file=paste("NumPlots",NumPlots,"_",values,".tiff",sep=""),width=1400,height=1000,units="px",pointsize=12,res=150)
#     
#     #par(mfrow=c(2,3))
#     ##### !!!!! Do y=Mean or Y=Sd
#     
#       if (values=="Mean"){
#   ylim<-c(-0.5,7)
#   }else{
#   ylim<-c(0.5,1.1)
#   }
# 
#     All.df<-droplevels(All.df0[which(All.df0$NumPlots==NumPlots),])
#     
#     plot<-ggplot(data=All.df,mapping=aes(x=Cycles,y=Mean))+
#       geom_point(aes(shape=nGP))+
#       geom_line(aes(group=Shorten,color=CycleTime,linetype=SelectSP)) +
#       theme_bw()+
#       theme(panel.grid.major.x=element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black"))+
#       labs(x="Year",y="GP genetic mean")+ 
#       facet_grid(rows=vars(Ne),cols=vars(varE))+
#       ylim(ylim)+
#       ggtitle(paste("number of Plots:",NumPlots,sep=""))+
#       scale_color_manual(breaks = c("1yr","2yr"),values=c("orangered", "gray17"))+
#       scale_shape_manual(values=c(3,16))
#   
#     print(plot)   ### Need to print(plot) otherwise cannot save out as tiff
#     dev.off()
#  
#  #### Add stderr bar   
#    #geom_errorbar(aes(ymin=Mean, ymax=Mean+StdErr), width=.2,position=position_dodge(.9))+ 
#   }
# }
# 
# ## how to add 2 newlines
# ## cat("\n\n", file = "tests.txt", append = TRUE)
# 



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
  