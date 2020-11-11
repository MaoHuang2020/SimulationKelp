rm(list=ls())
#setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/Run_Function/Different_varE_nDH/Final_combine_Run_1,2_year/Output_from_bioHPC/Files_Sum")

setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/201030Update/100Cycles/FileSum_GP")

filenames <- list.files(full.names=TRUE)  
filenames

# values<-"Mean"
# nDH<-24


for(values in c("Mean","Sd")){
  for (nDH in c(24,96)){
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
    All.df$nDH<-StringSplit[,4]
    
    All.df$Ne<-StringSplit[,6]
    All.df$varE<-StringSplit[,5]
    
    dim(All.df)
    head(All.df)
    
    #write.csv(All.df,paste0("nDH",nDH,"_All.df","_",values,".csv"))
    
    library(ggplot2)
    library(dplyr)
    
    All.df$TestSP<-as.factor(All.df$TestSP)
  #dev.set(i)
    tiff(file=paste("nDH",nDH,"_",values,".tiff",sep=""),width=1400,height=1000,units="px",pointsize=12,res=150)
    
    #par(mfrow=c(2,3))
    ##### !!!!! Do y=Mean or Y=Sd
    
      if (values=="Mean"){
  ylim<-c(-0.5,8)
  }else{
  ylim<-c(0.5,1.1)
  }

    plot<-ggplot(data=All.df,mapping=aes(x=Cycles,y=Mean))+
      geom_point(aes(shape=TestSP))+
      ylim(ylim)+
      geom_line(aes(group=Shorten,color=Year,linetype=SelectSP))+
      theme_bw()+
      labs(x="Year",y="GP genetic mean")+ 
      facet_grid(rows=vars(Ne),cols=vars(varE))+
      ggtitle(paste("number of GPs:",nDH,sep=""))+
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

# All.df %>%
#   group_by(SelectSP) %>%                        # Specify group indicator (TestSP/SelectSP/Year)
# summarise_at(vars(Mean),              # Specify column
#                list(name = mean))    ###!!! list(meant=mean)

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
  