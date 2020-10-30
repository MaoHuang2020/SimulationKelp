rm(list=ls())
setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/Run_Function/Different_varE_nDH/Final_combine_Run_1,2_year/Output_from_bioHPC/Files_Sum")


filenames <- list.files(full.names=TRUE)  
filenames

# values<-"Mean"
# nDH<-25


anova_list<-NULL

for(values in c("Mean","Sd")){
  for (nDH in c(25,96)){
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
    
    write.csv(All.df,paste0("nDH",nDH,"_All.df","_",values,".csv"))
    
    library(ggplot2)
    library(dplyr)
    
  #dev.set(i)
    tiff(file=paste("nDH",nDH,"_",values,".tiff",sep=""),width=1400,height=1000,units="px",pointsize=12,res=150)
    
    #par(mfrow=c(2,3))
    ##### !!!!! Do y=Mean or Y=Sd
    
    plot<-ggplot(data=All.df,mapping=aes(x=Cycles,y=Mean))+
      geom_point(aes(shape=Year))+
      geom_line(aes(group=Shorten,color=interaction(TestSP,Year),linetype=SelectSP))+
      theme_bw()+
      labs(x="Year",y="SP genetic mean")+ 
      facet_grid(rows=vars(Ne),cols=vars(varE))+
      ggtitle(paste("number of GPs:",nDH,sep=""))+
      scale_color_manual(breaks = c("1000.1yr", "400.1yr", "1000.2yr","400.2yr"),
                            values=c("magenta2","darkturquoise","orangered", "gray17")) 
    print(plot)   ### Need to print(plot) otherwise cannot save out as tiff
    dev.off()
    
    
    All.df<-na.omit(All.df)
    
    lmm<-lm(Mean~SelectSP+as.factor(TestSP)+Year+as.factor(Ne)+as.factor(varE),data=All.df)
    anova_test<-anova(lmm)
    anova_test
    cat("Anova Test\n", file = paste0("ANOVA_nDH",nDH,"_",values,".txt"), append = TRUE)
    # This add "Anova Test into the file in the first line
    capture.output(anova_test, file = paste0("ANOVA_nDH",nDH,"_",values,".txt"), append = TRUE)
  }
}

## how to add 2 newlines
## cat("\n\n", file = "tests.txt", append = TRUE)


####
values<-"Mean"
nDH<-25

for(values in c("Mean","Sd")){
  for (nDH in c(25,96)){
All.df<-read.csv(paste0("nDH",nDH,"_All.df","_",values,".csv"),sep=",",header=TRUE)

All.df %>%
  group_by(SelectSP) %>%                        # Specify group indicator (TestSP/SelectSP/Year)
summarise_at(vars(Mean),              # Specify column
               list(name = mean))    ###!!! list(meant=mean)
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
  