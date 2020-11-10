### Making plots
### Read in the files 1yr Rand, 2Yr Rand, 1 Yr Top, 2Yr Top
### read in and store in list, then merge to a table
# 
values<-"Sd"
nDH<-24
varE<-1
Ne<-60


rm(list=ls())
ls()
#setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/Run_Function/Different_varE_nDH/Combine_Run_1,2_year/Output_from_bioHPC")

setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/201030Update/100Cycles")

filenames <- list.files(full.names=TRUE)  
filenames

Mean_SP<-NULL
SchemePlot<-NULL

for (values in c("Mean","Sd")){
  for (nDH in c(24,96)){
    for (varE in c(1,4)){
      for (Ne in c(60,600)){
        MeanSP<-filenames[grepl(paste0("nDH",nDH,"_varE",varE,"_Ne",Ne,"_",values,"_GP.csv"),filenames)]  ## Used GP
        SchemePlot<-paste0("nDH",nDH,"_varE",varE,"_Ne",Ne)
        MeanSP
        All1 <- lapply(MeanSP,function(i){
          read.csv(i, header=TRUE,row.names=1)
        })
        
        class(All1)
        length(All1)
        dim(All1[[1]])
        All1[[1]][1:5,]
        
        cycles=7  ### 7 years of breeding cycles
        #Mean_SD<-matrix(nrow=8*cycles,ncol=3)
        
        Mean_SD<-data.frame(matrix(nrow=8*cycles,ncol=3))
        MeanSP ### order stored in All1 
        ### order<-c("Pheno_1000_1yr","Pheno_1000_2yr","Pheno_400_1yr","Pheno_400_2yr","Rand_1000_1yr","Rand_1000_2yr","Rand_400_1yr","Rand_400_2yr")
        
        Mean_SP<-NULL
        sd_SP<-NULL
        select_order<-NULL
        
        for (i in 1:length(All1)){
          Mean_SP<-c(Mean_SP,rowMeans(All1[[i]]))  # Cal Mean
          sd_SP<-c(sd_SP,apply(All1[[i]],1, function(x){sd(x)/sqrt(length(x))})) # Cal St_error
          select_order<-c(select_order,rep(MeanSP[i],cycles))
        }
        
        Mean_SD[,1]<-Mean_SP
        Mean_SD[,2]<-sd_SP
        Mean_SD[,3]<-as.factor(select_order)
        
        str(Mean_SD)
        
        colnames(Mean_SD)<-c("Mean","StdErr","Selection")  #StdErr
        Mean_SD[1:5,]
        tail(Mean_SD)
        
        # install.packages("stringr")
        # Mean_SD[,1]<-as.numeric(Mean_SD[,1])
        # Mean_SD[,2]<-as.numeric(Mean_SD[,2])
        
        str(Mean_SD)  ##### !!!! Ensure the values are correct
        
        library(stringr)
        Mean_SD$Shorten<-sub("./","",Mean_SD$Selection) 
        Mean_SD$Shorten<-sub(".csv","",Mean_SD$Shorten) 
        StringSplit<-str_split_fixed(string=Mean_SD$Shorten,"_",7) ### This text string becomes a 56x7 matrix
        Mean_SD$SelectSP<-as.factor(StringSplit[,1])
        Mean_SD$TestSP<-as.factor(StringSplit[,2])
        Mean_SD$Year<-as.factor(StringSplit[,3])
        
        ###OR Mean_SD$SelectSPtest<-str_split(string=Mean_SD$Shorten,"[^a-zA-Z0-9]")[[1:nrow(Mean_SD)]][1]
        
        Mean_SD[1:5,]
        Mean_SD$Cycles<-rep(1:7,length(All1))
        head(Mean_SD)
        tail(Mean_SD)
        str(Mean_SD)
        
        write.csv(Mean_SD,paste(SchemePlot,"_",values,".csv",sep=""))
        
        library(ggplot2)
        
        dev.set()
        tiff(file=paste(SchemePlot,".tiff",sep=""),width=1400,height=1000,units="px",pointsize=12,res=150)
        
        plot<-ggplot(Mean_SD,mapping=aes(x=Cycles,y=Mean))+
          geom_point(aes(shape=TestSP))+
          
          geom_line(aes(group=Shorten,color=Year,linetype=SelectSP))+
          theme_bw()+
          labs(x="Year",y="SP genetic mean")+ 
          ggtitle(SchemePlot)
        
        ######plot+geom_segment(arrow=arrows(x0 = Cycles-0.01, y0 = Mean, x1 = Cycles+0.01, y1 = Mean+SD, angle=90,colour = "gray"))
        
        plot+scale_color_manual(breaks = c("1yr","2yr"),
                                values=c("orangered", "gray17")) +
          scale_shape_manual(values=c(3,16))
        
        dev.off()
        
      
        }
    }
  }
}












# rm=list(ls())
# setwd("//Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/Run_Function/Different_varE_nDH/nDH25/varE5.67_Ne600/")
# SchemePlot<-"varE5.67_Ne600_nDH25" ### !!! Need change everytime
#MeanSP<-filenames[grepl("nDH25_varE1.22_Ne60_.*Mean_g.csv",filenames)]

#MeanSP<-grep("Mean_g.csv",filenames) ## !!!!! Just the Mean_g.csv
#values<-"Mean"

# MeanSP<-grep("Sd_g.csv",filenames) ## !!!!! Just the Sd_g.csv
# values<-"Sd"


########### Compile all the data files?

# library(RColorBrewer)
# plot+scale_color_brewer(palette="Set2")
#group=interaction(Year,SelectSP)
#scale_x_continuous(breaks=c(1:7)
#scale_shape_manual(values=c(1,2))+

# ### This is to add the error bars
#   geom_segment(aes(x=Cycles,y=Mean,xend=Cycles,yend=Mean+SD))+
#   geom_segment(aes(x=Cycles-0.02,y=Mean+SD,xend=Cycles+0.02,yend=Mean+SD))

# year<-"1yr"
# selection="rand";nPheno=400;nDH=96;varE=5.67;Ne=60;
# scheme1yr.1<-paste(selection,"_",nPheno,"_",year,"_nDH",nDH,"_varE",varE,"_","Ne",Ne,sep="")
# 
# selection="rand";nPheno=1000;nDH=96;varE=5.67;Ne=60;
# scheme1yr.2<-paste(selection,"_",nPheno,"_",year,"_nDH",nDH,"_varE",varE,"_","Ne",Ne,sep="")
# 
# selection="pheno";nPheno=400;nDH=96;varE=5.67;Ne=60;
# scheme1yr.3<-paste(selection,"_",nPheno,"_",year,"_nDH",nDH,"_varE",varE,"_","Ne",Ne,sep="")
# 
# selection="pheno";nPheno=1000;nDH=96;varE=5.67;Ne=60;
# scheme1yr.4<-paste(selection,"_",nPheno,"_",year,"_nDH",nDH,"_varE",varE,"_","Ne",Ne,sep="")
# 
# 
# selection="pheno";nPheno=400;nDH=96;varE=5.67;Ne=600;
# selection="rand";nPheno=400;nDH=96;varE=5.67;Ne=600;
# selection="rand";nPheno=1000;nDH=96;varE=5.67;Ne=600;
# selection="pheno";nPheno=1000;nDH=96;varE=5.67;Ne=600;
# 
# 
# selection="rand";nPheno=400;nDH=96;varE=5.67;Ne=60;
# selection="rand";nPheno=1000;nDH=96;varE=5.67;Ne=60;
# selection="pheno";nPheno=400;nDH=96;varE=5.67;Ne=60;
# selection="pheno";nPheno=1000;nDH=96;varE=5.67;Ne=60;
# 
# selection="pheno";nPheno=400;nDH=96;varE=5.67;Ne=600;
# selection="rand";nPheno=400;nDH=96;varE=5.67;Ne=600;
# selection="rand";nPheno=1000;nDH=96;varE=5.67;Ne=600;
# selection="pheno";nPheno=1000;nDH=96;varE=5.67;Ne=600;
# 
# 
# set<-c("C1_YLD","C2_YLD","C3_YLD","C1_HGT","C2_HGT","C1_HD","C2_HD","C3_HD")
# 
# file_list <- list.files(path="C:/Users/Luke/Documents/NBA")
# 
# #initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
# dataset <- data.frame()
# 
# MeanSP=length(set) ##
# All=NULL
# for (i in 1:length(set)){
#   filename<-paste(scheme,"_Mean_g.csv",sep="")
#   tmpfile <-read.table(file=filename,sep=",",header=TRUE)
#   assign(paste("Mean",i,sep="_"),tmpfile)
#   MeanSP[[i]]<-tmpfile
# }
