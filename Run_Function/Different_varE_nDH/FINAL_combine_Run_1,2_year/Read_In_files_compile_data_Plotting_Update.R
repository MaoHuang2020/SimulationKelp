###### STEP 1, Combine and get average values into one dataframe for each scenario

### Making plots
### Read in the files 1yr Rand, 2Yr Rand, 1 Yr Top, 2Yr Top
### read in and store in list, then merge to a table


rm(list=ls())
ls()
#setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/201030Update/100Cycles")

### meanG and varG
setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/210602Update/100Cycles")

filenames <- list.files(full.names=TRUE)  
filenames
Mean_SP<-NULL  ## Store the "average" value (over reps) for meanG or varG
SchemePlot<-NULL

for (values in c("Mean","Sd")){  #meanG, varG from AlphaSimR
  for (nDH in c(24,96)){
    for (varE in c(1,4)){
      for (Ne in c(60,600)){
        MeanSP<-filenames[grepl(paste0("nDH",nDH,"_varE",varE,"_Ne",Ne,"_",values,"_GP.csv"),filenames)]  ## Used GP
        SchemePlot<-paste0("nDH",nDH,"_varE",varE,"_Ne",Ne)

        All1 <- lapply(MeanSP,function(i){
          read.csv(i, header=TRUE,row.names=1)
        })
        
        cycles=10  ### 10 years of breeding cycles

        Mean_SD<-data.frame(matrix(nrow=8*cycles,ncol=3))  # 8 different combinations
        MeanSP ### order stored in All1 

        Mean_SP<-NULL
        sd_SP<-NULL
        select_order<-NULL
        
        for (i in 1:length(All1)){
          Mean_SP<-c(Mean_SP,rowMeans(All1[[i]]))  # Cal Mean Genetic values
          sd_SP<-c(sd_SP,apply(All1[[i]],1, function(x){sd(x)/sqrt(length(x))})) # Cal St_error for meanG
          select_order<-c(select_order,rep(MeanSP[i],cycles))
        }
        
        Mean_SD[,1]<-Mean_SP
        Mean_SD[,2]<-sd_SP
        Mean_SD[,3]<-as.factor(select_order)
        
         str(Mean_SD)
        
        colnames(Mean_SD)<-c("Mean","StdErr","Selection")  #StdErr
          Mean_SD[1:5,]
          tail(Mean_SD)

        library(stringr)
        Mean_SD$Shorten<-sub("./","",Mean_SD$Selection) 
        Mean_SD$Shorten<-sub(".csv","",Mean_SD$Shorten) 
        StringSplit<-str_split_fixed(string=Mean_SD$Shorten,"_",cycles) ### This text string becomes a 56xcycles matrix
        Mean_SD$SelectSP<-as.factor(StringSplit[,1])
        Mean_SD$TestSP<-as.factor(StringSplit[,2])
        Mean_SD$Year<-as.factor(StringSplit[,3])
        
        ###OR Mean_SD$SelectSPtest<-str_split(string=Mean_SD$Shorten,"[^a-zA-Z0-9]")[[1:nrow(Mean_SD)]][1]

        Mean_SD$Cycles<-rep(1:cycles,length(All1))

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




#### Add in the selection intensity on GP and SP, also GS cor on GP
setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/210602Update/100Cycles")
filenames <- list.files(full.names=TRUE)  
filenames
outputfdr<-"/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/210602Update/100Cycles/output/"

CalcMean_sd<-function(values="SPselInt",cycles=10,outputfdr=outputfdr){

for (nDH in c(24,96)){
  for (varE in c(1,4)){
    for (Ne in c(60,600)){
        
      FileInput<-filenames[grepl(paste0("nDH",nDH,"_varE",varE,"_Ne",Ne,"_",values,"_Reps.csv"),filenames)] # 8 cycles
      SchemePlot<-paste0("nDH",nDH,"_varE",varE,"_Ne",Ne)
      
      All1<-lapply(FileInput,function(i){
        read.csv(i,header=TRUE,row.names=1)
      })
      
      cycles=cycles    ### GPs Only had 8 years of breeding cycles

      Mean_sd<-data.frame(matrix(nrow=8*cycles,ncol=3))    # 8 different combinations
      FileInput   ### order stored in All1 
      ### order<-c("Pheno_1000_1yr","Pheno_1000_2yr","Pheno_400_1yr","Pheno_400_2yr","Rand_1000_1yr","Rand_1000_2yr","Rand_400_1yr","Rand_400_2yr")
      
      Mean_File<-NULL
      sd_File<-NULL
      select_order<-NULL
      
      for (i in 1:length(All1)){
        Mean_File<-c(Mean_File,rowMeans(All1[[i]]))  # Cal Mean Genetic values
        sd_File<-c(sd_File,apply(All1[[i]],1, function(x){sd(x)/sqrt(length(x))})) # Cal St_error
        select_order<-c(select_order,rep(FileInput[i],cycles))
      }
      
      Mean_sd[,1]<-Mean_File
      Mean_sd[,2]<-sd_File
      Mean_sd[,3]<-as.factor(select_order)
      
      str(Mean_sd)  ##### !!!! Ensure the values are correct
      
      colnames(Mean_sd)<-c("Mean","StdErr","Selection")  #average value, StdErr, category
      Mean_sd[1:5,]
      tail(Mean_sd)
      
      library(stringr)
      Mean_sd$Shorten<-sub("./","",Mean_sd$Selection) 
      Mean_sd$Shorten<-sub(".csv","",Mean_sd$Shorten) 
      StringSplit<-str_split_fixed(string=Mean_sd$Shorten,"_",cycles) ### This text string becomes a 56xcycles matrix
      Mean_sd$SelectSP<-as.factor(StringSplit[,1])
      Mean_sd$TestSP<-as.factor(StringSplit[,2])
      Mean_sd$Year<-as.factor(StringSplit[,3])
      
      Mean_sd[1:5,]
      Mean_sd$Cycles<-rep(1:cycles,length(All1))
      head(Mean_sd)
      tail(Mean_sd)
      str(Mean_sd)
      
      write.csv(Mean_sd,paste(outputfdr,SchemePlot,"_",values,".csv",sep=""))
      
       }
      }
    }
  }


CalcMean_sd(values="SPselInt",cycles=10,outputfdr=outputfdr)
CalcMean_sd(values="GPselInt",cycles=8,outputfdr=outputfdr)
CalcMean_sd(values="GPcor",cycles=8,outputfdr=outputfdr)

