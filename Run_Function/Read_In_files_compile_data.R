### Making plots

### Read in the files 1yr Rand, 2Yr Rand, 1 Yr Top, 2Yr Top

### read in and store in list, then merge to a table


#########
#####Starting from here:
setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/Run_Function/Different_varE_nDH/varE2_Ne60_nDH96/")

rm=list(ls())
setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SimulationKelp/Run_Function/Different_varE_nDH/varE2_Ne600_nDH96/")

filenames <- list.files(full.names=TRUE)  
str(filenames)

MeanSP<-grep("Mean_g.csv",filenames) ## Just the Mean_g.csv

All1 <- lapply(filenames[MeanSP],function(i){
  read.csv(i, header=TRUE,row.names=1)
})

class(All1)
length(All1)
dim(All1[[1]])
All1[[1]][1:5,1:6]

cycles=7
Mean_SD<-matrix(nrow=8*cycles,ncol=3)
colnames(Mean_SD)<-c("Mean","SD","Selection")

filenames[MeanSP] ### stored in All1
order<-c("Pheno_1000_1yr","Pheno_1000_2yr","Pheno_400_1yr","Pheno_400_2yr","Rand_1000_1yr","Rand_1000_2yr","Rand_400_1yr","Rand_400_2yr")

Mean_SP<-NULL
sd_SP<-NULL
select_order<-NULL

for (i in 1:length(All1)){
  Mean_SP<-c(Mean_SP,rowMeans(All1[[i]]))
  sd_SP<-c(sd_SP,apply(All1[[i]],1, sd, na.rm = TRUE))
  select_order<-c(select_order,rep(filenames[MeanSP][i],7))
}

Mean_SD[,1]<-as.numeric(Mean_SP)
Mean_SD[,2]<-as.numeric(sd_SP)
Mean_SD[,3]<-as.vector(select_order)

Mean_SD<-as.data.frame(Mean_SD)
Mean_SD[1:5,]
tail(Mean_SD)

install.packages("stringr")
library(stringr)
##!!!!!
Mean_SD$Shorten<-gsub("_nDH96_varE2_Ne60_Mean_g.csv","",Mean_SD$Selection) #### Need to update this!!!!
SchemePlot<-"nDH96_varE2_Ne60" ### !!! Need change too

Mean_SD$SelectSP<-gsub("_","",str_sub(Mean_SD$Shorten,3,7)) ### Substring: str_sub(data,starting point, ending point)
Mean_SD$TestSP<-gsub("_","",str_sub(Mean_SD$Shorten,-8,-5))  ### gsub("_","",data) remove "_" in data
Mean_SD$Year<-str_sub(Mean_SD$Shorten,-3,-1) 

Mean_SD$Cycles<-rep(1:7,length(All1))

library(ggplot2)


dev.set(1)
tiff(file=paste(SchemePlot,".tiff",sep=""),width=1400,height=1000,units="px",pointsize=12,res=150)

plot<-ggplot(Mean_SD,mapping=aes(x=Cycles,y=Mean))+
  geom_point(aes(shape=Year))+

  geom_line(aes(group=Shorten,color=interaction(TestSP,Year),linetype=SelectSP))+
  theme_bw()+
  labs(x="Year",y="SP genetic mean")
# plot+geom_segment(arrow=arrows(x0 = Cycles-0.01, y0 = Mean, x1 = Cycles+0.01, y1 = Mean+SD, angle=90,colour = "gray"))

plot+scale_color_manual(breaks = c("1000.1yr", "400.1yr", "1000.2yr","400.2yr"),
                        values=c("magenta2","darkturquoise","orangered", "gray17")) #",

dev.off()

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
