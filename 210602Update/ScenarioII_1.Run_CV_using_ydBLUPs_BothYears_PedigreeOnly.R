# Run Only Once
# Generate GCA and SCA within Both Years'
# Generate GCA and SCA var within 2019, then within 2020


######## Run GS CV Both Years !!
rm(list=ls())
library(BGLR)

WD<-"/local/workdir/mh865/GCA_SCA/PedigreeOnly_GCA_SCA/"  # run in terminal
datafdr<-paste0(WD,"data/")

load(paste0(datafdr,"dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot
#load("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP/dataNHim_withChk_3_sets_PhotoScore0123.rdata")  ## Indi
####
Y1<-dataNHpiBoth_C
colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM","AshFreedryWgtPerM")
Y2<-droplevels(Y1[,colKeep])
  head(Y2)
Y<-Y2 # Both Years
Y<-droplevels(Y[Y$popChk=="ES",])

mm.file<- paste0(datafdr,"A.csv")          # path to covariates file

Inputfiledir<-c("OneTime1920/GP1/","OneTime1920/GP2/","OneTime1920/GP1P2/") 

load(paste0(WD,Inputfiledir[1],"EVD.rda"))  
EVD1<-EVD
rm(EVD)
load(paste0(WD,Inputfiledir[2],"EVD.rda"))       
EVD2<-EVD
rm(EVD)
load(paste0(WD,Inputfiledir[3],"EVD.rda"))       
EVD3<-EVD
rm(EVD)
ETA<-list(
  list(V=EVD1$vectors,d=EVD1$values,model="RKHS"),
  list(V=EVD2$vectors,d=EVD2$values,model="RKHS"),
  list(V=EVD3$vectors,d=EVD3$values,model="RKHS")
)

#####!!!


load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_overTwoYears_AddBD.Rdata")) ##!!!
rownames(Both_dBLUPs)<-Both_dBLUPs$Row.names
Both_dBLUPs<-Both_dBLUPs[,-1]
traits<-colnames(Both_dBLUPs)
  print(traits)
CrossBLUE<-Both_dBLUPs ##!!!
##### !!!


  
load(paste0("/local/workdir/mh865/GCA_SCA/OneTime1920/data/sampleCV_250Indiv_0303_2021.Rdata"))

sampleCV<-sampleCV
folds   <- 1:10 
reps<-20 # !!!
ntraits<- length(traits) # !!!
cor<-matrix(nrow=reps,ncol=ntraits)
colnames(cor)<-traits


for (j in 1:length(traits)){
  Coltrait<-traits[j]

  Y$BLUE_Trait<-expss::vlookup(Y$Crosses,dict=CrossBLUE,result_column = paste0(Coltrait),lookup_column = "row.names") ### This is the DwPM
  head(Y)
  
  head(Y)
  dim(Y)
  
  Y<-Y   # Between Year !!!!
  y<-Y[,"BLUE_Trait"]  # phenotypes column  !!!

  yBLUE<-Y[,"BLUE_Trait"] # This is the BLUE for DwPM
  
  setwd(paste0(WD,"OneTime1920/Alldata_CV_output/"))


  for (i in 1:reps){
    setwd(paste0(WD,"OneTime1920/Alldata_CV_output/"))
    dir.create(paste0(Coltrait,"_ydBLUPsnoLocRep",i))
    savepath<-paste0("OneTime1920/Alldata_CV_output/",Coltrait,"_ydBLUPsnoLocRep",i,"/")  # the path within WD!
    
    tmp<-NULL
    for (fold in folds){
      
      yNA<-y
      testing<-which(sampleCV[,i]==fold)
      yNA[testing]<-NA
      
      fm<-BGLR(y=yNA,
               ETA=ETA,
               nIter=80000,
               burnIn=60000,
               saveAt=paste0(WD,savepath,"CVData1920_",fold,"thfold_rep",i),
               verbose=TRUE)
      save(fm,file=paste0(WD,savepath,"fm_",fold,"thfold_rep",i,".rda"))
      
      yPred<-fm$ETA[[3]]$u+fm$ETA[[2]]$u+fm$ETA[[1]]$u  #SCA+GCA2+GCA1
      predict<-data.frame(testing,
                          Crosses=Y[c(testing),]$Crosses,
                          yBLUE=yBLUE[testing],
                          yPred=yPred[testing],
                          yHat=fm$yHat[testing],
                          popChk=Y[c(testing),"popChk"],
                          Year=Y[c(testing),]$Year)
      predict<-droplevels(predict)
      
      tmp<-rbind(tmp,predict)
     
    }
    cor[i,j]<-cor(tmp$yBLUE,tmp$yPred,use="complete")
    
    write.table(tmp,file=paste0(WD,savepath,"predictions_rep",i,".csv"),row.names=FALSE,sep=",") 
    rm(fm) 
    unlink("*.dat")
  } 
  
}

standard_err<-function(x){sd(x,na.rm=TRUE)/sqrt(length(na.omit(x)))}
stderr<-apply(cor,2,standard_err)
cormean<-colMeans(cor)
cor_std<-rbind(cormean,stderr)
rownames(cor_std)<-c("corMean","StdErr")
colnames(cor_std)<-traits

write.csv(cor_std,paste0(paste0(WD,"OneTime1920/Alldata_CV_output/","cor_CV_no_Loc_BothYears_ydrBLUPs_data_",length(traits),"Traits_Mean_05272021.csv")))
write.csv(cor,paste0(WD,"OneTime1920/Alldata_CV_output/","cor_CV_no_Loc_BothYears_ydrBLUPs_data_",length(traits),"Traits_05272021.csv"))










#plot<-ggplot(data=Y,aes(BLUE_Trait,Year))+
#  geom_point(aes(color=as.factor(Year)))+ 
#  geom_line(aes(group=as.factor(Crosses)))
#print(plot)
  
### Run it only once
# # #the P1 column. column in phenotype file that gives the IDs that link observations to covariates or grouping factor
# CalGCA_noChk(Y=Y,mm.file=mm.file,colIDy = 3,colNam = "P1",savefiledir = "OneTime1920/GP1/250Individual/") 
# CalGCA_noChk(Y=Y,mm.file=mm.file,colIDy = 5,colNam = "P2",savefiledir = "OneTime1920/GP2/250Individual/") 
# # 
#  CalSCA(G1.file=paste0(WD,"OneTime1920/GP1/250Individual/","G.rda"),
#         G2.file=paste0(WD,"OneTime1920/GP2/250Individual/","G.rda"),
#         savefileDir="OneTime1920/GP1P2/250Individual/")
# # ############# Run it only once
# # 

# colMeans(cor)  # All 283 individuals from tmp   #noLOC 0.2696567

# 
# cor2<-matrix(nrow=reps,ncol=1)
# for (i in 1:reps){
#   savepath<-paste0("OneTime1920/Alldata_CV_output/y2020BLUEnoLocRep",i,"/")
#   predict<-read.csv(paste0(WD,savepath,"predictions_rep",i,".csv"),sep=",",header=TRUE)
#   cor2[i,]<-cor(predict[predict$popChk=="ES",]$yBLUE,predict[predict$popChk=="ES",]$yPred,use="complete")
# }
# colMeans(cor2) #0.1827861     #noLOC 0.2699811
# write.csv(cor2,"10foldCV_cor2_using_250ES_plots_NoFMLoc_y2020BLUEas_y.csv")
# 
# #
# cor2<-matrix(nrow=reps,ncol=1)
# for (i in 1:reps){
#   savepath<-paste0("OneTime1920/Alldata_CV_output/percDryWgt_yBLUEnoLocRep",i,"/")
#   predict<-read.csv(paste0(WD,savepath,"predictions_rep",i,".csv"),sep=",",header=TRUE)
#   cor2[i,]<-cor(predict[predict$popChk=="ES",]$yBLUE,predict[predict$popChk=="ES",]$yPred,use="complete")
# }
# 

