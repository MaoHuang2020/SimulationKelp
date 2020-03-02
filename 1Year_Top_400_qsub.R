
load("founderPop.Rdata")


scheme="1_Year_Top_400"  ####### !!!!!!!!!

library(AlphaSimR)

	nSnpPerChr<-100
	nQtlPerChr<-100
	
	mean_Trait<-0  ###!!! founder Pop trait
	var_Trait<-1   ###!!! founder Pop trait
	varE<-2        ###!!! Using varE instead of herit
	
	nCrosses<-400  ###!!! All set to this except for "Larger plots"
	
	SP = SimParam$new(founderPop)
	SP$addTraitA(nQtlPerChr=nQtlPerChr, mean=mean_Trait,var=var_Trait,corA=NULL)  ### Input corA? for multiple traits??
	
	SP$addSnpChip(nSnpPerChr=nSnpPerChr)  ###Randomly assign eligble SNPs to a SNP chip
	SP$setGender("yes_sys")	### !!! systemic assign indi to: male-female-male-female
	SP$setTrackRec(TRUE)
	
	pop <- newPop(founderPop, simParam=SP)
	pop<-setPheno(pop,varE=varE,simParam=SP) ### !!!
	genMean=meanG(pop) 
	  
	generation<- vector(length=6,mode="list") ### !!!
	generation[[1]]<-pop



nrep<-20  ### !!!
Mean_g_Rep<-matrix(nrow=length(generation),ncol=nrep)
Sd_g_Rep<-matrix(nrow=length(generation),ncol=nrep)

for (i in 1:nrep){	

GP0_DH<-makeDH(generation[[1]],nDH=1,simParam=SP)	### !!! When n_gp=200, use 1
	GP0_DH_all<-pullQtlGeno(GP0_DH)	
	dim(GP0_DH_all)
	
n_gp<-200 ## !!!



## Randomly select 200 as female, then randomly select another 200 as male
females<-selectInd(pop=GP0_DH,gender="F",nInd=n_gp,trait=1,use="rand",simParam=SP)
males<-selectInd(pop=GP0_DH,gender="M",nInd=n_gp,trait=1,use="rand",simParam=SP)	

GP_F<-row.names(pullQtlGeno(females))
GP_M<-row.names(pullQtlGeno(males))

GP_F_2<-GP_F[sample(length(GP_F))]	  					 
GP_M_2<-GP_M[sample(length(GP_M))] 					 

GP_Fs<-c(GP_F,GP_F_2)
GP_Ms<-c(GP_M,GP_M_2)

crossPlan<-cbind(GP_Fs,GP_Ms)

## Cross to create SP1 and phenotype them
SP1<-makeCross(GP0_DH,crossPlan=crossPlan,simParam=SP)	
SP1<-setPheno(SP1,varE=2,simParam=SP)


###	
## !!! OR			 
SP1_s<-selectInd(SP1,nInd=100,trait=1,selectTop=TRUE,use="pheno",simParam=SP)	
###	

## Make GP DH using SP1
GP1_DH<-makeDH(SP1_s,nDH=25,simParam=SP)	### !!!

## GS and GEBV on the GP1_DHs
GS1<-RRBLUP(SP1,traits=1,simParam=SP)
GP1_DH_GEBV<-setEBV(GP1_DH,GS1,simParam=SP)

#### Year 2
## Select GP1s, generate SP2s
females<-selectInd(pop=GP1_DH_GEBV,gender="F",nInd=n_gp,trait=1,use="ebv",selectTop=TRUE,simParam=SP)

males<-selectInd(pop=GP1_DH_GEBV,gender="M",nInd=n_gp,trait=1,use="ebv",selectTop=TRUE,simParam=SP)
	
GP_F<-row.names(pullQtlGeno(females))
GP_M<-row.names(pullQtlGeno(males))

GP_F_2<-GP_F[sample(length(GP_F))]	  					 
GP_M_2<-GP_M[sample(length(GP_M))] 					 

GP_Fs<-c(GP_F,GP_F_2)
GP_Ms<-c(GP_M,GP_M_2)

crossPlan<-cbind(GP_Fs,GP_Ms)

## Cross to create SP2, phenotype SP2
SP2<-makeCross(GP1_DH_GEBV,crossPlan=crossPlan,simParam=SP)	
SP2<-setPheno(SP2,varE=2,simParam=SP)


###	
## !!! OR			 
SP2_s<-selectInd(SP2,nInd=100,trait=1,selectTop=TRUE,use="pheno",simParam=SP)	
###	

## Make GP2 DH using SP2_s
GP2_DH<-makeDH(SP2_s,nDH=25,simParam=SP)	### !!!

## GS and GEBV on the GP2_DHs
popList=list(SP1,SP2)
SP1_SP2<-mergePops(popList)
GS2<-RRBLUP(SP1_SP2,traits=1,simParam=SP)
GP2_DH_GEBV<-setEBV(GP2_DH,GS2,simParam=SP)

##### Year 3

## Select GP2s
females<-selectInd(pop=GP2_DH_GEBV,gender="F",nInd=n_gp,trait=1,use="ebv",selectTop=TRUE,simParam=SP)

males<-selectInd(pop=GP2_DH_GEBV,gender="M",nInd=n_gp,trait=1,use="ebv",selectTop=TRUE,simParam=SP)
	
GP_F<-row.names(pullQtlGeno(females))
GP_M<-row.names(pullQtlGeno(males))

GP_F_2<-GP_F[sample(length(GP_F))]	  					 
GP_M_2<-GP_M[sample(length(GP_M))] 					 

GP_Fs<-c(GP_F,GP_F_2)
GP_Ms<-c(GP_M,GP_M_2)

crossPlan<-cbind(GP_Fs,GP_Ms)

## Cross to create SP3, phenotype SP3
SP3<-makeCross(GP2_DH_GEBV,crossPlan=crossPlan,simParam=SP)	
SP3<-setPheno(SP3,varE=2,simParam=SP)


###	
## !!! OR			 
SP3_s<-selectInd(SP3,nInd=100,trait=1,selectTop=TRUE,use="pheno",simParam=SP)	
###	

## Make GP3 DH using SP3_s
GP3_DH<-makeDH(SP3_s,nDH=25,simParam=SP)	### !!!

## GS and GEBV on the GP3_DHs
popList=list(SP1,SP2,SP3)
SP1_SP2_SP3<-mergePops(popList)
GS3<-RRBLUP(SP1_SP2_SP3,traits=1,simParam=SP)
GP3_DH_GEBV<-setEBV(GP3_DH,GS3,simParam=SP)

##### Year 4

## Select GP3s
females<-selectInd(pop=GP3_DH_GEBV,gender="F",nInd=n_gp,trait=1,use="ebv",selectTop=TRUE,simParam=SP)

males<-selectInd(pop=GP3_DH_GEBV,gender="M",nInd=n_gp,trait=1,use="ebv",selectTop=TRUE,simParam=SP)
	
GP_F<-row.names(pullQtlGeno(females))
GP_M<-row.names(pullQtlGeno(males))

GP_F_2<-GP_F[sample(length(GP_F))]	  					 
GP_M_2<-GP_M[sample(length(GP_M))] 					 

GP_Fs<-c(GP_F,GP_F_2)
GP_Ms<-c(GP_M,GP_M_2)

crossPlan<-cbind(GP_Fs,GP_Ms)

## Cross GP3s to create SP4, phenotype SP4
SP4<-makeCross(GP3_DH_GEBV,crossPlan=crossPlan,simParam=SP)	
SP4<-setPheno(SP4,varE=2,simParam=SP)


###	
## !!! OR			 
SP4_s<-selectInd(SP4,nInd=100,trait=1,selectTop=TRUE,use="pheno",simParam=SP)	
###	

## Make GP4 DH using SP4_s
GP4_DH<-makeDH(SP4_s,nDH=25,simParam=SP)	### !!!

## GS and GEBV on the GP3_DHs
popList=list(SP1,SP2,SP3,SP4)
SP1_SP2_SP3_SP4<-mergePops(popList)
GS4<-RRBLUP(SP1_SP2_SP3_SP4,traits=1,simParam=SP)
GP4_DH_GEBV<-setEBV(GP4_DH,GS4,simParam=SP)


###### Year 5
## Select GP4s
females<-selectInd(pop=GP4_DH_GEBV,gender="F",nInd=n_gp,trait=1,use="ebv",selectTop=TRUE,simParam=SP)

males<-selectInd(pop=GP4_DH_GEBV,gender="M",nInd=n_gp,trait=1,use="ebv",selectTop=TRUE,simParam=SP)
	
GP_F<-row.names(pullQtlGeno(females))
GP_M<-row.names(pullQtlGeno(males))

GP_F_2<-GP_F[sample(length(GP_F))]	  					 
GP_M_2<-GP_M[sample(length(GP_M))] 					 

GP_Fs<-c(GP_F,GP_F_2)
GP_Ms<-c(GP_M,GP_M_2)

crossPlan<-cbind(GP_Fs,GP_Ms)

## Cross GP4s to create SP5, phenotype SP5
SP5<-makeCross(GP4_DH_GEBV,crossPlan=crossPlan,simParam=SP)	
SP5<-setPheno(SP5,varE=2,simParam=SP)


###	
## !!! OR			 
SP5_s<-selectInd(SP5,nInd=100,trait=1,selectTop=TRUE,use="pheno",simParam=SP)	


## Make GP5 DH using SP5_s
GP5_DH<-makeDH(SP5_s,nDH=25,simParam=SP)	### !!!

# # # ## GS and GEBV on the GP5_DHs
# # # popList=list(SP1,SP2,SP3,SP4,SP5)
# # # SP1_SP2_SP3_SP4_SP5<-mergePops(popList)
# # # GS5<-RRBLUP(SP1_SP2_SP3_SP4_SP5,traits=1,simParam=SP)
# # # GP5_DH_GEBV<-setEBV(GP5_DH,GS5,simParam=SP)

### !!! ??? SP#_s 
generation[[2]]<-SP1
generation[[3]]<-SP2
generation[[4]]<-SP3
generation[[5]]<-SP4
generation[[6]]<-SP5

mean_g1<-unlist(lapply(generation,meanG))
sd_g1<-sqrt(unlist(lapply(generation,varG)))

Mean_g_Rep[,i]<-mean_g1
Sd_g_Rep[,i]<-sd_g1
}

Mean_g<-rowMeans(Mean_g_Rep)
Sd_g<-rowMeans(Sd_g_Rep)

Mean_Sd<-cbind(Mean_g,Sd_g)

write.csv(Mean_g_Rep,paste(scheme,"_Mean_g_WithRep.csv",sep=""))
write.csv(Sd_g_Rep,paste(scheme,"_Sd_g_WithRep.csv",sep=""))
write.csv(Mean_Sd,paste(scheme,"_Mean_g_Sd_AverageOverRep.csv",sep=""))



# # # # # # library(ggplot2)
# # # dev.set()
# # # pdf("PSRand400_1-Yr.pdf")
# # # plot_gain <- qplot(x = 0:5,
	                   # # # y = mean_g6,
	                   # # # ymin = mean_g6 - sd_g6,
	                   # # # ymax = mean_g6 + sd_g6,
	                   # # # geom = "pointrange",
	                   # # # main = "Genetic mean and standard deviation",
	                   # # # xlab = "Cycles", ylab = "Genetic mean")
# # # print(plot_gain)
# # # dev.off()	
	
