
# # # # # # 
# # # ## Estimating LD=0.0082
# # # LD=0.0082
# # # Ne=(1/(2*LD))-0.5
# # # Ne  ### 60.47561


#### Founder chromosomes
    ### The number of segregation sites to keep PER chr
    ### Sorghum #Physical lenght 1e8 base pairs
    ### Sorghum #Genetic length 1 Morgan
    ### mutRate=2.5/1E8*(4*Ne)      Sorghum  #2.5/1E8*(4*Ne), #Mutation rate adjusted for Ne
	### Sorghum:  1/1E8*(4*Ne), #Recombination rate adjusted for Ne
    ### Sorghum: 10/(4*Ne),100/Ne), #Modeling Ne=100 at 10 generations ago
    ### Within Reps for each founderPop/Rep ##!!!
    ### runMacs2() By default the species="Generic" #mutRate=mutRate,

#### Founder Pop trait
	#herit<-0.5     ###!!! founder SP heritability
    ###Using varE instead of herit   ### quick check on higher vs lower h2= 0.1, 0.3, 0.5 on patterns!!!
	
    ###!!! All set to this except for "Larger plots"
	
#### Set simulation parameters
    ### Input corA? for multiple traits??
    ###Randomly assign eligble SNPs to a SNP chip
	### !!! systemic assign indi to: male-female-male-female
    ### Create Founding sp population
    ### !!!


library(AlphaSimR)

#### Founder chromosomes
Ne=60
nInd<-800
nChr<-31
segSites<-500
SP_ploidy<-2L
bp<-1e+08
genLen=1

#### Founder Pop trait
nSnpPerChr<-100
nQtlPerChr<-100
mean_Trait<-0
var_Trait<-1
varE<-2  #### This is starting H=0.33  H=0.1 or H=0.5
nCrosses<-400


nrep<-20
Mean_g_Rep<-matrix(nrow=length(generation),ncol=nrep)
Sd_g_Rep<-matrix(nrow=length(generation),ncol=nrep)

for (i in 1:nrep){
    
  
    founderPop<-runMacs2(nInd=nInd,nChr=nChr,segSites=segSites,Ne=Ne,bp=bp,genLen=1,inbred=TRUE,ploidy=SP_ploidy, returnCommand = FALSE, nThreads = NULL)
    
    
    SP = SimParam$new(founderPop)
    SP$addTraitA(nQtlPerChr=nQtlPerChr, mean=mean_Trait,var=var_Trait,corA=NULL)
    SP$addSnpChip(nSnpPerChr=nSnpPerChr)
    SP$setGender("yes_sys")
    SP$setTrackRec(TRUE)
    
    pop <- newPop(founderPop, simParam=SP)
    pop<-setPheno(pop,varE=varE,simParam=SP)
    genMean=meanG(pop)
    
    generation<- vector(length=6,mode="list")
    generation[[1]]<-pop
    
    
schemes<-function(n_gp, selection)(
if(selection="rand"){
    if(n_gp=200){
    
 
        GP0_DH<-makeDH(generation[[1]],nDH=1,simParam=SP)
        females<-selectInd(pop=GP0_DH,gender="F",nInd=n_gp,trait=1,use="rand",simParam=SP)
        males<-selectInd(pop=GP0_DH,gender="M",nInd=n_gp,trait=1,use="rand",simParam=SP)
        
        GP_F<-row.names(pullQtlGeno(females))
        GP_M<-row.names(pullQtlGeno(males))
        
        GP_F_2<-GP_F[sample(length(GP_F))]
        GP_M_2<-GP_M[sample(length(GP_M))]
        
        GP_Fs<-c(GP_F,GP_F_2)
        GP_Ms<-c(GP_M,GP_M_2)
        
        crossPlan<-cbind(GP_Fs,GP_Ms)
        
        ## Cross to create SP1
        SP1<-makeCross(GP0_DH,crossPlan=crossPlan,simParam=SP)
        SP1<-setPheno(SP1,varE=2,simParam=SP)
        
        ## Random Select amongst SP1
        SP1_s<-selectInd(SP1,nInd=100,trait=1,use="rand",simParam=SP)
        
        ####
        ### !!! OR Select Top
        # SP1_s<-selectInd(SP1,nInd=100,trait=1,selectTop=TRUE,use="pheno",simParam=SP)
        ####
        
        ## Make GP DH using SP1
        GP1_DH<-makeDH(SP1_s,nDH=25,simParam=SP)	### !!!
        
        
        #### Year 2
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
        
        ## Cross to create SP2
        SP2<-makeCross(GP0_DH,crossPlan=crossPlan,simParam=SP)	
        SP2<-setPheno(SP2,varE=2,simParam=SP)
        
        ## Select amongst SP2 
        SP2_s<-selectInd(SP2,nInd=100,trait=1,use="rand",simParam=SP)					
        
        ### !!! OR			 
        # SP2_s<-selectInd(SP2,nInd=100,trait=1,selectTop=TRUE,use="pheno",simParam=SP)	
        ####	
        
        ## Make GP DH using SP2
        GP2_DH<-makeDH(SP2_s,nDH=25,simParam=SP)	### !!!
 
 GS<-NULL
 GEBV<-NULL
 popList<-list(SP1,SP2k)
 TP<-NULL
 
 SP_s<-NULL
 GP_DH<-NULL
 
 for (y in 3:5){
     popList=list(popList,SP[[y]])  ###
     
     #### Year 4
     ## Build GS model using SP1+SP2+SP3
     ## popList=list(SP1,SP2,SP3)
     
     
            TP[[y]]<-mergePops(popList,TP)
            GS[[y]]<-RRBLUP(TP,traits=1,simParam=SP)
            
            ## GEBVs on GP1s
            GEBV[[y]]<-setEBV(GP1_DH,GS1,simParam=SP)
            
            ## Select GP1s, generate SP3s
            females<-selectInd(pop=GEBV[[y]],gender="F",nInd=n_gp,selectTop=TRUE,trait=1,use="ebv",simParam=SP)
            
            males<-selectInd(pop=GEBV[[y]],gender="M",nInd=n_gp,selectTop=TRUE,trait=1,use="ebv",simParam=SP)
            
            GP_F<-row.names(pullQtlGeno(females))
            GP_M<-row.names(pullQtlGeno(males))
            
            GP_F_2<-GP_F[sample(length(GP_F))]
            GP_M_2<-GP_M[sample(length(GP_M))]
            
            GP_Fs<-c(GP_F,GP_F_2)
            GP_Ms<-c(GP_M,GP_M_2)
            
            crossPlan<-cbind(GP_Fs,GP_Ms)
            
            ## Cross to create SP3
            SP[[y]]<-makeCross(GP1_DH_GEBV,crossPlan=crossPlan,simParam=SP)
            SP[[y]]<-setPheno(SP[[y]],varE=2,simParam=SP)
            
            ## Select amongst SP3 
            SP_s[[y]]<-selectInd(SP[[y]],nInd=100,trait=1,use="rand",simParam=SP)
            
            ### !!! OR			 
            # SP3_s<-selectInd(SP3,nInd=100,trait=1,selectTop=TRUE,use="pheno",simParam=SP)	
            ####	
            
            ## Make GP DH using SP3
            GP_DH[[y]]<-makeDH(SP_s[[y]],nDH=25,simParam=SP)
            
        }
        
        
 }
        
    }else if(n_gp=500){
        
    }
    
}else if (selection="top"){
    if (n_gp=200){
        
    }else if(n_gp=500){
        
    }
}


)

		
###### 2-Year cycle
	
####### 400 Indi, 1GP/SP
	### !!! When n_gp=200, use 1
	GP0_DH_all<-pullQtlGeno(GP0_DH)	
	dim(GP0_DH_all)
####
n_gp<-200 ## !!!
####


# ####### 500 Indi, 2GP/SP
# GP0_DH<-makeDH(generation[[1]],nDH=2,simParam=SP)	### !!! When n_gp=5 needs to use 2			GP0_DH_all<-pullQtlGeno(GP0_DH)
#	dim(GP0_DH_all)
# ####
# n_gp<-500  ##!!!
# ####






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

## Cross to create SP1
SP1<-makeCross(GP0_DH,crossPlan=crossPlan,simParam=SP)	
SP1<-setPheno(SP1,varE=2,simParam=SP)

## Random Select amongst SP1
SP1_s<-selectInd(SP1,nInd=100,trait=1,use="rand",simParam=SP)					

####	
### !!! OR Select Top		 
# SP1_s<-selectInd(SP1,nInd=100,trait=1,selectTop=TRUE,use="pheno",simParam=SP)	
####	

## Make GP DH using SP1
GP1_DH<-makeDH(SP1_s,nDH=25,simParam=SP)	### !!!


#### Year 2
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

## Cross to create SP2
SP2<-makeCross(GP0_DH,crossPlan=crossPlan,simParam=SP)	
SP2<-setPheno(SP2,varE=2,simParam=SP)

## Select amongst SP2 
SP2_s<-selectInd(SP2,nInd=100,trait=1,use="rand",simParam=SP)					
	
### !!! OR			 
# SP2_s<-selectInd(SP2,nInd=100,trait=1,selectTop=TRUE,use="pheno",simParam=SP)	
####	

## Make GP DH using SP2
GP2_DH<-makeDH(SP2_s,nDH=25,simParam=SP)	### !!!

#### Year 3
## Build GS model using SP1+SP2
popList=list(SP1,SP2)
SP1_SP2<-mergePops(popList)
GS1<-RRBLUP(SP1_SP2,traits=1,simParam=SP)

## GEBVs on GP1s
GP1_DH_GEBV<-setEBV(GP1_DH,GS1,simParam=SP)

## Select GP1s, generate SP3s
females<-selectInd(pop=GP1_DH_GEBV,gender="F",nInd=n_gp,selectTop=TRUE,trait=1,use="ebv",simParam=SP)

males<-selectInd(pop=GP1_DH_GEBV,gender="M",nInd=n_gp,selectTop=TRUE,trait=1,use="ebv",simParam=SP)
	
GP_F<-row.names(pullQtlGeno(females))
GP_M<-row.names(pullQtlGeno(males))

GP_F_2<-GP_F[sample(length(GP_F))]	  					 
GP_M_2<-GP_M[sample(length(GP_M))] 					 

GP_Fs<-c(GP_F,GP_F_2)
GP_Ms<-c(GP_M,GP_M_2)

crossPlan<-cbind(GP_Fs,GP_Ms)

## Cross to create SP3
SP3<-makeCross(GP1_DH_GEBV,crossPlan=crossPlan,simParam=SP)	
SP3<-setPheno(SP3,varE=2,simParam=SP)

## Select amongst SP3 
SP3_s<-selectInd(SP3,nInd=100,trait=1,use="rand",simParam=SP)					
	
### !!! OR			 
# SP3_s<-selectInd(SP3,nInd=100,trait=1,selectTop=TRUE,use="pheno",simParam=SP)	
####	

## Make GP DH using SP3
GP3_DH<-makeDH(SP3_s,nDH=25,simParam=SP)	### !!!






#### Year 5 ## Build GS model using SP1+SP2+SP3+SP4

popList=list(SP1,SP2,SP3,SP4)
SP1_SP2_SP3_SP4<- mergePops(popList)
GS3<-RRBLUP(SP1_SP2_SP3_SP4,traits=1,simParam=SP)

## GEBVs on GP3s
GP3_DH_GEBV<-setEBV(GP3_DH,GS3,simParam=SP)

## Select GP3s, generate SP5s
females<-selectInd(pop=GP3_DH_GEBV,gender="F",nInd=n_gp,selectTop=TRUE,trait=1,use="ebv",simParam=SP)

males<-selectInd(pop=GP3_DH_GEBV,gender="M",nInd=n_gp,selectTop=TRUE,trait=1,use="ebv",simParam=SP)
	
GP_F<-row.names(pullQtlGeno(females))
GP_M<-row.names(pullQtlGeno(males))

GP_F_2<-GP_F[sample(length(GP_F))]	  					 
GP_M_2<-GP_M[sample(length(GP_M))] 					 

GP_Fs<-c(GP_F,GP_F_2)
GP_Ms<-c(GP_M,GP_M_2)

crossPlan<-cbind(GP_Fs,GP_Ms)

## Cross to create SP5
SP5<-makeCross(GP3_DH_GEBV,crossPlan=crossPlan,simParam=SP)

SP5<-setPheno(SP5,varE=2,simParam=SP)

## Select amongst SP5 
SP5_s<-selectInd(SP5,nInd=100,trait=1,use="rand",simParam=SP)					
	
### !!! OR			 
# SP3_s<-selectInd(SP3,nInd=100,trait=1,selectTop=TRUE,use="pheno",simParam=SP)	
####	

## Make GP4 DH using SP5
GP5_DH<-makeDH(SP5_s,nDH=25,simParam=SP)	### !!!

### !!! ??? SP#_s 
generation[[2]]<-SP1
generation[[3]]<-SP2
generation[[4]]<-SP3
generation[[5]]<-SP4
generation[[6]]<-SP5

mean_g1<-unlist(lapply(generation,meanG))  ### 2-Yr cycle, Select SP random, Pheno 400
sd_g1<-sqrt(unlist(lapply(generation,varG)))  ### 2-Yr cycle, Select SP random, Pheno 400

Mean_g_Rep[,i]<-mean_g1
Sd_g_Rep[,i]<-sd_g1
}

Mean_g<-rowMeans(Mean_g_Rep)
Sd_g<-rowMeans(Sd_g_Rep)

Mean_Sd<-cbind(Mean_g,Sd_g)

write.csv(Mean_g_Rep,paste(scheme,"_Mean_g_WithRep.csv",sep=""))
write.csv(Sd_g_Rep,paste(scheme,"_Sd_g_WithRep.csv",sep=""))
write.csv(Mean_Sd,paste(scheme,"_Mean_g_Sd_AverageOverRep.csv",sep=""))




# mean_g2<-unlist(lapply(generation,meanG))  ### 2-Yr cycle, Select SP random, Pheno 1000
# sd_g2<-sqrt(unlist(lapply(generation,varG)))  ### 2-Yr cycle, Select SP random, Pheno 1000


# # library(ggplot2)
# dev.set()
# pdf("PSRand1000_2-Yr.pdf")
# plot_gain <- qplot(x = 0:5,
	                   # y = mean_g2,
	                   # ymin = mean_g2 - sd_g2,
	                   # ymax = mean_g2 + sd_g2,
	                   # geom = "pointrange",
	                   # main = "Genetic mean and standard deviation",
	                   # xlab = "Cycles", ylab = "Genetic mean")
# print(plot_gain)
# dev.off()	
	
