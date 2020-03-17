
## Basic Parameters

library(AlphaSimR)
Ne=60       # Historical effective population size, estimated via E(r2)=1/(1+4*Ne*c)  ## If Ne=600, Pop str will inflat it
nInd<-1000   # Number of founders
nChr<-31    # Number of chr
segSites<-500   # Number of segregating sites per chromosome
Sporo_ploidy<-2L   # Ploidy
bp<-1e+08       # Base pair length of chromosome
genLen=1        # Genetic length of chromosome in Morgans

nSnpPerChr<-100  # Number of SNP per chromosome
nQtlPerChr<-100  # Number of QTL per chromosome
mean_Trait<-0    # Trait mean
var_Trait<-1     # Trait variance
varE<-5.67          # !!! This is starting H=0.33,  Try another H to represent Dry matter vs Wet biomass (0.15)

nDH0<-2           # Number of GPs per SP generated in the founder population to creat the inital GP0_DH

nDH<- 25         # !!! Or 96


## Define: Number of Sporo to Phenotype
###nPheno<-400  ##### !!!! Define
nPheno<-1000

n_gp<-nPheno/2

## Define: Reps and Cycles
nrep<-100  ## ?? Use shiny to render
cycles <-5


## Define: Sporo selection *Random* vs *Top*
### Select Random or Tops for the Sporophytes
s="Rand"    ######  !!!! Define 
## OR s="Top

Sel<-function (Population){
  if (s=="Rand"){
    Sp_select<-selectInd(Population,nInd=nPheno*0.1,trait=1,use="rand",simParam=SP)
  }else if (s =="Top"){
    Sp_select<-selectInd(Population,nInd=nPheno*0.1,trait=1,selectTop=TRUE,use="pheno",simParam=SP)	 
  }
  return(list=Sp_select) }

scheme<-paste(s,"_",nPheno,"_2yr",sep="")



## Running the reps and cycles
Mean_g_Rep<-matrix(nrow=cycles,ncol=nrep)
Sd_g_Rep<-matrix(nrow=cycles,ncol=nrep)

for (i in 1:nrep){
  founderPop<-runMacs2(nInd=nInd,nChr=nChr,segSites=segSites,Ne=Ne,bp=bp,genLen=1,inbred=TRUE,ploidy=Sporo_ploidy, returnCommand = FALSE, nThreads = NULL)
  # Founder Pop should be here
  
  SP<-SimParam$new(founderPop)
  SP$addTraitA(nQtlPerChr=nQtlPerChr, mean=mean_Trait,var=var_Trait,corA=NULL)
  SP$addSnpChip(nSnpPerChr=nSnpPerChr)
  SP$setGender("yes_sys")
  SP$setTrackRec(TRUE)
  
  pop <- newPop(founderPop, simParam=SP)
  pop<-setPheno(pop,varE=varE,simParam=SP)
  genMean=meanG(pop)
  
  generation<- vector(length=(cycles+1),mode="list")
  generation[[1]]<-pop   ### cycle 0, founder pop
  
  GP0_DH<-makeDH(generation[[1]],nDH=nDH0,simParam=SP)  ### Twice GPs available for making SPs	
  GP0_DH_all<-pullQtlGeno(GP0_DH)	
  dim(GP0_DH_all)
  
  GS<-NULL
  GEBV<-NULL
  
  Sporo<-NULL
  Sporo_s<-NULL
  GP_DH<-NULL
  
  for (j in 1:cycles){ 
    
    if (j<=2){
      ## Year 2, the same scheme as in Year 1
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
      
      ## Cross to create 400 Spj
      Spj<-makeCross(GP0_DH,crossPlan=crossPlan,simParam=SP)
      Spj<-setPheno(Spj,varE=2,simParam=SP)
      Sporo<-c(Sporo,Spj)
      
      ## Select amongst Spj
      Spj_s<-Sel(Sporo[[j]])
      
      Sporo_s<-c(Sporo_s,Spj_s)
      
      ## Make GP DH using Sporo_s
      GP_DHj<-makeDH(Sporo_s[[j]],nDH=nDH,simParam=SP)	### !!!
      GP_DH<-c(GP_DH,GP_DHj)
      print(j)
      
    } else if (j>2) {
      print (j)
      TP_j<-mergePops(Sporo)  
      GS_j<-RRBLUP(TP_j,traits=1,simParam=SP)  ### TP only has info from 1:(j-1)
      
      ## GEBVs on GP(y-2)s
      GEBV_j<-setEBV(GP_DH[[j-2]],GS_j,simParam=SP)  ### Use the GP_DH two years ago
      
      females<-selectInd(pop=GEBV_j,gender="F",nInd=n_gp,selectTop=TRUE,trait=1,use="ebv",simParam=SP)
      males<-selectInd(pop=GEBV_j,gender="M",nInd=n_gp,selectTop=TRUE,trait=1,use="ebv",simParam=SP)
      
      GP_F<-row.names(pullQtlGeno(females))
      GP_M<-row.names(pullQtlGeno(males))
      
      GP_F_2<-GP_F[sample(length(GP_F))]
      GP_M_2<-GP_M[sample(length(GP_M))]
      
      GP_Fs<-c(GP_F,GP_F_2)
      GP_Ms<-c(GP_M,GP_M_2)
      
      crossPlan<-cbind(GP_Fs,GP_Ms)
      
      ## Cross to create Sporo_j
      Sporo_j<-makeCross(GEBV_j,crossPlan=crossPlan,simParam=SP)
      Sporo_j<-setPheno(Sporo_j,varE=2,simParam=SP)
      
      Sporo<-c(Sporo,Sporo_j)
      
      ## Select amongst Sporoy 
      Sporo_js<-Sel(Sporo[[j]])
      
      Sporo_s<-c(Sporo_s,Sporo_js)
      
      ## Make GP_DH_j using Sporo_s
      GP_DH_j<-makeDH(Sporo_s[[j]],nDH=nDH,simParam=SP)
      GP_DH<-c(GP_DH,GP_DH_j)
      
      print(j)   
    }
    
    mean_g1<-unlist(lapply(Sporo,meanG))
    sd_g1<-sqrt(unlist(lapply(Sporo,varG)))
    
  }      
  
  Mean_g_Rep[,i]<-mean_g1
  Sd_g_Rep[,i]<-sd_g1
  
}

Mean_g<-rowMeans(Mean_g_Rep)
Sd_g<-rowMeans(Sd_g_Rep)

Mean_Sd<-cbind(Mean_g,Sd_g)

write.csv(Mean_g_Rep,paste(scheme,"_Mean_g_WithRep.csv",sep=""))
write.csv(Sd_g_Rep,paste(scheme,"_Sd_g_WithRep.csv",sep=""))
write.csv(Mean_Sd,paste(scheme,"_Mean_g_Sd_AverageOverRep.csv",sep=""))

