Simul1Yr<-function(selection,nPheno,nDH,varE,Ne,nrep,cycles){
  
  library(AlphaSimR)
  ## Define: Number of Sporo to Phenotype
  nPheno<-nPheno
  n_gp<-nPheno/2
  
  ## Define: Reps and Cycles
  nrep<-nrep  ## 100 Use shiny to render
  cycles <-cycles ## 7
  
  varE<-varE          # 5.67!!! This is starting H=0.33,  Try another H to represent Dry matter vs Wet biomass (0.15)
  nDH<-nDH          # 25!!! Or 96
  
  Ne=Ne       # Historical effective population size, estimated via E(r2)=1/(1+4*Ne*c)  ## If Ne=600, Pop str will inflat it
  
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
  nDH0<-2           # Number of GPs per SP generated in the founder population to creat the inital GP0_DH
  
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
  GP0_DH_all<-pullQtlGeno(GP0_DH,simParam=SP)	
  dim(GP0_DH_all)
  
  GS<-NULL
  GEBV<-NULL
  
  Sporo<-NULL
  Sporo_s<-NULL
  GP_DH<-NULL
  
  GScor<-NULL
  GPselInt<-NULL
  
  for (j in 1:cycles){ 
    
    if (j<=2){
      ## Year 2, the same scheme as in Year 1
      ## Randomly select 200 as female, then randomly select another 200 as male
      females<-selectInd(pop=GP0_DH,gender="F",nInd=n_gp,trait=1,use="rand",simParam=SP)
      males<-selectInd(pop=GP0_DH,gender="M",nInd=n_gp,trait=1,use="rand",simParam=SP)	
      
      GP_F<-row.names(pullQtlGeno(females,simParam=SP))   #get id of chosen females
      GP_M<-row.names(pullQtlGeno(males,simParam=SP))
      
      GP_F_2<-GP_F[sample(length(GP_F))]	  		#randomize the chosen females, for crossing	 		 
      GP_M_2<-GP_M[sample(length(GP_M))] 					 
      
      GP_Fs<-c(GP_F,GP_F_2)     # every GPs could be used twice for crossing 400 SPs
      GP_Ms<-c(GP_M,GP_M_2)
      
      crossPlan<-cbind(GP_Fs,GP_Ms)  
      
      ## Cross to create 400 Spj
      Spj<-makeCross(GP0_DH,crossPlan=crossPlan,simParam=SP)    # cross selected GP_Fs,GP_Ms
      Spj<-setPheno(Spj,varE=2,simParam=SP)
      Sporo<-c(Sporo,Spj)  ### This is a list now
      
      ## Select amongst Spj
      #Spj_s<-Sel(Sporo[[j]])
      
      Spj_s<-selectInd(Sporo[[j]],nInd=nPheno*0.1,trait=1,use=selection,simParam=SP)  # select Top 10%=40 SPs
      Sporo_s<-c(Sporo_s,Spj_s)
      
      ## Make GP DH using Sporo_s
      GP_DHj<-makeDH(Sporo_s[[j]],nDH=nDH,simParam=SP)	### !!! the number of GPs per SP (nDH) on jth yr; nDH=25/96; so 40SP*25/96 GPs 
      GP_DH<-c(GP_DH,GP_DHj)
      print(j)
      
    } else if (j>2) {
      print (j)
      TP_j<-mergePops(Sporo)  # all previous years SPs go into TP
      GS_j<-RRBLUP(TP_j,traits=1,simParam=SP)  ### TP only has info from 1:(j-1)
      
      ## GEBVs on GP(y-2)s
      GEBV_j<-setEBV(GP_DH[[j-1]],GS_j,simParam=SP)  ### Use the GP_DH 1 years ago, nGP=40(SP) *nDH(GP per SP) 
      
      ## GS cor
      GScorj<-cor(gv(GEBV_j),ebv(GEBV_j))  ## TRUE BVs= genetic values
      GScor<-c(GScor,GScorj)
      
      females<-selectInd(pop=GEBV_j,gender="F",nInd=n_gp,selectTop=TRUE,trait=1,use="ebv",simParam=SP)  # n_gp=nPheno/2
      males<-selectInd(pop=GEBV_j,gender="M",nInd=n_gp,selectTop=TRUE,trait=1,use="ebv",simParam=SP)
      
      ## Selection intensity on FG and MGs
      GEBV_j_df<-data.frame(GEBV_j@id,GEBV_j@ebv,GEBV_j@mother,GEBV_j@father)
      GEBV_j_df$id<-as.character(GEBV_j_df$GEBV_j.id)
      GEBV_j_df$ebv<-GEBV_j_df$GEBV_j.ebv
      GEBV_j_df$mean<-GEBV_j_df$ebv/sd(GEBV_j_df$ebv)
      
      GPFM<-c(females@id,males@id)
      GPselmean<-mean(GEBV_j_df[GEBV_j_df$id%in%GPFM,]$mean)
      
      GPselInt_j<-GPselmean-mean(GEBV_j_df$mean)
      GPselInt<-c(GPselInt,GPselInt_j)
      #dnorm(qnorm(1-0.4))/0.4  # p=400/1000  0.97
      
      GP_F<-row.names(pullQtlGeno(females,simParam=SP))
      GP_M<-row.names(pullQtlGeno(males,simParam=SP))
      
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
      Sporo_js<-selectInd(Sporo[[j]],nInd=nPheno*0.1,trait=1,use=selection,simParam=SP)
      
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

scheme<-paste(selection,"_",nPheno,"_1yr_nDH",nDH,"_varE",varE,"_","Ne",Ne,sep="") ## !!!

write.csv(Mean_g_Rep,paste(scheme,"_Mean_g.csv",sep=""))
write.csv(Sd_g_Rep,paste(scheme,"_Sd_g.csv",sep=""))
write.csv(Mean_Sd,paste(scheme,"_Mean_g_Sd_Average.csv",sep=""))

}