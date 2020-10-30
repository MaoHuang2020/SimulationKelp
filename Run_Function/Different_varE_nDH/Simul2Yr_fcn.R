Simul2Yr<-function(selection,nPheno,nDH,varE,Ne,nrep,cycles){
  
  library(AlphaSimR)
  ## Define: Number of Sporo to Phenotype
  nPheno<-nPheno
  n_gp<-nPheno/2
  
  ## Define: Reps and Cycles
  nrep<-nrep  ## 100 Use shiny to render
  cycles <-cycles ## 5
  
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
  
  ## Define: Sporo selection *Random* vs *Top*

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
    #GP0_DH_all<-pullQtlGeno(GP0_DH,simParam=SP)	### If not specify simParam=SP, then it is searching "SP" from global envir, which is isolated from the function envir and which is not defined and so it gives error "cannot find object SP"  
    #dim(GP0_DH_all) ## Use debugonce(Simul) to check out what is all going on
    
    GS<-NULL
    GEBV<-NULL
    
    Sporo<-NULL
    Sporo_s<-NULL
    GP_DH<-NULL
    
    Gameto_F<-NULL
    Gameto_M<-NULL
    
    for (j in 1:cycles){ 
      
      if (j<=2){
        ## Year 2, the same scheme as in Year 1
        ## Randomly select 200 as female, then randomly select another 200 as male
        females<-selectInd(pop=GP0_DH,gender="F",nInd=n_gp,trait=1,use="rand",simParam=SP)
        males<-selectInd(pop=GP0_DH,gender="M",nInd=n_gp,trait=1,use="rand",simParam=SP)	
        
        GP_F<-row.names(pullQtlGeno(females,simParam=SP))
        GP_M<-row.names(pullQtlGeno(males,simParam=SP))
        
        GP_F_2<-GP_F[sample(length(GP_F))]	  					 
        GP_M_2<-GP_M[sample(length(GP_M))] 					 
        
        GP_Fs<-c(GP_F,GP_F_2)  # make crosses for two times
        GP_Ms<-c(GP_M,GP_M_2)
        
        crossPlan<-cbind(GP_Fs,GP_Ms)
        
        Gameto_F<-c(Gameto_F,GP_Fs)  # store the GP_Fs
        Gameto_M<-c(Gameto_M,GP_Ms)  # store the GP_Ms
        Gameto_Both<-c(Gameto_F,Gameto_M)
          
        ## Cross to create 400 Spj
        Spj<-makeCross(GP0_DH,crossPlan=crossPlan,simParam=SP)
        Spj<-setPheno(Spj,varE=2,simParam=SP)
        Sporo<-c(Sporo,Spj)
        
        ## Select amongst Spj
        
        Spj_s<-selectInd(Sporo[[j]],nInd=nPheno*0.1,trait=1,use=selection,simParam=SP)
        
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
        
        GP_F<-row.names(pullQtlGeno(females,simParam=SP))
        GP_M<-row.names(pullQtlGeno(males,simParam=SP))
        
        GP_F_2<-GP_F[sample(length(GP_F))]
        GP_M_2<-GP_M[sample(length(GP_M))]
        
        GP_Fs<-c(GP_F,GP_F_2)
        GP_Ms<-c(GP_M,GP_M_2)
        
        crossPlan<-cbind(GP_Fs,GP_Ms)
        
        Gameto_F<-c(Gameto_F,GP_Fs) # store GP_Fs
        Gameto_M<-c(Gameto_M,GP_Ms) # store GP_Ms
        Gameto_Both<-c(Gameto_Both,Gameto_F,Gameto_M)
        
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
      
      mean_g1<-unlist(lapply(Sporo,meanG))   #Sporo contains the list of Sporophytes from each year
      sd_g1<-sqrt(unlist(lapply(Sporo,varG)))  
      
      
      mean_g2F<-unlist(lapply(Gameto_F,meanG))
      sd_g2F<-sqrt(unlist(lapply(Gemeto_F,varG)))
      
      mean_g2M<-unlist(lapply(Gameto_M,meanG))
      sd_g2M<-sqrt(unlist(lapply(Gemeto_M,varG)))
      
      mean_g2Both<-unlist(lapply(Gameto_Both,meanG))
      sd_g2Both<-sqrt(unlist(lapply(Gemeto_Both,varG)))
      
    }      
    
    Mean_g_Rep1[,i]<-mean_g1
    Sd_g_Rep1[,i]<-sd_g1
    
    Mean_g_Rep2F[,i]<-mean_g2F
    Sd_g_Rep2F[,i]<-sd_g2F
    
    Mean_g_Rep2M[,i]<-mean_g2M
    Sd_g_Rep2M[,i]<-sd_g2M
    
    Mean_g_Rep2Both[,i]<-mean_g2Both
    Sd_g_Rep2Both[,i]<-sd_g2Both
    
  }
  
  Mean_g1<-rowMeans(Mean_g_Rep1)
  Sd_g1<-rowMeans(Sd_g_Rep1)
  Mean_Sd1<-cbind(Mean_g1,Sd_g1)

  Mean_g2F<-rowMeans(Mean_g_Rep2F)
  Sd_g2F<-rowMeans(Sd_g_Rep2F)
  Mean_Sd2F<-cbind(Mean_g2F,Sd_g2F)
  
  Mean_g2M<-rowMeans(Mean_g_Rep2M)
  Sd_g2M<-rowMeans(Sd_g_Rep2M)
  Mean_Sd2M<-cbind(Mean_g2M,Sd_g2M)
  
  Mean_g2Both<-rowMeans(Mean_g_Rep2Both)
  Sd_g2Both<-rowMeans(Sd_g_Rep2Both)
  Mean_Sd2Both<-cbind(Mean_g2Both,Sd_g2Both)
  
  scheme<-paste(selection,"_",nPheno,"_2yr_nDH",nDH,"_varE",varE,"_","Ne",Ne,sep="") ## !!!
  
  write.csv(Mean_g_Rep1,paste(scheme,"_Mean_g.csv",sep=""))
  write.csv(Sd_g_Rep1,paste(scheme,"_Sd_g.csv",sep=""))
  write.csv(Mean_Sd1,paste(scheme,"_Mean_g_Sd_Average.csv",sep=""))
  
  write.csv(Mean_g_Rep2F,paste(scheme,"_Mean_g_GPsF.csv",sep=""))
  write.csv(Sd_g_Rep2F,paste(scheme,"_Sd_g_GPsF.csv",sep=""))
  write.csv(Mean_Sd2F,paste(scheme,"_Mean_g_Sd_Average_GPsF.csv",sep=""))
  
  write.csv(Mean_g_Rep2M,paste(scheme,"_Mean_g_GPsM.csv",sep=""))
  write.csv(Sd_g_Rep2M,paste(scheme,"_Sd_g_GPsM.csv",sep=""))
  write.csv(Mean_Sd2M,paste(scheme,"_Mean_g_Sd_Average_GPsM.csv",sep=""))
  
  write.csv(Mean_g_Rep2F,paste(scheme,"_Mean_g_GPsF.csv",sep=""))
  write.csv(Sd_g_Rep2F,paste(scheme,"_Sd_g_GPsF.csv",sep=""))
  write.csv(Mean_Sd2F,paste(scheme,"_Mean_g_Sd_Average_GPsF.csv",sep=""))
  
}