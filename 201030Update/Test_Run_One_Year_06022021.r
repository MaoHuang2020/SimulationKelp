library(parallel)


cycles <- 10
nrep<-1

runOneRep<-function(selection,nPheno,nDH,varE,Ne){

  library(AlphaSimR)
  n_gp<-nPheno/2
  founderPop<-runMacs2(nInd=nInd,nChr=nChr,segSites=segSites,Ne=Ne,bp=bp,genLen=1,inbred=TRUE,ploidy=Sporo_ploidy, returnCommand = FALSE, nThreads = NULL)
  # Founder Pop should be here
  
  SP<-SimParam$new(founderPop)
  SP$addTraitA(nQtlPerChr=nQtlPerChr, mean=mean_Trait,var=var_Trait,corA=NULL)
  SP$addSnpChip(nSnpPerChr=nSnpPerChr)
  SP$setSexes("yes_sys")
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
               SPselInt<-NULL
               for (j in 1:cycles){
                 
                 if (j<=2){
                   ## Year 2, the same scheme as in Year 1
                   ## Randomly select 200 as female, then randomly select another 200 as male
                   females<-selectInd(pop=GP0_DH,sex="F",nInd=n_gp,trait=1,use="rand",simParam=SP)
                   males<-selectInd(pop=GP0_DH,sex="M",nInd=n_gp,trait=1,use="rand",simParam=SP)
                   
                   GP_F<-row.names(pullQtlGeno(females,simParam=SP))
                   GP_M<-row.names(pullQtlGeno(males,simParam=SP))
                   
                   GP_F_2<-GP_F[sample(length(GP_F))]
                   GP_M_2<-GP_M[sample(length(GP_M))]
                   
                   GP_Fs<-c(GP_F,GP_F_2)
                   GP_Ms<-c(GP_M,GP_M_2)
                   
                   crossPlan<-cbind(GP_Fs,GP_Ms)
                   
                   ## Cross to create 400 Spj
                   Spj<-makeCross(GP0_DH,crossPlan=crossPlan,simParam=SP)
                   Spj<-setPheno(Spj,varE=varE,simParam=SP)
                   Sporo<-c(Sporo,Spj)  ### This is a list now
                   
                   ## Select amongst Spj
                   #Spj_s<-Sel(Sporo[[j]])
                   
                   Spj_s<-selectInd(Sporo[[j]],nInd=nPheno*0.1,trait=1,use=selection,simParam=SP)
                   Sporo_s<-c(Sporo_s,Spj_s)
                   
                   ## Selection Intensity on Sporo[[j]]
                   selectpop<-Spj_s@id
                   Refpop<-data.frame(Sporo[[j]]@id,Sporo[[j]]@pheno,Sporo[[j]]@mother,Sporo[[j]]@father)  # select on pheno
                   
                   ## Pull out selectInd, estimate the selection intensity 
                   colnames(Refpop)[1:2]<-c("id","trait")
                   Refpop$id<-as.character(Refpop$id)
                   Refpop$mean<-Refpop$trait/sd(Refpop$trait)
                   selectPopmean<-mean(Refpop[Refpop$id%in%selectpop,]$mean)
                   selectInt<-selectPopmean-mean(Refpop$mean)
                   
                   SPselInt<-c(SPselInt,selectInt) 
                   
                   ## Make GP DH using Sporo_s
                   GP_DHj<-makeDH(Sporo_s[[j]],nDH=nDH,simParam=SP)    ### !!!
                   GP_DH<-c(GP_DH,GP_DHj)
                   print(j)
                   
                 } else if (j>2) {
                   print (j)
                   TP_j<-mergePops(Sporo)
                   GS_j<-RRBLUP(TP_j,traits=1,simParam=SP)  ### TP only has info from 1:(j-1)
                   
                   ## GEBVs on GP(y-2)s
                   GEBV_j<-setEBV(GP_DH[[j-1]],GS_j,simParam=SP)  ### Use the GP_DH 1 years ago
                   
                   ## GS cor. TRUE BVs= genetic values
                   GScor<-c(GScor,cor(gv(GEBV_j),ebv(GEBV_j)))
                   
                   females<-selectInd(pop=GEBV_j,sex="F",nInd=n_gp,selectTop=TRUE,trait=1,use="ebv",simParam=SP)
                   males<-selectInd(pop=GEBV_j,sex="M",nInd=n_gp,selectTop=TRUE,trait=1,use="ebv",simParam=SP)
                   
                   ## Selection intensity on FG and MGs
                     selectpop<-c(females@id,males@id)
                     
                     Refpop<-data.frame(GEBV_j@id,GEBV_j@ebv,GEBV_j@mother,GEBV_j@father)
                     ## Pull out selectInd, estimate the selection intensity 
                     colnames(Refpop)[1:2]<-c("id","trait")
                     Refpop$mean<-Refpop$trait/sd(Refpop$trait)
                     selectPopmean<-mean(Refpop[Refpop$id%in%selectpop,]$mean)
                     selectInt<-selectPopmean-mean(Refpop$mean)
                   GPselInt<-c(GPselInt,selectInt) 
                   
                   ## another way pull out list of females and males
                   GP_F<-row.names(pullQtlGeno(females,simParam=SP))
                   GP_M<-row.names(pullQtlGeno(males,simParam=SP))
                   
                   GP_F_2<-GP_F[sample(length(GP_F))]
                   GP_M_2<-GP_M[sample(length(GP_M))]
                   
                   GP_Fs<-c(GP_F,GP_F_2)
                   GP_Ms<-c(GP_M,GP_M_2)
                   
                   crossPlan<-cbind(GP_Fs,GP_Ms)
                   
                   ## Cross to create Sporo_j
                   Sporo_j<-makeCross(GEBV_j,crossPlan=crossPlan,simParam=SP)
                   Sporo_j<-setPheno(Sporo_j,varE=varE,simParam=SP)
                   
                   Sporo<-c(Sporo,Sporo_j)
                   
                   ## Select amongst Sporoy from jth year
                   Sporo_js<-selectInd(Sporo[[j]],nInd=nPheno*0.1,trait=1,use=selection,simParam=SP)  # select 10% SPs
                   
                   Sporo_s<-c(Sporo_s,Sporo_js)
                   
                   ## Selection intensity on Sporo[[j]]
                   selectpop<-Sporo_js@id
                   Refpop<-data.frame(Sporo[[j]]@id,Sporo[[j]]@pheno,Sporo[[j]]@mother,Sporo[[j]]@father)  # select on pheno
                  
                   ## Pull out selectInd, estimate the selection intensity 
                   colnames(Refpop)[1:2]<-c("id","trait")
                   Refpop$id<-as.character(Refpop$id)
                   Refpop$mean<-Refpop$trait/sd(Refpop$trait)
                   selectPopmean<-mean(Refpop[Refpop$id%in%selectpop,]$mean)
                   selectInt<-selectPopmean-mean(Refpop$mean)
                   
                  SPselInt<-c(SPselInt,selectInt) 
                   
                   ## Make GP_DH_j using Sporo_s
                   GP_DH_j<-makeDH(Sporo_s[[j]],nDH=nDH,simParam=SP)
                   GP_DH<-c(GP_DH,GP_DH_j)
                   
                   print(j)
                 }
                 mean_g1<-unlist(lapply(Sporo,meanG))
                 sd_g1<-sqrt(unlist(lapply(Sporo,varG)))
                 
                 mean_g2<-unlist(lapply(GP_DH,meanG))
                 sd_g2<-sqrt(unlist(lapply(GP_DH,varG)))
                 
               }    
               return(list(mean_g1=mean_g1, sd_g1=sd_g1,mean_g2=mean_g2,sd_g2=sd_g2,GPselInt=GPselInt,SPselInt=SPselInt)) 
}   
  
#END runOneRep


for (selection in c("rand","pheno")){
  for (nPheno in c(400,1000)){
    for (nDH in c(24,96)){
      for (varE in c(1,4)){
        for (Ne in c(60,600)){
          
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
 
         allRep <- mclapply(1:nrep, function(dummy) {print(dummy); runOneRep(selection, nPheno, nDH, varE, Ne)},mc.cores=15)
          
         Mean_SP_Rep<-matrix(nrow=cycles,ncol=nrep)
         Sd_SP_Rep<-matrix(nrow=cycles,ncol=nrep)
         
         Mean_GP_Rep<-matrix(nrow=cycles,ncol=nrep)
         Sd_GP_Rep<-matrix(nrow=cycles,ncol=nrep)
        
         GPselInt_Rep<-matrix(nrow=cycles-2,ncol=nrep)  
         SPselInt_Rep<-matrix(nrow=cycles,ncol=nrep)  
         
          for (i in 1:nrep){
            Mean_SP_Rep[,i] <- allRep[[i]]$mean_g1
            Sd_SP_Rep[,i] <- allRep[[i]]$sd_g1
            
            Mean_GP_Rep[,i]<-allRep[[i]]$mean_g2
            Sd_GP_Rep[,i]<-allRep[[i]]$sd_g2
            
            GPselInt_Rep[,i]<-allRep[[i]]$GPselInt
            SPselInt_Rep[,i]<-allRep[[i]]$SPselInt
          }
          
          Mean_SP<-rowMeans(Mean_SP_Rep)
          Sd_SP<-rowMeans(Sd_SP_Rep)
          Mean_Sd_SP<-cbind(Mean_SP,Sd_SP)
          
          Mean_GP<-rowMeans(Mean_GP_Rep)
          Sd_GP<-rowMeans(Sd_GP_Rep)
          Mean_Sd_GP<-cbind(Mean_GP,Sd_GP)
          
          GPselInt_Rep<-rowMeans(GPselInt_Rep)
          SPselInt_Rep<-rowMeans(SPselInt_Rep)
          
          scheme<-paste(selection,"_",nPheno,"_1yr_nDH",nDH,"_varE",varE,"_","Ne",Ne,sep="") ## !!!
          
          write.csv(Mean_SP_Rep,paste(scheme,"_Mean_SP.csv",sep=""))
          write.csv(Sd_SP_Rep,paste(scheme,"_Sd_SP.csv",sep=""))
          write.csv(Mean_Sd_SP,paste(scheme,"_Mean_g_Sd_Average_SP.csv",sep=""))
          
          write.csv(Mean_GP_Rep,paste(scheme,"_Mean_GP.csv",sep=""))
          write.csv(Sd_GP_Rep,paste(scheme,"_Sd_GP.csv",sep=""))
          write.csv(Mean_Sd_GP,paste(scheme,"_Mean_g_Sd_Average_GP.csv",sep=""))
          
          write.csv(GPselInt_Rep,paste(scheme,"_GPselInt.csv",sep=""))
          write.csv(SPselInt_Rep,paste(scheme,"_SPselInt.csv",sep=""))
          
        }
      }
    }
  }
}

