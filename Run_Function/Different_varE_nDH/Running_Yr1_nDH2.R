
source("Simul1Yr_fcn.R")
### nDH2_varE1
Simul1Yr(selection="rand",nPheno=400,nDH=96,varE=5.67,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="rand",nPheno=400,nDH=96,varE=5.67,Ne=600,nrep=100,cycles=7) # Ne change

Simul1Yr(selection="rand",nPheno=1000,nDH=96,varE=5.67,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="rand",nPheno=1000,nDH=96,varE=5.67,Ne=600,nrep=100,cycles=7) # nPheno change

Simul1Yr(selection="pheno",nPheno=400,nDH=96,varE=5.67,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="pheno",nPheno=400,nDH=96,varE=5.67,Ne=600,nrep=100,cycles=7) #selection 

Simul1Yr(selection="pheno",nPheno=1000,nDH=96,varE=5.67,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="pheno",nPheno=1000,nDH=96,varE=5.67,Ne=600,nrep=100,cycles=7) # nPheno 


### nDH2_varE2
Simul1Yr(selection="rand",nPheno=400,nDH=96,varE=2,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="rand",nPheno=400,nDH=96,varE=2,Ne=600,nrep=100,cycles=7) # Ne change

Simul1Yr(selection="rand",nPheno=1000,nDH=96,varE=2,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="rand",nPheno=1000,nDH=96,varE=2,Ne=600,nrep=100,cycles=7) # nPheno change

Simul1Yr(selection="pheno",nPheno=400,nDH=96,varE=2,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="pheno",nPheno=400,nDH=96,varE=2,Ne=600,nrep=100,cycles=7) #selection 

Simul1Yr(selection="pheno",nPheno=1000,nDH=96,varE=2,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="pheno",nPheno=1000,nDH=96,varE=2,Ne=600,nrep=100,cycles=7) # nPheno 

