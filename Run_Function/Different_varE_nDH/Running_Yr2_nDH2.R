
source("Simul2Yr_fcn.R")
### nDH2_varE1
Simul2Yr(selection="rand",nPheno=400,nDH=96,varE=5.67,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="rand",nPheno=400,nDH=96,varE=5.67,Ne=600,nrep=100,cycles=7) # Ne change

Simul2Yr(selection="rand",nPheno=1000,nDH=96,varE=5.67,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="rand",nPheno=1000,nDH=96,varE=5.67,Ne=600,nrep=100,cycles=7) # nPheno change

Simul2Yr(selection="pheno",nPheno=400,nDH=96,varE=5.67,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="pheno",nPheno=400,nDH=96,varE=5.67,Ne=600,nrep=100,cycles=7) #selection

Simul2Yr(selection="pheno",nPheno=1000,nDH=96,varE=5.67,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="pheno",nPheno=1000,nDH=96,varE=5.67,Ne=600,nrep=100,cycles=7) # nPheno


### nDH2_varE2
Simul2Yr(selection="rand",nPheno=400,nDH=96,varE=2,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="rand",nPheno=400,nDH=96,varE=2,Ne=600,nrep=100,cycles=7) # Ne change

Simul2Yr(selection="rand",nPheno=1000,nDH=96,varE=2,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="rand",nPheno=1000,nDH=96,varE=2,Ne=600,nrep=100,cycles=7) # nPheno change

Simul2Yr(selection="pheno",nPheno=400,nDH=96,varE=2,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="pheno",nPheno=400,nDH=96,varE=2,Ne=600,nrep=100,cycles=7) #selection

Simul2Yr(selection="pheno",nPheno=1000,nDH=96,varE=2,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="pheno",nPheno=1000,nDH=96,varE=2,Ne=600,nrep=100,cycles=7) # nPheno

