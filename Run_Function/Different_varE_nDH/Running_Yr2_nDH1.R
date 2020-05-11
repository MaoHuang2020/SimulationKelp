#### Test

#Simul2Yr(selection="rand",nPheno=400,nDH=25,varE=5.67,Ne=60,nrep=10,cycles=5)
#Simul2Yr(selection="rand",nPheno=400,nDH=25,varE=5.67,Ne=60,nrep=10,cycles=5)

source("Simul2Yr_fnc.R")
### nDH1_varE1
Simul2Yr(selection="rand",nPheno=400,nDH=25,varE=5.67,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="rand",nPheno=400,nDH=25,varE=5.67,Ne=600,nrep=100,cycles=7) # Ne change

Simul2Yr(selection="rand",nPheno=1000,nDH=25,varE=5.67,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="rand",nPheno=1000,nDH=25,varE=5.67,Ne=600,nrep=100,cycles=7) # nPheno change

Simul2Yr(selection="pheno",nPheno=400,nDH=25,varE=5.67,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="pheno",nPheno=400,nDH=25,varE=5.67,Ne=600,nrep=100,cycles=7) #selection 

Simul2Yr(selection="pheno",nPheno=1000,nDH=25,varE=5.67,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="pheno",nPheno=1000,nDH=25,varE=5.67,Ne=600,nrep=100,cycles=7) # nPheno 


### nDH1_varE2
Simul2Yr(selection="rand",nPheno=400,nDH=25,varE=2,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="rand",nPheno=400,nDH=25,varE=2,Ne=600,nrep=100,cycles=7)

Simul2Yr(selection="rand",nPheno=1000,nDH=25,varE=2,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="rand",nPheno=1000,nDH=25,varE=2,Ne=600,nrep=100,cycles=7)

Simul2Yr(selection="pheno",nPheno=400,nDH=25,varE=2,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="pheno",nPheno=400,nDH=25,varE=2,Ne=600,nrep=100,cycles=7)

Simul2Yr(selection="pheno",nPheno=1000,nDH=25,varE=2,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="pheno",nPheno=1000,nDH=25,varE=2,Ne=600,nrep=100,cycles=7)


### nDH1_varE3
Simul2Yr(selection="rand",nPheno=400,nDH=25,varE=1.22,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="rand",nPheno=400,nDH=25,varE=1.22,Ne=600,nrep=100,cycles=7)

Simul2Yr(selection="rand",nPheno=1000,nDH=25,varE=1.22,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="rand",nPheno=1000,nDH=25,varE=1.22,Ne=600,nrep=100,cycles=7)

Simul2Yr(selection="pheno",nPheno=400,nDH=25,varE=1.22,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="pheno",nPheno=400,nDH=25,varE=1.22,Ne=600,nrep=100,cycles=7)

Simul2Yr(selection="pheno",nPheno=1000,nDH=25,varE=1.22,Ne=60,nrep=100,cycles=7)
Simul2Yr(selection="pheno",nPheno=1000,nDH=25,varE=1.22,Ne=600,nrep=100,cycles=7)
