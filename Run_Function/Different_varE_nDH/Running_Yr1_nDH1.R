#### Test

#Simul1Yr(selection="rand",nPheno=400,nDH=25,varE=5.67,Ne=60,nrep=10,cycles=5)
#Simul2Yr(selection="rand",nPheno=400,nDH=25,varE=5.67,Ne=60,nrep=10,cycles=5)

source("Simul1Yr_fnc.R")
### nDH1_varE1    Var_Trait=1, H=0.15=1/(1+Ve) -> Ve=5.67
Simul1Yr(selection="rand",nPheno=400,nDH=25,varE=5.67,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="rand",nPheno=400,nDH=25,varE=5.67,Ne=600,nrep=100,cycles=7) # Ne change

Simul1Yr(selection="rand",nPheno=1000,nDH=25,varE=5.67,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="rand",nPheno=1000,nDH=25,varE=5.67,Ne=600,nrep=100,cycles=7) # nPheno change

Simul1Yr(selection="pheno",nPheno=400,nDH=25,varE=5.67,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="pheno",nPheno=400,nDH=25,varE=5.67,Ne=600,nrep=100,cycles=7) #selection 

Simul1Yr(selection="pheno",nPheno=1000,nDH=25,varE=5.67,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="pheno",nPheno=1000,nDH=25,varE=5.67,Ne=600,nrep=100,cycles=7) # nPheno 


### nDH1_varE2  Var_Trait=1, H=0.33=1/(1+Ve) -> Ve=2
Simul1Yr(selection="rand",nPheno=400,nDH=25,varE=2,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="rand",nPheno=400,nDH=25,varE=2,Ne=600,nrep=100,cycles=7)

Simul1Yr(selection="rand",nPheno=1000,nDH=25,varE=2,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="rand",nPheno=1000,nDH=25,varE=2,Ne=600,nrep=100,cycles=7)

Simul1Yr(selection="pheno",nPheno=400,nDH=25,varE=2,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="pheno",nPheno=400,nDH=25,varE=2,Ne=600,nrep=100,cycles=7)

Simul1Yr(selection="pheno",nPheno=1000,nDH=25,varE=2,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="pheno",nPheno=1000,nDH=25,varE=2,Ne=600,nrep=100,cycles=7)

### nDH1_varE3: Var_Trait=1, H=0.45=1/(1+Ve) -> Ve=1.22
Simul1Yr(selection="rand",nPheno=400,nDH=25,varE=1.22,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="rand",nPheno=400,nDH=25,varE=1.22,Ne=600,nrep=100,cycles=7)

Simul1Yr(selection="rand",nPheno=1000,nDH=25,varE=1.22,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="rand",nPheno=1000,nDH=25,varE=1.22,Ne=600,nrep=100,cycles=7)

Simul1Yr(selection="pheno",nPheno=400,nDH=25,varE=1.22,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="pheno",nPheno=400,nDH=25,varE=1.22,Ne=600,nrep=100,cycles=7)

Simul1Yr(selection="pheno",nPheno=1000,nDH=25,varE=1.22,Ne=60,nrep=100,cycles=7)
Simul1Yr(selection="pheno",nPheno=1000,nDH=25,varE=1.22,Ne=600,nrep=100,cycles=7)
