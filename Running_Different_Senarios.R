#### varE1
varE=5.67;nDH=25;nPheno=400;nrep=1;selection="rand";cycles=3


source("Scheme2Yr_fcn.R")
source("Scheme1Yr_fcn.R")


Simul(varE=5.67,nDH=25,nPheno=400,nrep=1,selection="rand",cycles=3)



varE=5.67;nDH=25;nPheno=1000;nrep=100;s="Rand";cycles=5
source("Scheme2Yr_fcn.R")

varE=5.67;nDH=25;nPheno=400;nrep=100;s="Top";cycles=5
source("Scheme2Yr_fcn.R")

varE=5.67;nDH=25;nPheno=1000;nrep=100;s="Top";cycles=5
source("Scheme2Yr_fcn.R")

### Change varE2
varE=2;nDH=25;nPheno=400;nrep=1;s="Rand";cycles=3
source("Scheme2Yr_fcn.R")

varE=2;nDH=25;nPheno=1000;nrep=100;s="Rand";cycles=5
source("Scheme2Yr_fcn.R")

varE=2;nDH=25;nPheno=400;nrep=100;s="Top";cycles=5
source("Scheme2Yr_fcn.R")

varE=2;nDH=25;nPheno=1000;nrep=100;s="Top";cycles=5
source("Scheme2Yr_fcn.R")

