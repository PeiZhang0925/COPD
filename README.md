# COPD   
     
This repository contains all the codes and data used to produce the results of the HK-COPD paper.    
     
## Analyses     
     
The Analyses folder contains all the codes used in this study.      
     
(1) PCMCI+.py        
(2) DYNOTEARS.py   
(3) GAM.R     
      
### 1. PCMCI+.py       
     
**PCMCI+.py** file is used for the PCMCI+ analysis.     
     
PCMCI+ belongs to the so-called constraint-based causal discovery methods family, which is based on conditional independence test. Here “PC” refers to the developers Peter and Clark, “MCI” means that the momentary conditional independence (MCI) test idea is added to the traditional PC algorithm, and “+” reminds users that it extends the earlier version of PCMCI to include the discovery of contemporaneous links. Like other causal graphic models, PCMCI+ works under the general assumptions of the causal markov condition and faithfulness. On top of the general assumptions, two specific assumptions are also requested: causal stationarity (i.e., the causal links hold for all the studied time points) and causal sufficiency (i.e. measured variables include all of the common causes).      
     
### 2. DYNOTEARS.py     
  
**DYNOTEARS.py** file is for the analyses of DYNOTEARS.    
     
Deriving from NOTEARS, DYNOTEARS is a score-based causal discovery algorithm which leverages Structural Equation Model (SEM) to illustrate causal relationships among time series variables. It converts the combinatorial optimization problem in conventional score-based methods into a continuous programming task, thereby enabling facile resolution through standard optimization methods and obtaining better computation efficiency.          
     
### 3. GAM.R          
        
**GAM.R** contains all the codes for GAM analysis.     
   
GAM is the expansion of the generalized linear model, aiming at modeling non-linear interactions. In environmental health domain, GAM has dominated time series studies for decades. To explore the associations between environmental factors and COPD exacerbation, a GAM with Quasi-Poisson distribution was adopted in this study.         

