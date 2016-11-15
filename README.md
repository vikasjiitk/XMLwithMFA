Extreme Multi-label Learning Using Mixture of Factor Analyzers
=======
Undergraduate Project, Autumn 2016  
[Vikas Jain](http://home.iitk.ac.in/~vikasj)  
Under the guidance of [Prof Piyush Rai](http://cse.iitk.ac.in/users/piyush)  

***
##Project Idea 
To develop a probabilistic model using mixture of factor analyzers for extreme multi-label learning. The MFA model used is adapted from the work by Montrari et al. "Dimensionally reduced mixtures of regression models" [(paper link)](http://www.sciencedirect.com/science/article/pii/S0378375810005240).  
The purpose of the project is to scale the learning of the model and incorporating 0-1 labelling instead of continuous values. To achieve this, online EM algorithm and EM algorithm for binary logistic regression are used based on the work by James G. Scott et al. "Expectation-maximization for logistic regression" [(paper link)](https://arxiv.org/abs/1306.0040).

##Code Source
From [https://cran.r-project.org/web/packages/FactMixtAnalysis/](https://cran.r-project.org/web/packages/FactMixtAnalysis/)

##Modifications
####Online EM
Added a new function `onlinefma()` for online EM algorithm.  
How to use:
```
fit = onlinefma(y,k,r,x.z=NULL,x.w=NULL,it=1,eps=0.0001,seed=4,scaling=FALSE,init=NULL,no_iter=100,batch_size=1000)
```
Refer to the document [here](https://cran.r-project.org/web/packages/FactMixtAnalysis/) for details of the variables.

#### EM algorithm for binary logistic regression
TODO

##Experiments
TODO
