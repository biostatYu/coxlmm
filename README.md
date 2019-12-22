# coxlmm: linear mixed effect Cox model with clinical covariates and gene expressions levels

## Introduction
linear mixed effect Cox model with clinical covariates and gene expressions levels (coxlmm) is a [**R**](https://cran.r-project.org/) procedure to access the prognosis of different cancers and estimate the partition of clinical variance and genetic variance.

Specifically, denote the observed survival time by *t<sub>i</sub>* and the true survival time by *T<sub>i</sub>* with di indicating the censored status, *X<sub>i</sub>* is a *p*-dimensional vector for available clinical covariates (e.g. disease stage, age and gender) for individual *i*, and *G<sub>i</sub>* be an *m*-dimensional vector for a set of gene expression levels for individual. To link the survival risk with the clinical information and all the genetic information, we employ Cox model within the framework of linear mixed models:
<img src="https://latex.codecogs.com/gif.latex?h({t_i}|{X_i},\;{G_i})\;&space;=&space;\;{h_0}({t_i}){e^{X_i^Ta\;&space;&plus;&space;\;G_i^Tb}},\;{\rm{&space;}}{b_j}\;\~\;N(0,\;{\rm{\sigma&space;}}_b^2)\" title="h({t_i}|{X_i},\;{G_i})\; = \;{h_0}({t_i}){e^{X_i^Ta\; + \;G_i^Tb}},\;{\rm{ }}{b_j}\;\~\;N(0,\;{\rm{\sigma }}_b^2)\" />


***` y = Xa + Gb + e; b ~ N(0, sigam_g2), e ~　N(0, sigam_e2)　`***

Above, sigam_g2 is the genetic variance, and sigam_e2 is residual variance. 

LRT examines the association of ***G*** with ***y*** (while controling for ***X***) by testing for:

***` H0: b = 0; <==>, H0: sigam_g2 = 0 `***

The LRT statistic is defined as the difference of log-likelihood value between the null model (i.e. y = Xa + e) and the alternative model (i.e. y = Xa + Gb + e). Due to the limited parameter space of sigam_g2, obtaining the null distribution of the LRT statistic is not trivial. We use the method of spectral decomposition proposed by Crainiceanu & Ruppert [(2004)](http://onlinelibrary.wiley.com/wol1/doi/10.1111/j.1467-9868.2004.00438.x/abstract) to generate the exact finite sample null distribution (implemented via the [**R**](https://cran.r-project.org/) package [RLRsim](https://github.com/fabian-s/RLRsim)).


## Note
The LRT procedure was finished in about 2013. At that time, the author (i.e., [Ping Zeng](https://github.com/biostatpzeng)) was a newer and was outside the door of statistical genetics. Thus, you will find that this procedure was not well designed. I put it here for a beautiful recall of that time.

#### A new LRT R function, named [ReLRT](https://github.com/biostatpzeng/LRT/blob/master/ReLRT.R), was recently rewritten. We applied [ReLRT](https://github.com/biostatpzeng/LRT/blob/master/ReLRT.R) to detect eGenes in gene expression data using cis-SNPs.  [ReLRT](https://github.com/biostatpzeng/LRT/blob/master/ReLRT.R) esimates the linear mixed models basde on lme function in R package [nlme](https://cran.r-project.org/web/packages/nlme/index.html) and thus is more efficient than pervious [LRT](https://github.com/biostatpzeng/LRT/blob/master/LRT.R). [ReLRT](https://github.com/biostatpzeng/LRT/blob/master/ReLRT.R) also perfroms the approximate LRT (aLRT) using a mixture null distribution. We call the original LRT via the simulation-based algrithmn exact LRT (eLRT). 

## References
+ [Ping Zeng](https://github.com/biostatpzeng), Yang Zhao, Jin Liu, Liya Liu, Liwei Zhang, Ting Wang, Shuiping Huang and Feng Chen. Likelihood Ratio Tests in Rare Variant Detection for Continuous Phenotypes. Annals of Human Genetics, 2014, 78(5): 320-332. [DOI: 10.1111/ahg.12071](http://onlinelibrary.wiley.com/wol1/doi/10.1111/ahg.12071/abstract) 

+ [**Ping Zeng**](https://github.com/biostatpzeng), Ting Wang and Shuiping Huang (2017). Cis-SNPs Set Testing and PrediXcan Analysis for Gene Expression Data using Linear Mixed Model. **Scientific Reports** (in press).

+ Ciprian M. Crainiceanu and David Ruppert. Likelihood ratio tests in linear mixed models with one variance component. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 2004, 66(1): 165–185. [DOI: 10.1111/j.1467-9868.2004.00438.x](http://onlinelibrary.wiley.com/wol1/doi/10.1111/j.1467-9868.2004.00438.x/abstract) 

+ Fabian Scheipl, Sonja Greven and Helmut Küchenhoff. Size and power of tests for a zero random effect variance or polynomial regression in additive and linear mixed models. Computational Statistics & Data Analysis, 2008, 52(7): 3283-3299. [DOI: 10.1016/j.csda.2007.10.022](http://www.sciencedirect.com/science/article/pii/S0167947307004306)

+ R Core Team. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria, 2013 (https://www.R-project.org/). 

## Cite
#### [Ping Zeng](https://github.com/biostatpzeng), Yang Zhao, Jin Liu, Liya Liu, Liwei Zhang, Ting Wang, Shuiping Huang and Feng Chen. Likelihood Ratio Tests in Rare Variant Detection for Continuous Phenotypes. Annals of Human Genetics, 2014, 78(5): 320-332. [DOI: 10.1111/ahg.12071](http://onlinelibrary.wiley.com/wol1/doi/10.1111/ahg.12071/abstract) 


## Contact
I am very grateful to any questions, comments, or bugs reports; and please contact [Ping Zeng](https://github.com/biostatpzeng) via zpstat@xzhmu.edu.cn or ~~pingzeng@umich.edu~~.
