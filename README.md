# coxlmm: linear mixed effect Cox model with clinical covariates and gene expressions levels

## Introduction
linear mixed effect Cox model with clinical covariates and gene expressions levels (coxlmm) is a [**R**](https://cran.r-project.org/) procedure to access the prognosis of different cancers and estimate the partition of clinical variance and genetic variance.

Specifically, denote the observed survival time by *t<sub>i</sub>* and the true survival time by *T<sub>i</sub>* with di indicating the censored status, *X<sub>i</sub>* is a *p*-dimensional vector for available clinical covariates (e.g. disease stage, age and gender) for individual *i*, and *G<sub>i</sub>* be an *m*-dimensional vector for a set of gene expression levels for individual. To link the survival risk with the clinical information and all the genetic information, we employ Cox model within the framework of linear mixed models:

<div align=center><img src="https://latex.codecogs.com/gif.latex?h({t_i}|{X_i},\;{G_i})\;&space;=&space;\;{h_0}({t_i}){e^{X_i^Ta\;&space;&plus;&space;\;G_i^Tb}},\;{\rm{&space;}}{b_j}\;\~\;N(0,\;{\rm{\sigma&space;}}_b^2)" title="h({t_i}|{X_i},\;{G_i})\; = \;{h_0}({t_i}){e^{X_i^Ta\; + \;G_i^Tb}},\;{\rm{ }}{b_j}\;\~\;N(0,\;{\rm{\sigma }}_b^2)" />
</div>

where <img src="https://latex.codecogs.com/gif.latex?{\rm{\sigma&space;}}_b^2" title="{\rm{\sigma }}_b^2" /> is the variance for gene expressions 
The coxlmm model (6) is fitted with the R coxme (version 2.2-10) package (https://CRAN.R-project.org/package=coxme), in which the Laplace approximation method is implemented based on the second order Taylor series. However, fitting coxlmm with coxme directly is time-consuming due to the high-dimensional problem. Instead, we fit coxlmm in an efficient alternative way. Specifically, note that the genetic component <img src="https://latex.codecogs.com/gif.latex?G_i^Tb" title="G_i^Tb" /> can be re-expressed as

<div align=center><img src="https://latex.codecogs.com/gif.latex?{u_i}\;&space;=&space;\;G_i^Tb\;\~\;N(0,\;G_i^T{G_i}{\rm{\sigma&space;}}_b^2)" title="{u_i}\; = \;G_i^Tb\;\~\;N(0,\;G_i^T{G_i}{\rm{\sigma }}_b^2)" />
</div>

Based on the above relationship, we construct an equivalent Cox mixed model

<div align=center><img src="https://latex.codecogs.com/gif.latex?h({t_i}|{X_i},\;{G_i})&space;&&space;\;&space;=&space;\;{h_0}({t_i}){e^{X_i^Ta\;&space;&plus;&space;\;Z_i^T\delta&space;}},\;{\delta&space;_j}&space;&&space;\;\~\;N(0,\;{\rm{\sigma&space;}}_b^2)," title="h({t_i}|{X_i},\;{G_i}) & \; = \;{h_0}({t_i}){e^{X_i^Ta\; + \;Z_i^T\delta }},\;{\delta _j} & \;\~\;N(0,\;{\rm{\sigma }}_b^2)," />
</div>

where **Z**<sub>*i*</sub> is the ith row vector of **K**<sup>1/2</sup> with **K** = **GG**<sup>*T*</sup>, which is often referred to as the genetic relationship matrix (GRM) in genetic prediction; ***δ*** is an n-dimensional vector for the effect sizes of **Z**<sub>*i*</sub>.


## References

+ Ripatti S, Palmgren J. Estimation of multivariate frailty models using penalized partial likelihood[J]. Biometrics, 2000, 56(4): 1016-1022. [DOI: 10.1111/j.0006-341X.2000.01016.x](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.0006-341X.2000.01016.x)

+ Therneau T M, Grambsch P M, Pankratz V S. Penalized survival models and frailty[J]. Journal of computational and graphical statistics, 2003, 12(1): 156-175. [DOI: 10.1198/1061860031365](https://amstat.tandfonline.com/doi/abs/10.1198/1061860031365#.XgARlkczYdU)

+ R Core Team. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria, 2013 (https://www.R-project.org/). 

## Cite






## Contact
I am very grateful to any questions, comments, or bugs reports; and please contact [Xinghao Yu](https://github.com/biostatyu) via xinghaoyu@stu.xzhmu.edu.cn or [Ping Zeng](https://github.com/biostatpzeng) via zpstat@xzhmu.edu.cn.
