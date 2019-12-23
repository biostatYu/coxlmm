# Example of LRT

## Input data
**[coxlmm](https://github.com/biostatYu/coxlmm/blob/master/coxlmm.R)** requires ***y*** (survival information (https://github.com/biostatpzeng/LRT/blob/master/phenotype.fam)), ***X*** (clinical variable (https://github.com/biostatYu/coxlmm/blob/master/example/x.txt)) and ***G*** (gene expression levels (https://github.com/biostatYu/coxlmm/blob/master/example/g.txt)). It first computes the ***K*** matrix: ***K*** = GG^T (a n by n matrix measureed the similarity of individuals). For the input data, no missing data is allowed. So, missing data should be removed before data analysis.

## model fit of prognosis with *coxlmm* (linear mixed effect Cox model with clinical covariates and gene expressions levels)
x = read.table("x.txt", header = T) #survival time and status
y = read.table("y.txt", header = T) #clinical covariates
g = read.table("g.txt", header = T) #gene expression levels
source("coxlmm.R")
fit = coxlmm(y,x,g)

$fit #model fit with coxme
Cox mixed-effects model fit by maximum likelihood

  events, n = 3, 5
  Iterations= 2 52 
                    NULL    Integrated        Fitted
Log-likelihood -4.094345 -1.201705e-10 -4.170175e-11

                  Chisq df        p  AIC  BIC
Integrated loglik  8.19  4 0.084906 0.19 3.79
 Penalized loglik  8.19  3 0.042269 2.19 4.89

Model:  y_surv ~ x + (K12 | 1) 
Fixed coefficients
              coef    exp(coef) se(coef)   z p
xage    0.09293339 1.097389e+00 123806.3   0 1
xsex   36.56749820 7.604360e+15      0.0 Inf 0
xstage 23.53996603 1.672158e+10 115138.0   0 1

Random effects
 Group Variable    Std Dev    Variance  
 1     (Shrinkage) 0.12930908 0.01672084

$beta_x # beta of x
       xage        xsex      xstage 
 0.09293339 36.56749820 23.53996603 

$beta_g # beta of g
                      [,1]
ARHGEF10L    -1.629072e-13
HIF3A        -2.103226e-13
RNF17         6.609706e-13
RNF11        -6.164177e-13
RNF13        -4.116933e-13
GTF2IP1       2.273647e-13
REM1          3.273254e-14
...

$lp #linear prediction
           [,1]
[1,]  34.778347
[2,]   9.779199
[3,] -15.374094
[4,]  59.643957
[5,] -40.368868

## prediction of prognosis with *coxlmm*

predict = coxlmm_pred(y,x,g,fit)
$lpnew # linear prediction with test data
           [,1]
[1,]  34.778347
[2,]   9.779199
[3,] -15.374094
[4,]  59.643957
[5,] -40.368868

$cindex #concordance index
[1] 1

## estimation of the partition of clinical variance and genetic variance

PVE.est = PVE(y,x,g,fit)
$G1
         [,1]
[1,] 1564.732

$G2
             [,1]
[1,] 4.898075e-22

$PCE # proportion of clinical variable
          [,1]
[1,] 0.9989498

$PGE # proportion of gene expression levels
             [,1]
[1,] 3.127008e-25



