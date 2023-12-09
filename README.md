# EnsembleTests
This package implemented the ensemble subset chi-squared tests. For the implementation of ensemble Burden, SKAT and MORST tests, please see the function `EnsembleSetTests` in the [MORST](https://github.com/yaowuliu/MORST/blob/master/doc/MORST_manual.pdf) package.

## Installation
```
library(devtools)
devtools::install_github("yaowuliu/EnsembleTests")
```
## Usage
Here is a simple example to run the ensemble subset chi-squared test and to plot the ensemble p-value path.

```{r}
library(MASS)
set.seed(1234)

## simulation setting
p = 20
Sigma = matrix(0.1,nrow = p, ncol = p) ## correlation matrix of z-scores
diag(Sigma) = 1
X = mvrnorm(1,rep(0,p),Sigma)
X[1:3] = X[1:3]+3.5   ## z-scores or z statistics

## run the ensemble subset chi-squared test
res = EnSubsetChisq(X,Sigma =Sigma,B=1000, is.pvals.path = T)
res$pval.ensemble.test  ## the ensemble p-value

## plot the ensemble p-value path
pvals.path = res$pval.path.ensemble  ## p-value path plot
plot(c(1:length(pvals.path)), -log10(pvals.path),xlab = "Number of base tests",ylab = "-log10(p-value)",type = "l")


```
