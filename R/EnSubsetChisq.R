
#' Ensemble Subset Chi-squared test
#'
#' @param X a numeric vector of z-scores.
#' @param Sigma the correlation matrix of the z-scores.
#' @param B the number of base tests.
#' @param m the length of the subset of X. If not provided, m is set to be floor(sqrt(length(X))).
#' @param is.Sigma.identity logical. If \code{TRUE}, \code{Sigma} will be set to be an identity matrix and the implementation would be more efficient. 
#' @param is.pvals.path logical. If \code{TRUE}, the p-value path of the ensemble test as \emph{B} increases will be provided in the output and can be used to draw the ensemble p-value path plot (See examples).
#' 
#' @return A list with components:
#' \describe{
#' \item{\code{pval.ensemble.test}}{The p-value of the ensemble test}
#' \item{\code{pval.base.test}}{The p-values of the individual base tests}
#' \item{\code{pval.path.ensemble}}{The p-value path of the ensemble test when \emph{is.pvals.path == TRUE}.}
#' }
#'
#'
#' @author Yaowu Liu
#' 
#' @references Liu, Y., Liu, Z., and Lin, X. (2024+) Ensemble methods for testing a global null.
#' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}. Accepted.
#' 
#' @examples #### Example 1, Sigma is an identity matrix
#' @examples set.seed(1234); p = 20
#' @examples X = rnorm(p)
#' @examples X[1:3] = X[1:3]+3.5   ## z-scores
#' @examples res = EnSubsetChisq(X,Sigma = diag(p),B=1000,is.Sigma.identity = T, is.pvals.path = T)
#' @examples res$pval.ensemble.test   ## the ensemble p-value
#' 
#' @examples pvals.path = res$pval.path.ensemble   ## p-value path plot
#' @examples plot(c(1:length(pvals.path)), -log10(pvals.path),xlab = "Number of base tests",ylab = "-log10(p-value)",type = "l")
#' 
#' @examples #### Example 2, Sigma is not an identity matrix
#' @examples library(MASS)
#' @examples set.seed(1234); p = 20
#' @examples Sigma = matrix(0.1,nrow = p, ncol = p); diag(Sigma) = 1;
#' @examples X = mvrnorm(1,rep(0,p),Sigma)
#' @examples X[1:3] = X[1:3]+3.5   ## z-scores
#' @examples res = EnSubsetChisq(X,Sigma =Sigma,B=1000, is.pvals.path = T)
#' @examples res$pval.ensemble.test  ## the ensemble p-value
#'
#' @examples pvals.path = res$pval.path.ensemble  ## p-value path plot
#' @examples plot(c(1:length(pvals.path)), -log10(pvals.path),xlab = "Number of base tests",ylab = "-log10(p-value)",type = "l")
#' 
#' @export
EnSubsetChisq <- function (X,Sigma,B,m = NULL,is.Sigma.identity = FALSE,is.pvals.path = FALSE){
    #### check
    if (!is.vector(X)){
        stop("X should be number vector!")
    }
    p = length(X)
    
    if (is.Sigma.identity){
        Sigma = diag(p)
    }else{
        if (nrow(Sigma)!=p || ncol(Sigma)!=p){
            stop("Sigma should be a square matrix whose dimension equals to the length of X!")
        }
        
        if (sum(diag(Sigma) != 1)>0){
            stop("Sigma should be a correlation matrix whose diagonals are all 1.")
        }
    }
    
    if (is.null(m)){
        m = floor(sqrt(p))
    }else if (m >= p){
        stop("m should be less than the length of X!")
    }
    
    #### test
    CHI.stat<-Get_chistat(X=as.matrix(X),Sigma=Sigma,B=B,m_select=m,is_Sigma_identity=is.Sigma.identity)
    
    out = list()
    out[["pval.ensemble.test"]] = NA
    out[["pval.base.test"]] = pchisq(as.vector(CHI.stat),df=m,lower.tail = F)
    out[["pval.ensemble.test"]] = ACAT(out[["pval.base.test"]])
    
    if (is.pvals.path){
        out[["pval.path.ensemble"]] = ACAT.path(out[["pval.base.test"]])
    }else{
        out[["pval.path.ensemble"]] = NULL
    }
    
    return(out)
}

ACAT.path<-function(Pvals){
    if (is.vector(Pvals)){
        Pvals = as.matrix(Pvals,ncol = 1)
    }
    
    if (!is.matrix(Pvals)){
        stop("The input Pvals must be a matrix or vector!")
    }
    
    v.cauchy = tan((0.5 - Pvals) * pi)
    for (i in 2:nrow(v.cauchy)){
        v.cauchy[i,] = v.cauchy[i,]+v.cauchy[(i-1),]
        v.cauchy[(i-1),] = v.cauchy[(i-1),]/(i-1)
    }
    v.cauchy[i,] = v.cauchy[i,]/i
    
    acat.pvals = pcauchy(v.cauchy, lower.tail = F)
    
    return(acat.pvals)
}


