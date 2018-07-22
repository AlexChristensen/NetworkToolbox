#' Convert Correlation Matrix to Covariance Matrix
#' @description Converts a correlation matrix to a covariance matrix
#' 
#' @param cormat A correlation matrix
#' 
#' @param data The dataset the correlation matrix is from
#' 
#' @return Returns a covariance matrix
#' 
#' @examples
#' cormat <- cor(neoOpen)
#' 
#' covmat <- cor2cov(cormat,neoOpen)
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Corrlation to covariance----
cor2cov <- function (cormat, data)
{
    sds<-apply(data,2,sd)
    
    b<-sds%*%t(sds)
    
    S<-cormat*b
    
    return(S)
}
#----