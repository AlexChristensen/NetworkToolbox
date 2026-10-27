#' Determines if Network is Graphical
#' @description Tests for whether the network is graphical.
#' Input must be a partial correlation network.
#' Function assumes that partial correlations were computed from a multivariate normal distribution
#' 
#' @param A A partial correlation network (adjacency matrix)
#' 
#' @return Returns a TRUE/FALSE for whether network is graphical
#' 
#' @examples
#' \dontrun{
#' A <- LoGo(neoOpen, normal = TRUE, partial = TRUE)
#' 
#' is.graphical(A)
#' }
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Is network graphical?----
is.graphical <- function (A)
{
    A <- as.matrix(A)
    
    #partial correlations must be symmetric
    if(!isSymmetric(A, check.attributes = FALSE))
    {return(FALSE)}
    
    #make diagonal 0
    diag(A)<-0
    
    #convert partial correlations to standardized inverse covariance
    I<-diag(1, dim(A)[1])
    inv<-I-A
    
    #reduce small values to 0
    check<-zapsmall(inv)
    
    #inverse covariance matrix must be positive (semi-)definite
    values<-try(eigen(check, symmetric = TRUE, only.values = TRUE)$values,silent=TRUE)
    if(is.character(values))
    {ret<-FALSE
    }else if(any(values<0))
    {ret<-FALSE
    }else{ret<-TRUE}
    
    return(ret)
}
#----
