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
    #make diagonal 1
    diag(A)<-1
    
    #covert partial correlations to covariance
    I<-diag(1, dim(A)[1])
    
    #compute covariance matrix
    error<-try(solve(I-A)%*%t(solve(I-A)),silent=TRUE)
    if(is.character(error))
    {ret<-FALSE
    }else{covmat<-solve(I-A)%*%t(solve(I-A))}
    
    #covert to inverse covariance
    error<-try(solve(covmat),silent=TRUE)
    if(is.character(error))
    {ret<-FALSE
    }else{inv<-solve(covmat)
    
    #reduce small values to 0
    check<-zapsmall(inv)
    
    
    error<-try(any(eigen(check)$values<0),silent=TRUE)
    if(is.character(error))
    {ret<-FALSE
    }else if(any(eigen(check)$values<0))
    {ret<-FALSE
    }else{ret<-TRUE}
    }
    return(ret)
}
#----