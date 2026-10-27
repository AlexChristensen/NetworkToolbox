#' Transitivity
#' @description Computes transitivity of a network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param weighted Is the network weighted?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} for a weighted measure of transitivity
#' 
#' @return Returns a value of transitivity
#' 
#' @details
#' Translated from the MATLAB functions \code{transitivity_bu.m} (Alexandros
#' Goulas) and \code{transitivity_wu.m} (Mika Rubinov) in the Brain Connectivity
#' Toolbox (Rubinov & Sporns, 2010).
#' 
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' trans <- transitivity(A, weighted=TRUE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{NeuroImage}, \emph{52}, 1059-1069.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Transitivity----
transitivity <- function (A, weighted = FALSE)
{
    A <- as.matrix(A)

    if(nrow(A)!=ncol(A))
    {stop("Input not an adjacency matrix")}

    # self-loops do not form triangles
    diag(A) <- 0

    if(!weighted)
    {
        A<-ifelse(A!=0,1,0)
        trans<-sum(diag(A%*%A%*%A))/((sum(A%*%A))-sum(diag(A%*%A)))
    }else if(weighted){
        K<-colSums(ifelse(A!=0,1,0))
        W<-A^(1/3)
        cyc<-diag(W%*%W%*%W)
        trans<-sum(cyc)/sum(K*(K-1))
    }
    
    return(trans)
}
#----