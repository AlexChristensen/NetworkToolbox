#' Clustering Coefficient
#' @description Computes global clustering coefficient (CC) and local clustering coefficient (CCi)
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param weighted Is the network weighted?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} for weighted measures of CC and CCi
#' @return Returns a list containing:
#' 
#' \item{CC}{Global clustering coefficient. The average clustering coefficient for each node in the network}
#' 
#' \item{CCi}{Local clustering coefficient. The clustering coefficient for each node in the network}
#' 
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' #Unweighted CC
#' CCu <- clustcoeff(A)
#' 
#' #Weighted CC
#' CCw <- clustcoeff(A, weighted=TRUE)
#' 
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{NeuroImage}, \emph{52}, 1059-1069.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Clustering Coefficient----
clustcoeff <- function (A, weighted = FALSE)
{
    diag(A) <- 0
    
    if(nrow(A)!=ncol(A))
    {stop("Input not an adjacency matrix")}
    if(!weighted)
    {n<-ncol(A)
    A<-ifelse(A!=0,1,0)
    C<-matrix(0,nrow=n,ncol=1)
    
    for(i in 1:n)
    {
        v<-which(A[i,]!=0)
        k<-length(v)
        if(k >= 2)
        {
            S<-A[v,v]
            C[i]<-sum(S)/(k^2-k)
        }}
    
    C <- round(as.vector(C),3)
    names(C)<-colnames(A)
    
    CCi<-C
    CC <- mean(C)
    
    }else{
        K<-colSums(A!=0)
        m<-A^(1/3)
        cyc<-diag(m%*%m%*%m)
        K[cyc==0]<-Inf
        C <- as.vector(round(cyc/(K*(K-1)),3))
        names(C)<-colnames(A)
        CCi<-C
        CC<-mean(C)
    }
    
    return(list(CC=CC,CCi=CCi))
}
#----