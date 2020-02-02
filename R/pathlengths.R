#' Characteristic Path Lengths
#' @description Computes global average shortest path length,
#' local average shortest path length, eccentricity,
#' and diameter of a network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param weighted Is the network weighted?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} for weighted measures
#' 
#' @return Returns a list containing:
#' 
#' \item{ASPL}{Global average shortest path length}
#' 
#' \item{ASPLi}{Local average shortest path length}
#' 
#' \item{ecc}{Eccentricity (i.e., maximal shortest path length between a node and any other node)}
#' 
#' \item{D}{Diameter of the network (i.e., the maximum of eccentricity)}
#' 
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' #Unweighted
#' PL <- pathlengths(A)
#' 
#' #Weighted
#' PL <- pathlengths(A, weighted = TRUE)
#' 
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}, 1059-1069.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Path Lengths----
pathlengths <- function (A, weighted = FALSE)
{
    if(nrow(A)!=ncol(A))
    {stop("Input not an adjacency matrix")}
    if(!weighted)
    {D<-distance(A,weighted=FALSE)}else if(weighted){D<-distance(A,weighted=TRUE)}
    n<-nrow(D)
    for(i in 1:ncol(D))
        for(j in 1:nrow(D))
            if(is.infinite(D[j,i]))
            {D[j,i]<-0}
    if(any(colSums(D)==0))
    {
        newD <- D
        newD[,(which(colSums(D)==0))] <- rep(NA,length(which(colSums(D)==0)))
    }else{newD <- D}
    
    aspli<-colSums(newD*(newD!=0))/(ncol(newD)-1)
    aspl<-mean(aspli,na.rm=TRUE)
    
    Emat<-(D*(D!=0))
    
    ecc<-matrix(nrow=nrow(Emat),ncol=1)
    
    for(i in 1:nrow(Emat))
    {ecc[i,]<-max(Emat[i,])}
    
    d<-max(ecc)
    
    ecc <- as.vector(ecc)
    names(ecc) <- colnames(A)
    
    aspli <- as.vector(aspli)
    names(aspli) <- colnames(A)
    
    return(list(ASPL=aspl,ASPLi=aspli,ecc=ecc,diameter=d))
}
#----