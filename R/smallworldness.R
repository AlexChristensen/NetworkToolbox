#' Small-worldness Measure
#' @description Computes the small-worldness measure of a network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param iter Number of random (or lattice) networks to generate,
#' which are used to calculate the mean random ASPL and CC (or lattice)
#' 
#' @param progBar Defaults to \code{FALSE}.
#' Set to \code{TRUE} to see progress bar
#' 
#' @param method Defaults to \code{"HG"} (Humphries & Gurney, 2008).
#' Set to \code{"rand"} for the CC to be calculated using a random network or
#' set to \code{"TJHBL"} for (Telesford et al., 2011) where CC is calculated from a lattice network
#' 
#' @return Returns a list containing:
#' 
#' \item{swm}{Small-worldness value}
#' 
#' \item{rASPL}{Global average shortest path length from random network}
#' 
#' \item{lrCCt}{When \code{"rand"}, clustering coefficient from a random network.
#' When \code{"HG"}, transitivity from a random network.
#' When \code{"TJHBL"}, clustering coefficient from a lattice network}
#' 
#' @details
#' For \code{"rand"}, values > 1 indicate a small-world network.
#' For \code{"HG"}, values > 3 indicate a small-world network.
#' For \code{"TJHBL"}, values near 0 indicate a small-world network,
#' while < 0 indicates a more regular network and > 0 indicates a more random network
#' 
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#'
#' swmHG <- smallworldness(A, method="HG")
#' 
#' swmRand <- smallworldness(A, method="rand")
#' 
#' swmTJHBL <- smallworldness(A, method="TJHBL")
#' 
#' @references 
#' Humphries, M. D., & Gurney, K. (2008).
#' Network 'small-world-ness': A quantitative method for determining canonical network equivalence.
#' \emph{PLoS one}, \emph{3}, e0002051.
#' 
#' Telesford, Q. K., Joyce, K. E., Hayasaka, S., Burdette, J. H., & Laurienti, P. J. (2011).
#' The ubiquity of small-world networks.
#' \emph{Brain Connectivity}, \emph{1}(5), 367-375.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Small-worldness Measure----
smallworldness <- function (A, iter = 100, progBar = FALSE, method = c("HG","rand","TJHBL"))
{
    if(missing(method))
    {method<-"HG"
    }else{method<-match.arg(method)}
    
    mat<-matrix(0,nrow=nrow(A),ncol=ncol(A)) #Initialize bootstrap matrix
    asamps<-matrix(0,nrow=iter) #Initialize sample matrix
    csamps<-matrix(0,nrow=iter) #Initialize sample matrix
    if(progBar)
    {pb <- txtProgressBar(max=iter, style = 3)}
    for(i in 1:iter) #Generate array of bootstrapped samples
    {
        f<-round(runif(i,min=1,max=1000000),0)
        set.seed(f[round(runif(i,min=1,max=length(f)),0)])
        rand<-randnet(ncol(A),sum(ifelse(A!=0,1,0))/2)
        if(method=="TJHBL")
        {latt<-lattnet(ncol(A),sum(ifelse(A!=0,1,0))/2)}
        asamps[i,]<-pathlengths(rand)$ASPL
        if(method=="rand")
        {csamps[i,]<-clustcoeff(rand)$CC
        }else if(method=="HG"){csamps[i,]<-transitivity(rand)
        }else if(method=="TJHBL"){csamps[i,]<-clustcoeff(latt)$CC}else{stop("Method not available")}
        if(progBar)
        {setTxtProgressBar(pb, i)}
    }
    if(progBar)
    {close(pb)}
    
    nodes<-ncol(A)
    ASPL<-pathlengths(A)$ASPL
    CC<-clustcoeff(A)$CC
    trans<-transitivity(A)
    rASPL<-mean(asamps)
    
    if(method=="rand")
    {rCC<-mean(csamps)
    swm<-(CC/rCC)/(ASPL/rASPL)
    lrCCt<-rCC
    }else if(method=="HG")
    {rtrans<-mean(csamps)
    swm<-(trans/rtrans)/(ASPL/rASPL)
    lrCCt<-rtrans
    }else if(method=="TJHBL")
    {lCC<-mean(csamps)
    swm<-(rASPL/ASPL)-(CC/lCC)
    lrCCt<-lCC
    }
    
    return(list(swm=swm, rASPL=rASPL, lrCCt=lrCCt))
}
#----