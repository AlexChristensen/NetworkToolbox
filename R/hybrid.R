#' Hybrid Centrality
#' @description Computes hybrid centrality of each node in a network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param BC How should the betweenness centrality be computed?
#' Defaults to \code{"random"}.
#' Set to \code{"standard"} for standard \code{\link[NetworkToolbox]{betweenness}}.
#' 
#' @param beta Beta parameter to be passed to the \code{\link[NetworkToolbox]{rspbc}} function
#' Defaults to .01
#' 
#' @return A vector of hybrid centrality values for each node in the network
#' (higher values are more central, lower values are more peripheral)
#' 
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' HC <- hybrid(A)
#' @references 
#' Christensen, A. P., Kenett, Y. N., Aste, T., Silvia, P. J., & Kwapil, T. R. (2018).
#' Network structure of the Wisconsin Schizotypy Scales-Short Forms:
#' Examining psychometric network filtering approaches.
#' \emph{Behavior Research Methods}, 1-20.
#' doi: \href{https://doi.org/10.3758/s13428-018-1032-9}{10.3758/s13428-018-1032-9}
#' 
#' Pozzi, F., Di Matteo, T., & Aste, T. (2013).
#' Spread of risk across financial markets: Better to invest in the peripheries. 
#' \emph{Scientific Reports}, \emph{3}, 1655.
#' doi: \href{https://doi.org/10.1038/srep01665}{10.1038/srep01665}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Hybrid Centality----
hybrid <- function (A, BC = c("standard","random"), beta)
{
    A <- abs(A)
    A <- as.matrix(A)
    
    if(missing(BC))
    {BC<-"random"
    }else{BC<-match.arg(BC)}
    
    if(missing(beta))
    {beta<-.01
    }else{beta<-beta}
    
    if(nrow(A)!=ncol(A))
    {stop("Input not an adjacency matrix")}
    if(BC=="standard")
    {
        BCu<-betweenness(A,weighted=FALSE)
        BCw<-betweenness(A)
    }else if(BC=="random")
    {
        BCu<-rspbc(binarize(A),beta=beta)
        BCw<-rspbc(A,beta=beta)
    }

    CCu<-closeness(A,weighted=FALSE)
    CCw<-closeness(A)
    if(isSymmetric(A, check.attributes = FALSE))
    {Deg<-degree(A)
    }else{Deg<-degree(A)$outDeg}
    if(isSymmetric(A, check.attributes = FALSE))
    {Str<-strength(A)
    }else{Str<-strength(A)$outStr}
    ECu<-eigenvector(A,weighted=FALSE)
    ECw<-eigenvector(A)
    #levu<-leverage(A,weighted=FALSE)
    #levw<-leverage(A)
    #Eu<-PathLengths(A,weighted=FALSE)$ecc
    #Ew<-PathLengths(A)$ecc
    
    hyb<-((rank(BCu,ties.method="max")+
                  rank(BCw,ties.method="max")+
                  rank(CCu,ties.method="max")+
                  rank(CCw,ties.method="max")+
                  rank(Deg,ties.method="max")+
                  rank(Str,ties.method="max")+
                  rank(ECu,ties.method="max")+
                  rank(ECw,ties.method="max")-
                  #rank(levu,ties.method="max")+
                  #rank(levw,ties.method="max")-
                  #rev(rank(Eu,ties.method="max"))+
                  #rev(rank(Ew,ties.method="max"))-
                  8)/(8*((ncol(A))-1)))
    
    hyb<-round(as.vector(hyb),3)
    
    names(hyb) <- colnames(A)
    
    return(hyb)
}
#----