#' Hybrid Centrality
#' @description Computes hybrid centrality of each node in a network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param BC How should the betweenness centrality be computed?
#' Defaults to "standard".
#' Set to "random" for rspbc.
#' Set to "average" for the average of "standard" and "random"
#' 
#' @param beta Beta parameter to be passed to the \emph{rspbc} function
#' 
#' @return A vector of hybrid centrality values for each node in the network
#' (higher values are more central, lower values are more peripheral)
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
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
hybrid <- function (A, BC = c("standard","random","average"), beta)
{
    A <- abs(A)
    
    if(missing(BC))
    {BC<-"standard"
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
    }else if(BC=="average")
    {
        BCu<-rowMeans(cbind(betweenness(A,weighted=FALSE),rspbc(binarize(A))))
        BCw<-rowMeans(cbind(betweenness(A),rspbc(A)))
    }
    CCu<-closeness(A,weighted=FALSE)
    CCw<-closeness(A)
    if(isSymmetric(A))
    {Deg<-degree(A)
    }else{Deg<-degree(A)$outDeg}
    if(isSymmetric(A))
    {Str<-strength(A)
    }else{Str<-strength(A)$outStr}
    ECu<-eigenvector(A,weighted=FALSE)
    ECw<-eigenvector(A)
    #levu<-leverage(A,weighted=FALSE)
    #levw<-leverage(A)
    #Eu<-PathLengths(A,weighted=FALSE)$ecc
    #Ew<-PathLengths(A)$ecc
    
    hybrid<-((rank(BCu,ties.method="max")+
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
    
    hybrid<-round(as.vector(hybrid),3)
    
    names(hybrid) <- colnames(A)
    
    return(hybrid)
}
#----