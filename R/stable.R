#' Stabilizing Nodes
#' @description Computes the within-community centrality for each node in the network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param comm Can be a vector of community assignments or community detection algorithms
#' ("walktrap" or "louvain") can be used to determine the number of factors.
#' Defaults to "walktrap".
#' Set to "louvain" for louvain community detection
#' 
#' @param cent Centrality measure to be used.
#' Defaults to "strength".
#' 
#' @param ... Additional arguments for community detection algorithms
#' 
#' @return A matrix containing the within-community centrality value for each node
#' 
#' @examples
#' A<-TMFG(neoOpen)$A
#' 
#' stabilizing <- stable(A, comm = "walktrap")
#' 
#' @references 
#' Blanken, T. F., Deserno, M. K., Dalege, J., Borsboom, D., Blanken, P., Kerkhof, G. A., & Cramer, A. O. (2018).
#' The role of stabilizing and communicating symptoms given overlapping communities in psychopathology networks.
#' \emph{Scientific Reports}, \emph{8}, 5854.
#' doi: \href{https://doi.org/10.1038/s41598-018-24224-2}{10.1038/s41598-018-24224-2}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Stabilizing----
stable <- function (A, comm = c("walktrap","louvain"),
                    cent = c("betweenness","rspbc","strength","degree","hybrid"), ...)
{
    #nodes
    n <- ncol(A)
    
    #set diagonal to zero
    diag(A) <- 0
    
    if(missing(cent))
    {cent<-"strength"
    }else{cent<-match.arg(cent)}
    
    #set communities
    if(missing(comm))
    {comm<-"walktrap"
    }else{comm<-comm}
    
    if(!is.numeric(comm))
    {
        if(comm=="walktrap")
        {facts<-igraph::walktrap.community(convert2igraph(A),...)$membership
        }else if(comm=="louvain")
        {facts<-louvain(A,...)$community}
    }else{facts<-comm}
    
    keepord<-cbind(rep(1:ncol(A)),facts)
    ord<-keepord[order(keepord[,2]),]
    
    fact<-list()
    
    for(i in 1:max(facts))
    {
        Ah<-A[which(facts==i),which(facts==i)]
        
        if(cent=="betweenness")
        {stab<-betweenness(Ah)
        }else if(cent=="rspbc")
        {stab<-rspbc(Ah)
        }else if(cent=="strength")
        {stab<-strength(Ah)
        }else if(cent=="degree")
        {stab<-degree(Ah)}
        
        fact[[i]]<-stab
    }
    
    stabil<-unlist(fact)
    
    bind<-cbind(ord,stabil)
    
    stabord<-bind[order(bind[,1]),]
    
    stabmat<-matrix(stabord[,3],nrow=nrow(stabord),ncol=1)
    
    stabmat <- as.vector(stabmat)
    
    names(stabmat)<-colnames(A)
    
    return(stabmat)
}
#----