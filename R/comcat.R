#' Communicating Nodes
#' @description Computes the between-community strength for each node in the network
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
#' @return A vector containing the between-community strength value for each node
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' communicating <- comcat(A, comm = "walktrap")
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
#Communicating----
comcat <- function (A, comm = c("walktrap","louvain"),
                    cent = c("strength","degree"), ...)
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
    {factors<-"walktrap"
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
        Ah<-A[which(facts!=i),which(facts==i)]
        
        if(cent=="degree")
        {com<-colSums(ifelse(Ah!=0,1,0))
        }else if(cent=="strength")
        {com<-colSums(Ah)}
        
        fact[[i]]<-com
    }
    
    commn<-unlist(fact)
    
    bind<-cbind(ord,commn)
    
    commord<-bind[order(bind[,1]),]
    
    commmat<-matrix(commord[,3],nrow=nrow(commord),ncol=1)
    
    commmat <- as.vector(commmat)
    
    names(commmat) <- colnames(A)
    
    return(commmat)
}
#----