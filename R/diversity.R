#' Diversity Coefficient
#' @description Computes the diversity coefficient for each node. The diversity
#' coefficient measures a node's connections to communitites outside of its
#' own community. Nodes that have many connections to other communities will
#' have higher diversity coefficient values. Positive and negative signed weights
#' for diversity coefficients are computed separately.
#' 
#' @param A Network adjacency matrix
#' 
#' @param comm A vector of corresponding to each item's community.
#' Defaults to \code{"walktrap"} for the \code{\link[igraph]{cluster_walktrap}} community detection algorithm.
#' Set to \code{"louvain"} for the \code{\link[NetworkToolbox]{louvain}} community detection algorithm.
#' Can also be set to user-specified communities (see examples)
#' 
#' @return Returns a list containing:
#' 
#' \item{overall}{Diversity coefficient without signs considered}
#' 
#' \item{positive}{Diversity coefficient with only positive sign}
#' 
#' \item{negative}{Diversity coefficient wih only negative sign}
#' 
#' @details 
#' Values closer to 1 suggest greater between-community connectivity and 
#' values closer to 0 suggest greater within-community connectivity
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' #theoretical communities
#' comm <- c(rep(1,8), rep(2,8), rep(3,8), rep(4,8), rep(5,8), rep(6,8))
#' 
#' gdiv <- diversity(A, comm = comm)
#' 
#' #walktrap communities
#' wdiv <- diversity(A, comm = "walktrap")
#' 
#' @references
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}, 1059-1069.
#' doi: \href{https://doi.org/10.1016/j.neuroimage.2009.10.003}{10.1016/j.neuroimage.2009.10.003}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Diversity Coefficient----
diversity <- function (A, comm = c("walktrap","louvain"))
{
    #nodes
    n <- ncol(A)
    
    #set diagonal to zero
    diag(A) <- 0
    
    #set communities
    if(missing(comm))
    {comm<-"walktrap"
    }else{comm<-comm}
    
    if(is.numeric(comm))
    {facts <- comm
    }else{
        if(comm=="walktrap")
        {facts <- igraph::walktrap.community(convert2igraph(A))$membership
        }else if(comm=="louvain")
        {facts <- louvain(A)$community}
    }
    
    #number of communities
    m <- max(facts)
    
    ent <- function (A, facts)
    {
        S <- colSums(A)
        Snm <- matrix(0,nrow=n,ncol=m)
        
        for(i in 1:m)
        {Snm[,i] <- rowSums(A[,facts==i])}
        
        pnm <- Snm/(S*matrix(1,nrow=n,ncol=m))
        pnm[is.na(pnm)]<-0
        pnm[!pnm]<-1
        H <- -rowSums(pnm*log(pnm))/log(m)
        
        return(H)
    }
    
    Hpos <- ent(ifelse(A>0,A,0),facts)
    Hneg <- ent(ifelse(A<0,A,0),facts)
    Hov <- ent(A,facts)
    
    return(list(overall=Hov,positive=Hpos,negative=Hneg))
}
#----