#' Gateway Coefficient
#' @description Computes the gateway coefficient for each node. The gateway
#' coefficient measures a node's connections between its community and other communities.
#' Nodes that are solely responsible for inter-community connectivity will
#' have higher gateway coefficient values. Positive and negative signed weights
#' for gateway coefficients are computed separately.
#' 
#' @param A Network adjacency matrix
#' 
#' @param comm A vector of corresponding to each item's community.
#' Defaults to \code{"walktrap"} for the \code{\link[igraph]{cluster_walktrap}} community detection algorithm.
#' Set to \code{"louvain"} for the \code{\link[NetworkToolbox]{louvain}} community detection algorithm.
#' Can also be set to user-specified communities (see examples)
#' 
#' @param cent Centrality to community gateway coefficient.
#' Defaults to \code{"strength"}.
#' Set to \code{"betweenness"} to use the \code{\link[NetworkToolbox]{betweenness}} centrality
#' 
#' @return Returns a list containing:
#' 
#' \item{overall}{Gateway coefficient without signs considered}
#' 
#' \item{positive}{Gateway coefficient with only positive sign}
#' 
#' \item{negative}{Gateway coefficient with only negative sign}
#' 
#' @examples
#' #theoretical communities
#' comm <- rep(1:8, each = 6)
#' 
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' gw <- gateway(A, comm = comm)
#' 
#' #walktrap communities
#' wgw <- gateway(A, comm = "walktrap")
#' 
#' @references
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}, 1059-1069.
#' 
#' Vargas, E. R., & Wahl, L. M. (2014).
#' The gateway coefficient: A novel metric for identifying critical connections in modular networks.
#' \emph{The European Physical Journal B}, \emph{87}, 1-10.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Gateway Coefficient----
gateway <- function (A, comm = c("walktrap","louvain"),
                     cent = c("strength","betweenness"))
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
    
    if(missing(cent))
    {cent<-"strength"
    }else{cent<-match.arg(cent)}
    
    gate <- function (W, t)
    {
        S <- colSums(W)
        Gc <- (W!=0)%*%diag(facts)
        Sc2 <- vector(mode="numeric",length=n)
        ksm <- vector(mode="numeric",length=n)
        centm <- vector(mode="numeric",length=n)
        
        if(t=="strength")
        {cents <- as.vector(S)
        }else if(t=="betweenness")
        {cents <- as.vector(as.matrix(betweenness(W)))}
        
        for(i in 1:max(facts))
        {
            ks <- rowSums(W*(Gc==i))
            Sc2 <- Sc2 + (ks^2)
            for(j in 1:max(facts))
            {
                ksm[facts==j] <- ksm[facts==j] + (ks[facts==j]/(sum(ks[facts==j])))
            }
            centm[facts==i] <- sum(cents[facts==i])
        }
        
        centm <- centm/max(centm)
        
        gs <- (1-(ksm*centm))^2
        
        GW <- (vector(mode="numeric",n)+1) - (Sc2/(S^2)*gs)
        GW[is.na(GW)]<-0
        GW[!GW]<-0
        
        return(GW)
    }
    
    GWpos <- gate(W=ifelse(A>0,A,0),t=cent)
    GWneg <- gate(W=ifelse(A<0,A,0),t=cent)
    GW <- gate(W=A,t=cent)
    
    return(list(overall=GW,positive=GWpos,negative=GWneg))
}
#----