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
#' \item{negative}{Diversity coefficient with only negative sign}
#' 
#' @details 
#' Values closer to 1 suggest greater between-community connectivity and 
#' values closer to 0 suggest greater within-community connectivity
#' 
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' #theoretical communities
#' comm <- rep(1:8, each = 6)
#' 
#' gdiv <- diversity(A, comm = comm)
#' 
#' #walktrap communities
#' wdiv <- diversity(A, comm = "walktrap")
#' 
#' @references
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{NeuroImage}, \emph{52}, 1059-1069.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Diversity Coefficient----
# Updated 01.05.2022
diversity <- function (A, comm = c("walktrap","louvain"))
{
    #convert to matrix
    A <- as.matrix(A)
    
    #nodes
    n <- ncol(A)
    
    #set diagonal to zero
    diag(A) <- 0
    
    #set communities
    if(missing(comm))
    {comm<-"walktrap"
    }else{comm<-comm}
    
    
    #check if comm is character
    if(is.character(comm))
    {
        if(length(comm) == 1)
        {
            facts <- switch(comm,
                            walktrap = suppressWarnings(igraph::walktrap.community(convert2igraph(A))$membership),
                            louvain = suppressWarnings(louvain(A)$community)
            )
        }else{
            
            uni <- unique(comm)
            
            facts <- comm
            
            for(i in 1:length(uni))
            {facts[which(facts==uni[i])] <- i}
            
        }
        
    }else{facts <- comm}
    
    #ensure communities are numeric
    facts <- as.numeric(facts)
    
    #number of communities
    m <- max(facts)
    
    ent <- function (A, facts)
    {
        S <- colSums(A)
        Snm <- matrix(0,nrow=n,ncol=m)
        
        for(i in 1:m)
        {
            if(is.vector(A[,which(facts==i)]))
            {Snm[,i] <- sum(A[,which(facts==i)])
            }else{Snm[,i] <- rowSums(A[,which(facts==i)])}
        }
        
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