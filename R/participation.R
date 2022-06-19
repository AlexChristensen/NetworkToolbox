#' Participation Coefficient
#' @description Computes the participation coefficient for each node. The participation
#' coefficient measures the strength of a node's connections within its community. Positive
#' and negative signed weights for participation coefficients are computed separately.
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
#' \item{overall}{Participation coefficient without signs considered}
#' 
#' \item{positive}{Participation coefficient with only positive sign}
#' 
#' \item{negative}{Participation coefficient with only negative sign}
#' 
#' @details 
#' Values closer to 0 suggest greater within-community connectivity and 
#' values closer to 1 suggest greater between-community connectivity
#' 
#' @examples
#' #theoretical factors
#' comm <- rep(1:8, each = 6)
#' 
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' pc <- participation(A, comm = comm)
#' 
#' # Walktrap factors
#' wpc <- participation(A, comm = "walktrap")
#' 
#' @references
#' Guimera, R., & Amaral, L. A. N. (2005).
#' Functional cartography of complex metabolic networks.
#' \emph{Nature}, \emph{433}, 895-900.
#' 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{NeuroImage}, \emph{52}, 1059-1069.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Participation Coefficient----
# Updated 19.06.2022
participation <- function (A, comm = c("walktrap","louvain"))
{
    #make sure its a matrix
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
    if(is.character(comm)){
      
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
    
    #participation coefficient
    pcoef <- function (A, facts){
      
        k <- colSums(A) #strength
        Gc <- facts #communities
        Kc2 <- vector(mode="numeric",length=n)  
        
        for(i in 1:max(Gc)){
          
            #strength within communities squared
            if(is.vector(A*(Gc==i)))
            {Kc2 <- Kc2 + sum(A*(Gc==i))^2
            }else{Kc2 <- Kc2 + colSums(A*(Gc==i))^2}
        } 
        
        ones <- vector(mode="numeric",length=n) + 1
        
        P <- ones - Kc2 / (k^2)
        
        P[!k] <- 0
        
        return(P)
    }
    
    overall <- pcoef(A, facts) #overall participation coefficient
    
    #signed participation coefficient
    poswei <- ifelse(A>=0,A,0) #positive weights
    negwei <- ifelse(A<=0,A,0) #negative weights
    
    pos <- pcoef(poswei, facts) #positive  participation coefficient
    if(all(pos==1))
    {pos<-1-pos}
    neg <- pcoef(negwei, facts) #negative participation coefficient
    if(all(neg==1))
    {neg<-1-neg}
    
    return(list(overall=overall,positive=pos,negative=neg))
}
#----