#' Community Eigenvector Centrality
#' @description Computes the \link[NetworkToolbox]{flow.frac} for each
#' community in the network. The values are equivalent to the community's
#' eigenvector centrality
#' 
#' @param A An adjacency matrix
#' 
#' @param comm A vector or matrix corresponding to the
#' community each node belongs to
#' 
#' @return A vector of community eigenvector centrality values for
#' each specified community in the network
#' (larger values suggest more central positioning)
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' comm <- igraph::walktrap.community(convert2igraph(abs(A)))$membership
#' 
#' result <- comm.eigen(A, comm)
#' 
#' @references 
#' Giscard, P. L., & Wilson, R. C. (2018).
#' A centrality measure for cycles and subgraphs II.
#' \emph{Applied Network Science}, \emph{3}, 9.
#' doi: \href{https://doi.org/10.1007/s41109-018-0064-5}{10.1007/s41109-018-0064-5}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Community EC
comm.eigen <- function (A, comm)
{
    if(is.null(comm))
    {stop("comm must be input")}
    
    comm <- as.vector(comm)
    
    if(ncol(A)!=length(comm))
    {stop("length of comm does not match nodes in matrix")}
    
    uniq <- unique(comm)
    len <- length(uniq)
    
    commEC <- vector("numeric",length=len)
    
    for(i in 1:len)
    {commEC[i] <- flow.frac(A, which(comm==uniq[i]))}
    
    names(commEC) <- uniq
    
    return(commEC)
}
#----