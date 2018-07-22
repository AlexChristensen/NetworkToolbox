#' Community Centrality
#' @description Computes the community centrality measure of each
#' community in a network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param comm A vector or matrix corresponding to the
#' community each node belongs to
#' 
#' @return A vector of community centrality values for each specified community in the network
#' (larger values suggest more central positioning)
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' comm <- igraph::walktrap.community(convert2igraph(abs(A)))$membership
#' 
#' result <- comm.cent(A, comm)
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Community Centrality----
comm.cent <- function (A, comm)
{
    if(is.null(comm))
    {stop("comm must be input")}
    
    comm <- as.vector(comm)
    
    if(ncol(A)!=length(comm))
    {stop("length of comm does not match nodes in matrix")}
    
    uniq <- unique(comm)
    len <- length(uniq)
    
    allP <- pathlengths(A)$ASPLi
    mean.allP <- mean(allP)
    remove <- matrix(0,nrow=len,ncol=1)
    
    for(j in 1:len)
    {
        rem <- which(comm==uniq[j])
        
        remove[j,] <- (mean(allP[-rem]))-mean.allP
    }
    
    norm <- remove
    
    for(k in 1:len)
    {norm[k,] <- (remove[k,] - min(remove))/(max(remove)-min(remove))}
    
    norm <- as.vector(round(norm,3))
    
    names(norm) <- uniq
    
    return(norm)
}
#----