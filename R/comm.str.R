#' Community Strength/Degree Centrality
#' @description Computes the community strength/degree centrality measure of each
#' community in a network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param comm A vector corresponding to the
#' community each node belongs to
#' 
#' @param weighted Is the network weighted?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} for weighted measures
#' 
#' @return A vector of community strength/degree centrality values for each specified
#' community in the network
#' (larger values suggest more central positioning)
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' comm <- igraph::walktrap.community(convert2igraph(abs(A)))$membership
#' 
#' #Strength
#' result <- comm.str(A, comm)
#' 
#' #Degree
#' result <- comm.str(A, comm, weighted = FALSE)
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Community Strength/Degree Centrality
comm.str <- function (A, comm, weighted = TRUE)
{
    if(is.null(comm))
    {stop("comm must be input")}
    
    comm <- as.vector(comm)
    
    if(ncol(A)!=length(comm))
    {stop("length of comm does not match nodes in matrix")}
    
    uniq <- unique(comm)
    len <- length(uniq)
    
    remove <- matrix(0,nrow=len,ncol=1)
    
    for(j in 1:len)
    {
        rem <- which(comm==uniq[j])
        
        if(length(rem)!=1)
        {
            if(weighted)
            {remove[j,] <- sum(colSums(A[,rem]))
            }else{remove[j,] <- sum(colSums(binarize(A)[,rem]))}
        }else{remove[j,] <- 0}
    }
    
    norm <- as.vector(round(remove,3))
    
    names(norm) <- uniq
    
    return(norm)
}
#----