#' Community Closeness Centrality
#' @description Computes the community closeness centrality measure of each
#' community in a network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param comm A vector or matrix corresponding to the
#' community each node belongs to
#' 
#' @param weighted Is the network weighted?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} for weighted measures
#' 
#' @return A vector of community closeness centrality values for each specified
#' community in the network
#' (larger values suggest more central positioning)
#' 
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' comm <- igraph::walktrap.community(convert2igraph(abs(A)))$membership
#' 
#' #Weighted
#' result <- comm.close(A, comm)
#' 
#' #Unweighted
#' result <- comm.close(A, comm, weighted = FALSE)
#'
#' @references 
#' Christensen, A. P. (in press).
#' NetworkToolbox: Methods and measures for brain, cognitive, and psychometric network analysis in R.
#' \emph{The R Journal}, \emph{10}, 422-439.
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Community Closeness Centrality----
comm.close <- function (A, comm, weighted = FALSE)
{
    if(is.null(comm))
    {stop("comm must be input")}
    
    comm <- as.vector(comm)
    
    if(ncol(A)!=length(comm))
    {stop("length of comm does not match nodes in matrix")}
    
    uniq <- unique(comm)
    uniq <- uniq[order(uniq)]
    len <- length(uniq)
    
    allP <- pathlengths(A, weighted = weighted)$ASPLi
    mean.allP <- mean(allP)
    remove <- vector("numeric",length=len)
    
    for(j in 1:len)
    {
        rem <- which(comm==uniq[j])
        
        remove[j] <- 1/mean(allP[rem])
    }
    
    names(remove) <- uniq
    
    return(remove)
}
#----