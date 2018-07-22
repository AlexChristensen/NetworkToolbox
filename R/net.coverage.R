#' Network Coverage
#' @description Computes the mean distance across a subset of nodes in a network.
#' This measure can be used to identify the effectiveness of a subset of nodes'
#' coverage of the network space
#' 
#' @param A An adjacency matrix
#' 
#' @param nodes Subset of nodes to examine the coverage of the network
#' 
#' @return Returns a list containing:
#' 
#' \item{mean}{The average distance from the subset of nodes to all other nodes in the network}
#' 
#' \item{sd}{The standard deviation of distance from the subset of nodes to all other nodes in the network}
#' 
#' \item{range}{The range of distance from the subset of nodes to all other nodes in the network}
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' nodes <- c(1,3,5,7,11,13,17,19,23,29,31,37,41,43,47)
#' 
#' result <- net.coverage(A, nodes)
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' Mathias Benedek <mathias.benedek@uni-graz.at>
#' 
#' @export
#Network Coverage----
net.coverage <- function (A, nodes)
{
    fnames <- colnames(A)
    
    if(is.numeric(nodes))
    {inames <- colnames(A)[nodes]
    }else{inames <- nodes}
    
    dist <- distance(A)
    
    diag(dist) <- Inf
    
    sdist <- dist[fnames,inames]
    
    res <- list()
    
    res$mean <- mean(apply(sdist,1,min))
    res$sd <- sd(apply(sdist,1,min))
    res$range <- range(apply(sdist,1,min))
    
    return(res)
}
#----