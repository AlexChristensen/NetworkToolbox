#' Network Coverage
#' @description Computes the mean distance across a subset of nodes in a network.
#' This measure can be used to identify the effectiveness of a subset of nodes'
#' coverage of the network space
#' 
#' @param A An adjacency matrix
#' 
#' @param nodes Subset of nodes to examine the coverage of the network
#' 
#' @param weighted Is the network weighted?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} for weighted measures
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
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' nodes <- seq(1,48,2)
#' 
#' result <- network.coverage(A, nodes)
#' 
#' @references 
#' Christensen, A. P., Cotter, K. N., Silvia, P. J., & Benedek, M. (2018)
#' Scale development via network analysis: A comprehensive and concise measure of Openness to Experience
#' \emph{PsyArXiv}, 1-40.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com> and Mathias Benedek <mathias.benedek@uni-graz.at>
#' 
#' @export
#Network Coverage----
network.coverage <- function (A, nodes, weighted = FALSE)
{
    fnames <- colnames(A)
    
    if(is.numeric(nodes))
    {inames <- colnames(A)[nodes]
    }else{inames <- nodes}
    
    dist <- distance(A, weighted = weighted)
    
    diag(dist) <- Inf
    
    sdist <- dist[fnames,inames]
    
    res <- list()
    
    res$mean <- mean(apply(sdist,1,min))
    res$sd <- sd(apply(sdist,1,min))
    res$range <- range(apply(sdist,1,min))
    res$dist <- apply(sdist,1,min)
    
    return(res)
}
#----