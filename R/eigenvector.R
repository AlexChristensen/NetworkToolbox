#' Eigenvector Centrality
#' @description Computes eigenvector centrality of each node in a network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param weighted Is the network weighted?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} for unweighted measure of eigenvector centrality
#' 
#' @return A vector of eigenvector centrality values for each node in the network
#' 
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' #Weighted
#' EC <- eigenvector(A)
#' 
#' #Unweighted 
#' EC <- eigenvector(A, weighted = FALSE)
#' 
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{NeuroImage}, \emph{52}, 1059-1069.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Eigenvector----
eigenvector <- function (A, weighted = TRUE)
{
    if(nrow(A)!=ncol(A))
    {stop("Input not an adjacency matrix")}
    
    # A <- abs(A)
    A <- as.matrix(A)
    
    if(!weighted)
    {A <- binarize(A)}
    ec <- abs(eigen(A)$vectors[,1])
    ec <- as.vector(round(ec,3))
    names(ec) <- colnames(A)

    return(ec)
}
#----