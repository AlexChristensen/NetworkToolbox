#' Eigenvector Centrality
#' @description Computes eigenvector centrality of each node in a network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param weighted Is the network weighted?
#' Defaults to TRUE.
#' Set to FALSE for unweighted measure of eigenvector centrality
#' 
#' @return A vector of eigenvector centrality values for each node in the network
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
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
#' \emph{Neuroimage}, \emph{52}, 1059-1069.
#' doi: \href{https://doi.org/10.1016/j.neuroimage.2009.10.003}{10.1016/j.neuroimage.2009.10.003}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Eigenvector----
eigenvector <- function (A, weighted = TRUE)
{
    if(nrow(A)!=ncol(A))
    {stop("Input not an adjacency matrix")}
    
    A <- abs(A)
    A <- as.matrix(A)
    
    if(!weighted)
    {A <- binarize(A)}
    ec <- abs(eigen(A)$vectors[,1])
    ec <- as.vector(round(ec,3))
    names(ec) <- colnames(A)

    return(ec)
}
#----