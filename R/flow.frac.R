#' Flow Fraction
#' @description Computes \code{\link[NetworkToolbox]{eigenvector}} centrality over nodes in a subset of nodes
#' in the network. This measure generalizes across any subset of nodes and
#' is not specific to communities
#' 
#' @param A An adjacency matrix
#' 
#' @param nodes A subset of nodes in the network
#' 
#' @return Returns a flow fraction value
#' 
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' nodes <- seq(1,48,2)
#' 
#' result <- flow.frac(A, nodes)
#' 
#' @references 
#' Giscard, P. L., & Wilson, R. C. (2018).
#' A centrality measure for cycles and subgraphs II.
#' \emph{Applied Network Science}, \emph{3}, 9.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Flow fraction
flow.frac <- function (A, nodes)
{
    #make network diagonal 0
    diag(A) <- 0
    #grab first eigenvalue
    eig <- eigen(A)$values[1]
    #grab number of nodes
    n <- ncol(A)
    
    #remove nodes from network
    A[nodes,] <- 0
    A[,nodes] <- 0
    
    #compute flow.frac
    eye <- matrix(0,nrow=n,ncol=n)
    diag(eye) <- 1
    res <- det(eye - 1/eig*A)
    
    return(res)
}
#----
