#' Binarize Network
#' @description Converts weighted adjacency matrix to a binarized adjacency matrix
#' 
#' @param A An adjacency matrix of network data (or an array of matrices)
#' 
#' @return Returns an adjacency matrix of 1's and 0's
#' 
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' neoB <- binarize(A)
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Binarize function----
binarize <- function (A)
{
    A <- as.matrix(A)
    bin <- ifelse(A!=0,1,0)
    row.names(bin) <- row.names(A)
    colnames(bin) <- colnames(A)
    
    return(bin)
}
#----