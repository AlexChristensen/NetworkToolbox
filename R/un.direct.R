#' Convert Directed Network to Undirected Network
#' @description Converts a directed network to an undirected network
#' 
#' @param A Matrix or data frame.
#' Adjacency matrix (network matrix)
#' 
#' @param diagonal Numeric.
#' Number to be placed on the diagonal.
#' Defaults to \code{0}
#' 
#' @return Returns a symmetric adjacency matrix
#' 
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' # create a directed network
#' dir <- A * sample(c(0,1), size = length(A), replace = TRUE)
#' 
#' # undirect the directed network
#' undir <- un.direct(dir)
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Convert matrices to igraph format----
un.direct <- function (A, diagonal = 0)
{
  # turn network into matrix
  mat <- as.matrix(A)
  
  # check if already symmetric
  if(isSymmetric(mat))
  {
    message("Adjacency matrix is already symmetric")
    return(A)
  }
  
  # triangles
  lower <- mat[lower.tri(mat)]
  upper <- mat[upper.tri(mat)]
  
  # initialize symmetricize vector
  sym <- numeric(length(lower))
  
  # symmetricize
  sym[which(lower == upper)] <- lower[which(lower == upper)]
  sym[which(lower > upper)] <- lower[which(lower > upper)]
  sym[which(upper > lower)] <- upper[which(upper > lower)]
  
  # initialize new matrix
  mat2 <- matrix(0, nrow = nrow(A), ncol = ncol(A))
  
  # put in values
  mat2[lower.tri(mat2)] <- sym
  mat2 <- t(mat2) + mat2
  
  # check if symmetric
  if(isSymmetric(mat2))
  {
    row.names(mat2) <- colnames(A)
    colnames(mat2) <- colnames(A)
    diag(mat2) <- diagonal
    return(mat2)
  }else(stop("Adjacency matrix could not be made symmetric"))
}
#----