#' isSymmetric Wrapper Function
#' @description Checks if matrix is symmetric when \code{isSymmetric} fails
#' 
#' @param matrix A matrix
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' isSym(A)
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Hybrid Centality----
isSym <- function (matrix)
{
    matrix <- as.matrix(matrix)
    n <- nrow(matrix)
    symmat <- matrix(0,nrow=n,ncol=n)
    
    for(i in 1:n)
        for(j in 1:n)
            if(i!=j)
            {
                if(matrix[i,j]==matrix[j,i])
                {symmat[i,j] <- 0
                }else{symmat[i,j] <- 1}
            }
    
    if(sum(symmat)==0){return(TRUE)}else{return(FALSE)}
}
#----