#' Root Mean Square Error
#' @description Computes the root mean square error of a sparse model to a full model
#' 
#' @param base Base (or full) model to be evaulated against
#' 
#' @param test Reduced (or testing) model (e.g., a sparse correlation or covariance matrix)
#' 
#' @return RMSE value (lower values suggest more similarity between the full and sparse model)
#' 
#' @examples
#' A1 <- solve(cov(neoOpen))
#' 
#' A2 <- LoGo(neoOpen)
#' 
#' root <- rmse(A1, A2)
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Root Mean Square Error----
rmse <- function (base, test)
{
    base <- as.vector(base)
    test <- as.vector(test)
    
    error <- base - test
    
    root <- sqrt(mean(error^2))
    
    return(root)
}
#----