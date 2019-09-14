#' Kullback-Leibler Divergence
#' @description Estimates the Kullback-Leibler Divergence which measures how one probability distribution
#' diverges from the original distribution (equivalent means are assumed)
#' Matrices \strong{must} be positive definite inverse covariance matrix for accurate measurement.
#' This is a \strong{relative} metric
#' 
#' @param base Full or base model
#' 
#' @param test Reduced or testing model
#' 
#' @return A value greater than 0.
#' Smaller values suggest the probability distribution of the reduced model is near the full model
#' 
#' @examples
#' A1 <- solve(cov(neoOpen))
#' 
#' \dontrun{
#' A2 <- LoGo(neoOpen)
#' 
#' kld_value <- kld(A1, A2)
#' }
#' 
#' @references 
#' Kullback, S., & Leibler, R. A. (1951).
#' On information and sufficiency.
#' \emph{The Annals of Mathematical Statistics}, \emph{22}, 79-86.
#' doi: \href{https://doi.org/10.1214/aoms/1177729694}{10.1214/aoms/1177729694}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Kullback Leibler Divergence----
kld <- function (base, test)
{
    if(nrow(base)!=ncol(base))
    {base <- solve(cov(base))}
    
    if(nrow(test)!=ncol(test))
    {stop("Test must be an adjacency matrix")}
    
    kl1 <- sum(diag(solve(base)%*%test)) - log(det(solve(base)%*%test)) - ncol(base)
    kl2 <- sum(diag(solve(test)%*%base)) - log(det(solve(test)%*%base)) - ncol(test)
    
    kl <- log(kl1 + kl2)
    
    return(kl)
}
#----