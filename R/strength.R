#' Node Strength
#' @description Computes strength of each node in a network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param absolute Should network use absolute weights?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} for signed weights
#' 
#' @return A vector of strength values for each node in the network.
#' 
#' If directed network, returns a list containing:
#' 
#' \item{inStrength}{Strength of incoming edges (pointing to the node)}
#' 
#' \item{outStrength}{Strength of outgoing edges (pointing away from the node)}
#' 
#' \item{relInf}{Relative degree of incoming and outgoing edges.
#' Positive values indicate more outgoing strength relative to incoming strength.
#' Negative values indicate more incoming strength relative to outgoing strength}
#' 
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' str <- strength(A)
#' 
#' #Directed network
#' \dontrun{
#' dep <- depend(neoOpen)
#' 
#' Adep <- TMFG(dep, depend = TRUE)$A
#' 
#' str <- strength(Adep)
#' }
#' 
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{NeuroImage}, \emph{52} 1059-1069.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Node Strength----
strength <- function (A, absolute = TRUE)
{
    if(is.vector(A))
    {return(0)
    }else if(nrow(A)!=ncol(A))
    {stop("Input not an adjacency matrix")}
    
    if(absolute)
    {A <- abs(A)}
    A <- as.matrix(A)
    
    if(isSymmetric(A, check.attributes = FALSE))
    {
        Str <- round(as.vector(colSums(A)),2)
        names(Str) <- colnames(A)
        return(Str)
    }else{
        #In-strength
        inStr <- as.vector(colSums(A))
        names(inStr) <- colnames(A)
        #Out-strength
        outStr <- as.vector(rowSums(A))
        names(outStr) <- colnames(A)
        #Relative influence
        relinf <- as.vector((outStr-inStr)/(outStr+inStr))
        names(relinf) <- colnames(A)
            
            if(all(relinf<.001))
            {Str <- round(as.vector(colSums(A)),2)
            names(Str) <- colnames(A)
            return(Str)
            }else{return(list(inStrength=inStr,outStrength=outStr,relInf=relinf))}
        }
}
#----