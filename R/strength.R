#' Node Strength
#' @description Computes strength of each node in a network
#' 
#' @param A An adjacency matrix of network data
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
#' #Undirected network
#' A <- TMFG(neoOpen)$A
#' 
#' str <- strength(A)
#' 
#' #Directed network
#' dep <- depend(neoOpen)
#' 
#' Adep <- TMFG(dep, depend = TRUE)$A
#' 
#' str <- strength(Adep)
#' 
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52} 1059-1069.
#' doi: \href{https://doi.org/10.1016/j.neuroimage.2009.10.003}{10.1016/j.neuroimage.2009.10.003}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Node Strength----
strength <- function (A)
{
    if(nrow(A)!=ncol(A))
    {stop("Input not an adjacency matrix")}
    
    A <- abs(A)
    A <- as.matrix(A)
    
    if(isSym(A)==TRUE)
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
            
            if(all(inStr==outStr))
            {Str <- round(as.vector(colSums(A)),2)
            names(Str) <- colnames(A)
            return(Str)
            }else{return(list(inStrength=inStr,outStrength=outStr,relInf=relinf))}
        }
}
#----