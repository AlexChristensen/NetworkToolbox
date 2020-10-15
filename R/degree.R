#' Degree
#' @description Computes degree of each node in a network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @return A vector of degree values for each node in the network.
#' 
#' If directed network, returns a list containing:
#' 
#' \item{inDegree}{Degree of incoming edges (pointing to the node)}
#' 
#' \item{outDegree}{Degree of outgoing edges (pointing away from the node)}
#' 
#' \item{relInf}{Relative degree of incoming and outgoing edges.
#' Positive values indicate more outgoing degree relative to incoming degree.
#' Negative values indicate more incoming degree relative to outgoing degree}
#' 
#' @examples
#' #Undirected network
#' ## Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' deg <- degree(A)
#' 
#' #Directed network
#' \dontrun{
#' dep <- depend(neoOpen)
#' 
#' Adep <- TMFG(dep, depend = TRUE)$A
#' 
#' deg <- degree(Adep)
#' }
#' 
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{NeuroImage}, \emph{52}, 1059-1069.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Degree----
# Updated 18.05.2020
degree <- function (A)
{
    if(nrow(A)!=ncol(A))
    {stop("Input not an adjacency matrix")}
    
    A <- as.matrix(A)
    
    A <- binarize(A)
    
    if(isSymmetric(A, check.attributes = FALSE))
    {
        Deg <- as.vector(colSums(A))
        names(Deg) <- colnames(A)
        return(Deg)
    }else
    {
        #In-degree
        inDeg <- as.vector(colSums(A, na.rm = TRUE))
        names(inDeg) <- colnames(A)
        
        #Out-degree
        outDeg <- as.vector(rowSums(A, na.rm = TRUE))
        names(outDeg) <- colnames(A)
        
        #Relative influence
        relinf <- as.vector((outDeg-inDeg)/(outDeg+inDeg))
        names(relinf) <- colnames(A)
        
        #Reciprocal degree
        for(i in 1:nrow(A))
            for(j in 1:ncol(A))
            {
                A[i,j] <- ifelse(A[i,j] == 1 & A[j,i] == 1, 1, 0)
                A[j,i] <- ifelse(A[i,j] == 1 & A[j,i] == 1, 1, 0)
            }
        
        recipDeg <- colSums(A, na.rm = TRUE)
        names(recipDeg) <- colnames(A)
    
    if(all(relinf<.001))
    {Deg <- as.vector(inDeg)
    names(Deg) <- colnames(A)
    return(Deg)
    }else{return(list(inDegree=inDeg,
                      outDegree=outDeg,
                      recipDegree=recipDeg,
                      relInf=relinf))}
    }
}
#----