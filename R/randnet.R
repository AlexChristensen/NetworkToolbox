#' Generates a Random Network
#' 
#' @description Generates a random binary network
#' 
#' @param nodes Numeric.
#' Number of nodes in random network
#' 
#' @param edges Numeric.
#' Number of edges in random network
#' 
#' @param A Matrix or data frame.
#' An adjacency matrix (i.e., network) to be used to estimated a random network with
#' fixed edges (allows for asymmetric network estimation)
#' 
#' @return Returns an adjacency matrix of a random network
#' 
#' @examples
#' rand <- randnet(10, 27)
#' 
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{NeuroImage}, \emph{52}, 1059-1069.
#' 
#' Csardi, G., & Nepusz, T. (2006).
#' The \emph{igraph} software package for complex network research.
#' \emph{InterJournal, Complex Systems}, 1695.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
# Random Network----
# Updated 20.04.2020
randnet <- function (nodes = NULL, edges = NULL, A = NULL)
{
    if(is.null(A))
    {
        # Initialize matrix
        mat <- matrix(1, nrow = nodes, ncol = nodes)
        
        # Set diagonal to zero
        diag(mat) <- 0
        
        # Indices of upper diagonal
        ind <- ifelse(upper.tri(mat) == TRUE, 1, 0)
        i <- which(ind == 1)
        
        # Sample zeros and ones
        rp <- sample(length(i))
        # Get indices
        irp <- i[rp]
        
        # Initialize random matrix
        rand <- matrix(0, nrow = nodes, ncol = nodes)
        
        # Insert edges
        rand[irp[1:edges]] <- 1
        
        # Make symmetric
        rand <- rand + t(rand)
        
    }else{
        
        # Make diagonal of network zero
        diag(A) <- 0
        
        # Compute degree
        degrees <- NetworkToolbox::degree(A)
        
        # Get degrees based on directed or undirected
        # Use igraph
        if(is.list(degrees))
        {rand <- as.matrix(igraph::as_adj(igraph::sample_degseq(out.deg = degrees$outDegree, in.deg = degrees$inDegree, method = "vl")))
        }else{rand <- as.matrix(igraph::as_adj(igraph::sample_degseq(out.deg = degrees, method = "vl")))}
    }
    
    return(rand)
}
#----