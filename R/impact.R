#' Node Impact
#' @description Computes impact measure or how much the average distance in the
#' network changes with that node removed of each node in a network
#' (\strong{Please see and cite Kenett et al., 2011})
#' 
#' @param A An adjacency matrix of network data
#' 
#' @return A vector of node impact values for each node in the network
#' (impact > 0, greater ASPL when node is removed; impact < 0,
#' lower ASPL when node is removed)
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' nodeimpact <- impact(A)
#' 
#' @references 
#' Cotter, K. N., Christensen, A. P., & Silvia, P. J. (in press).
#' Understanding inner music: A dimensional approach to musical imagery.
#' \emph{Psychology of Aesthetics, Creativity, and the Arts}.
#' 
#' Kenett, Y. N., Kenett, D. Y., Ben-Jacob, E., & Faust, M. (2011).
#' Global and local features of semantic networks: Evidence from the Hebrew mental lexicon.
#' \emph{PloS one}, \emph{6}, e23912.
#' doi: \href{https://doi.org/10.1371/journal.pone.0023912}{10.1371/journal.pone.0023912}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Node Impact----
impact <- function (A)
{
        allP <- pathlengths(A)$ASPL
        remove <- matrix(0,nrow=nrow(A),ncol=1)
        
        pb <- txtProgressBar(max=ncol(A), style = 3)
        
        for(j in 1:ncol(A))
        {
            remove[j,]<-(pathlengths(A[-j,-j])$ASPL)-allP
            
            setTxtProgressBar(pb, j)
        }
        
        close(pb)
        
        remove <- as.vector(round(remove,3))
        names(remove) <- colnames(A)
        
    return(remove)
}
#----