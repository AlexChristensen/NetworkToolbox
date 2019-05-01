#' Latent Network Loadings
#'
#' @description Computes the within- and between-community
#' \code{\link[NetworkToolbox]{strength}} of each item
#' for each community
#'
#' @param A An adjacency matrix of network data
#'
#' @param comm A vector of community assignments
#' 
#' @param absolute Should network use absolute weights?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} for signed weights
#' 
#' @param rm.zero Should zeros be removed from the resulting matrix?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} to reduce the noise in the results
#' 
#' @return Returns a matrix of the unstandardized within- and between-community
#' strength values for each node
#' 
#' @details Simulation studies have demonstrated that a node's strength
#' centrality is roughly equivalent to factor loadings
#' (Christensen, Golino, & Silvia, 2019; Hallquist, Wright, & Molenaar, 2019).
#' Hallquist and colleagues (2019) found that node strength represented a
#' combination of dominant and cross-factor loadings. This function computes
#' each node's strength within each specified dimension, providing a rough
#' equivalent to factor loadings (including cross-loadings).
#'
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' wc <- louvain(A)$community
#' 
#' str.dim <- lnet.loads(A, wc, rm.zero = TRUE)
#' 
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com> and Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @export
#Strength by Dimension function
lnet.loads <- function(A, comm, absolute = TRUE,
                         rm.zero = FALSE)
{
    comc <- comcat(A,comm,metric="each",absolute=absolute)
    stab <- stable(A,comm,absolute=absolute)
    
    for(q in 1:nrow(comc))
    {comc[q,which(is.na(comc[q,]))] <- stab[q]}
    
    if(ncol(comc)!=1)
    {
        comm.str <- comc[,order(colnames(comc))]
        comm.str <- round(comm.str,3)
    }else{comm.str <- comc}
    
    if(rm.zero)
    {comm.str <- as.data.frame(ifelse(comm.str==0,"",comm.str))}
    
    return(comm.str)
}
#----