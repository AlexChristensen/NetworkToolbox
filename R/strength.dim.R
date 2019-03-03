#' Node Strength by Dimension
#'
#' @description Computes the within- and between-community strength of each item
#' for each community
#'
#' @param A An adjacency matrix of network data
#'
#' @param comm A vector of community assignments
#' 
#' @param rm.zero Should zeros be removed from the resulting matrix?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} to reduce the noise in the results
#' 
#' @return Returns a matrix of the unstandardized within- and between-community
#' strength values for each node
#'
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' wc <- louvain(A)$community
#' 
#' str.dim <- strength.dim(A, wc, rm.zero = TRUE)
#' 
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com> and Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @export
#Strength by Dimension function
strength.dim <- function(A, comm, rm.zero = FALSE)
{
    comc <- comcat(A,comm,metric="each")
    stab <- stable(A,comm)
    
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