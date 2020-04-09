#' Communicating Nodes
#' @description Computes the between-community strength for each node in the network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param comm Can be a vector of community assignments or community detection algorithms
#' (\code{"walktrap"} or \code{"louvain"}) can be used to determine the number of factors.
#' Defaults to \code{"walktrap"}.
#' Set to \code{"louvain"} for \code{\link[NetworkToolbox]{louvain}} community detection
#' 
#' @param cent Centrality measure to be used.
#' Defaults to \code{"strength"}.
#' 
#' @param absolute Should network use absolute weights?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} for signed weights
#' 
#' @param metric Whether the metric should be compute for across all of the communities
#' (a single value) or for each community (a value for each community).
#' Defaults to \code{"across"}.
#' Set to \code{"each"} for values for each community
#' 
#' @param diagonal Sets the diagonal values of the \code{A} input.
#' Defaults to \code{0}
#' 
#' @param ... Additional arguments for \code{\link[igraph]{cluster_walktrap}}
#' and \code{\link[NetworkToolbox]{louvain}} community detection algorithms
#' 
#' @return A vector containing the between-community strength value for each node
#' 
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' communicating <- comcat(A, comm = "walktrap", cent = "strength", metric = "across")
#' 
#' @references 
#' Blanken, T. F., Deserno, M. K., Dalege, J., Borsboom, D., Blanken, P., Kerkhof, G. A., & Cramer, A. O. (2018).
#' The role of stabilizing and communicating symptoms given overlapping communities in psychopathology networks.
#' \emph{Scientific Reports}, \emph{8}, 5854.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Communicating----
#Updated 18.03.2020
comcat <- function (A, comm = c("walktrap","louvain"),
                    cent = c("strength","degree"),
                    absolute = TRUE,
                    metric = c("across","each"),
                    diagonal = 0, ...)
{
    ###########################
    #### MISSING ARGUMENTS ####
    ###########################
    
    if(missing(comm))
    {comm <- "walktrap"
    }else{comm <- comm}
    
    if(missing(cent))
    {cent <- "strength"
    }else{cent <- match.arg(cent)}
    
    
    if(missing(metric))
    {metric <- "each"
    }else{metric <- match.arg(metric)}
    
    if(missing(diagonal))
    {diagonal <- 0
    }else{diagonal <- diagonal}
    
    #######################
    #### MAIN FUNCTION ####
    #######################
    
    # Set diagonal
    diag(A) <- diagonal
    
    # Change edges to absolute
    if(absolute)
    {A <- abs(A)}
    
    # Convert to communities
    if(any(eval(formals(NetworkToolbox::stable)$comm) %in% comm))
    {
        facts <- switch(comm,
                        walktrap = igraph::cluster_walktrap(NetworkToolbox::convert2igraph(A), ...)$membership,
                        louvain = igraph::cluster_louvain(NetworkToolbox::convert2igraph(A), ...)$membership
                        )
    }else{facts <- comm}
    
    # Convert facts to characters
    facts <- paste(facts)
    
    # Check for names of nodes
    if(is.null(colnames(A)))
    {colnames(A) <- paste("V", 1:ncol(A), sep = "")}
    
    names(facts) <- colnames(A)
    
    # Unique communities
    uniq <- unique(facts)
    
    # Initialize list
    fact <- list()
    
    # Check for unidimensionality
    if(length(na.omit(uniq)) != 1)
    {
        if(metric=="across") # All communities
        {
            
            # Loop over all communities not in target community
            for(i in 1:length(uniq))
            {
                # Connections outside of target community
                Ah <- A[which(facts != uniq[i]), which(facts == uniq[i])]
                
                # Centrality
                com <- switch(cent,
                              degree = colSums(NetworkToolbox::binarize(Ah)),
                              strength = colSums(Ah)
                              )
                
                # Input into list
                fact[[i]] <- com
            }
            
            # Unlist to vector
            commn <- unlist(fact)
            
            # Reorder to be consist with labels
            commat <- commn[names(facts)]
            
            # Label vector
            names(commat) <- colnames(A)
            
        }else if(metric=="each") # Individual communities
        {
            # Initialize item list
            item <- list()
            
            # Initialize matrix
            commat <- matrix(NA, nrow = nrow(A), ncol = length(uniq))
            
            # Add column names
            colnames(commat) <- paste(uniq)
            
            # Loop through each node
            for(i in 1:ncol(A))
            {
                # Connections for node 'i'
                Ah <- A[,i]
                
                # Communities outside of target community
                uniq.no <- uniq[which(uniq!=facts[i])]
                
                # Loop through each community
                for(j in 1:length(uniq.no))
                {
                    # Edges in target outside community
                    Aha <- Ah[which(facts==uniq.no[j])]
                    
                    # Centrality
                    com <- switch(cent,
                                  degree = sum(ifelse(Aha!=0,1,0)),
                                  strength = sum(Aha)
                                  )
                    
                    # Input into matrix
                    commat[i,paste(uniq.no[j])] <- com
                }
            }
            
            # Order and label matrix
            colnames(commat) <- uniq
            row.names(commat) <- colnames(A)
        }
    }else{
        # Initialize matrix
        commat <- as.matrix(rep(0,ncol(A)), ncol = 1)
        # Label matrix
        colnames(commat) <- paste(unique(comm))
        row.names(commat) <- colnames(A)
    }
    
    return(commat)
}
#----