#' Stabilizing Nodes
#' @description Computes the within-community centrality for each node in the network
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
#' @param ... Additional arguments for \code{\link[igraph]{cluster_walktrap}}
#' and \code{\link[NetworkToolbox]{louvain}} community detection algorithms
#' 
#' @param diagonal Sets the diagonal values of the \code{A} input.
#' Defaults to \code{0}
#' 
#' @return A matrix containing the within-community centrality value for each node
#' 
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' stabilizing <- stable(A, comm = "walktrap")
#' 
#' @references 
#' Blanken, T. F., Deserno, M. K., Dalege, J., Borsboom, D., Blanken, P., Kerkhof, G. A., & Cramer, A. O. (2018).
#' The role of stabilizing and communicating symptoms given overlapping communities in psychopathology networks.
#' \emph{Scientific Reports}, \emph{8}, 5854.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Stabilizing----
#Updated 06.03.2020
stable <- function (A, comm = c("walktrap","louvain"),
                    cent = c("betweenness","rspbc","closeness",
                             "strength","degree","hybrid"), 
                    absolute = TRUE, diagonal = 0, ...)
{
    ###########################
    #### MISSING ARGUMENTS ####
    ###########################
    
    if(missing(comm))
    {comm <- "walktrap"
    }else{comm <- comm}
    
    if(missing(diagonal))
    {diagonal <- 0
    }else{diagonal <- diagonal}
    
    if(missing(cent))
    {cent <- "strength"
    }else{cent <- match.arg(cent)}
    
    #######################
    #### MAIN FUNCTION ####
    #######################
    
    # Set diagonal
    diag(A) <- diagonal
    
    # Make weights absolute
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
    
    # Check for names of nodes
    if(is.null(colnames(A)))
    {colnames(A) <- paste("V", 1:ncol(A), sep = "")}
    
    names(facts) <- colnames(A)
    
    # Unique communities
    uniq <- unique(facts)
    
    # Initialize community list
    fact <- list()
    
    # Loop through computing within-community centrality
    for(i in 1:length(uniq))
    {
        # Nodes only in community 'i'
        Ah <- A[which(facts == uniq[i]), which(facts == uniq[i])]
        
        # Centrality measure
        stab <- switch(cent,
                       betweenness = NetworkToolbox::betweenness(Ah),
                       rspbc = NetworkToolbox::rspbc(Ah),
                       closeness = NetworkToolbox::closeness(Ah),
                       strength = colSums(Ah),
                       degree = colSums(NetworkToolbox::binarize(Ah))
                       )
        
        # Input into list
        fact[[i]]<-stab
    }
    
    # Unlist for vector
    stabil <- unlist(fact)
    
    # Reorder to be consist with labels
    stabil <- stabil[names(facts)]
    
    # Check for missing values (change to 0)
    stabil <- ifelse(is.na(stabil),0,stabil)
    
    return(stabil)
}
#----
