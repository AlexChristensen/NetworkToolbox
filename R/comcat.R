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
#' doi: \href{https://doi.org/10.1038/s41598-018-24224-2}{10.1038/s41598-018-24224-2}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Communicating----
comcat <- function (A, comm = c("walktrap","louvain"),
                    cent = c("strength","degree"),
                    absolute = TRUE,
                    metric = c("across","each"),
                    diagonal = 0, ...)
{
    #nodes
    n <- ncol(A)
    
    #change diagonal values if necessary
    if(missing(diagonal))
    {diagonal <- 0
    }else{diagonal <- diagonal}
    
    if(missing(cent))
    {cent<-"strength"
    }else{cent<-match.arg(cent)}
    
    #set communities
    if(missing(comm))
    {comm<-"walktrap"
    }else{comm<-comm}
    
    if(missing(metric))
    {metric <- "across"
    }else{metric <- match.arg(metric)}
    
    if(!is.numeric(comm))
    {
        if(comm=="walktrap")
        {facts<-igraph::walktrap.community(convert2igraph(A), ...)$membership
        }else if(comm=="louvain")
        {facts<-louvain(A,...)$community}
    }else{facts<-comm}
    
    keepord<-cbind(rep(1:ncol(A)),facts)
    ord<-keepord[order(keepord[,2]),]
    
    fact<-list()
    
    if(length(unique(comm))!=1)
    {
        if(metric=="across")
        {
            
            for(i in 1:max(facts, na.rm = TRUE))
            {
                Ah <- A[which(facts!=i),which(facts==i)]
                
                if(cent=="degree")
                {com<-colSums(binarize(Ah))
                }else if(cent=="strength")
                {com<-colSums(Ah,absolute)}
                
                fact[[i]]<-com
            }
            
            commn<-unlist(fact)
            
            bind<-cbind(ord,commn)
            
            commord<-bind[order(bind[,1]),]
            
            commmat<-matrix(commord[,3],nrow=nrow(commord),ncol=1)
            
            commmat <- as.vector(commmat)
            
            names(commmat) <- colnames(A)
            
            return(commmat)
        }else if(metric=="each")
        {
            
            uniq <- unique(facts)
            
            item <- list()
            
            commat <- matrix(NA,nrow=nrow(A),ncol=length(uniq))
            
            colnames(commat) <- paste(uniq)
            
            for(i in 1:ncol(A))
            {
                Ah <- A[,i]
                
                uniq.no <- uniq[which(uniq!=facts[i])]
                
                for(j in 1:length(uniq.no))
                {
                    Aha <- Ah[which(facts==uniq.no[j])]
                    
                    if(cent=="degree")
                    {com<-sum(ifelse(Aha!=0,1,0))
                    }else if(cent=="strength")
                    {com<-sum(Aha)}
                    
                    commat[i,paste(uniq.no[j])] <- com
                }
            }
            
            commat[,paste(uniq[order(uniq)])]
            row.names(commat) <- colnames(A)
            
            commat <- commat[,order(colnames(commat))]
            
            return(commat)
        }
    }else{
        commat <- as.matrix(rep(0,ncol(A)),ncol=1)
        colnames(commat) <- paste(unique(comm))
        
        return(commat)
    }
}
#----