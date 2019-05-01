#' Latent Network Scores
#' @description This function computes latent network scores for
#' factor analysis models. Latent scores are computed based on 
#' each node's \code{\link[NetworkToolbox]{strength}} within each
#' community (i.e., factor) in the network. These values are used
#' as network "factor loadings" for the weights of each item. Notably,
#' network analysis allows nodes to load onto more than one factor.
#' These loadings are considered in the factor scores. In addition,
#' if the construct is a hierarchy (e.g., personality questionnaire;
#' items in facet scales in a trait domain), then an overall latent
#' score can be computed (see argument \code{general}). These overall
#' latent scores are computed using \code{\link[NetworkToolbox]{comm.close}}
#' as weights, which are roughly similar to general factor loadings in a
#' CFA model (see Christensen, Golino, & Silvia, 2019). The latent score
#' estimates are roughly equivalent to the Maximum Likelihood method in
#' \code{lavaan}'s \code{\link[lavaan]{cfa}} function. An important difference
#' is that the latent network scores account for cross-loadings in their
#' estimation of latent scores.
#' 
#' @param data Must be a dataset
#' 
#' @param A Network. An adjacency matrix
#' 
#' @param comm Vector. Can be a vector of community assignments or community detection algorithms
#' (\code{"walktrap"} or \code{"louvain"}) can be used to determine the number of communities.
#' Defaults to \code{1} community.
#' Set to \code{"walktrap"} for the \code{\link[igraph]{cluster_walktrap}} algortihm.
#' Set to \code{"louvain"} for \code{\link[NetworkToolbox]{louvain}} community detection
#' 
#' @param ... Additional arguments for \code{\link[igraph]{cluster_walktrap}}
#' and \code{\link[NetworkToolbox]{louvain}} community detection algorithms
#' 
#' @return Returns a list containing:
#' 
#' \item{scores}{The standardized latent network scores for each participant
#' and community (including the overall score)}
#' 
#' \item{commCor}{Correlations between the specified or identified communities}
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' #one community
#' sumadj <- lnet.score(neoOpen, A)
#' 
#' #theoretical communities
#' knowncomm <- lnet.score(neoOpen, A,
#' comm = c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,8),rep(6,8)))
#' 
#' #walktrap communities
#' walkadj <- lnet.score(neoOpen, A, comm = "walktrap")
#' 
#' @references
#' Christensen, A. P. (2018).
#' NetworkToolbox: Methods and measures for brain, cognitive, and psychometric network analysis in R.
#' \emph{The R Journal}, \emph{10}, 422-439.
#' doi: \href{https://doi.org/10.32614/RJ-2018-065}{10.32614/RJ-2018-065}
#' 
#' Christensen, A. P., Golino, H. F., & Silvia, P. J. (2019).
#' A psychometric network perspective on the measurement and assessment of personality traits.
#' \emph{PsyArXiv}.
#' doi: \href{https://doi.org/10.31234/osf.io/ktejp}{10.31234/osf.io/ktejp}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @importFrom utils packageDescription
#' 
#' @export
#Network Adjusted Mean/Sum----
lnet.score <- function (data, A,
                  comm = c("walktrap","louvain"), ...)
{
    ####Missing arguments checks####
    if(missing(data))
    {stop("Data is required for analysis")}
    
    if(missing(A))
    {stop("Adjacency matrix is required for analysis")}
    
    #Default to single latent variable
    if(missing(comm))
    {comm <- rep(1,ncol(data))}
    ####Missing arguments checks####
    
    #Determine factors
    if(length(comm)==1)
    {
        if(comm=="walktrap") #walktrap
        {facts <- igraph::walktrap.community(convert2igraph(A),...)$membership
        }else if(comm=="louvain") #louvain
        {facts <- louvain(A,...)$community}
    }else{facts <- comm}
    
    #Number of factors
    nfacts <- length(unique(facts))
    
    #Initialize factor result matrix
    if(nfacts > 1)
    {fact.res <- as.data.frame(matrix(0, nrow = nrow(data), ncol = (nfacts + 1)))
    }else{fact.res <- as.data.frame(matrix(0, nrow = nrow(data), ncol = nfacts))}
    
    #Compute network loadings
    if(nfacts > 1)
    {
        P <- lnet.loads(A, comm = facts, absolute = FALSE)
        #Standardize
        P <- t(t(P) / sqrt(colSums(P)))
    }else{
        P <- strength(A, absolute = FALSE)
        #Standardize
        P <- t(t(P) / sqrt(sum(P)))
    }
    
    ####LATENT SCORE FUNCTION####
    lv.score <- function(loads, data)
    {
        #Initialize participant latent scores
        lat.sco <- matrix(0, nrow = nrow(data), ncol = ncol(loads))
        
        #Compute latent factor scores (ML)
        for(i in 1:ncol(loads))
        {
            #Network loadings for each factor
            f.load <- loads[which(loads[,i]!=0),i]
            
            #Grab items associated with factor
            dat <- data[,names(f.load)]
            
            #Grab std dev of items associated with factor 
            f.sds <- apply(dat,2,sd)
            
            #Obtain relative weights
            rel <- f.load / f.sds
            rel.wei <- rel / sum(rel)
            
            #Compute latent scores
            lat.sco[,i] <- as.vector(rowSums(t(t(dat) * rel.wei)))
        }
        
        colnames(lat.sco) <- colnames(loads)
        
        return(lat.sco)
    }
    ####LATENT SCORE FUNCTION####
    
    #Populate factor result matrix
    lv.sco <- lv.score(P, data)
    fact.res[,1:nfacts] <- lv.sco
    
    if(nfacts > 1)
    {colnames(fact.res)[1:nfacts] <- colnames(P)
    }else{colnames(fact.res) <- "1"}

    #Compute correlations between latent factors
    C <- cor(lv.sco)
    
    if(nfacts > 1)
    {
        #Compute general network loadings
        Pg <- comm.close(A = A, comm = facts)
        
        #Overall score
        G <- rowSums(t(t(lv.sco) * Pg))
        fact.res[,(nfacts + 1)] <- G
        colnames(fact.res)[nfacts + 1] <- "Overall"
        
        #Re-compute correlations between latent factors
        C <- cor(fact.res)
    }
    
    #Results
    res <- list()
    res$scores <- round(apply(fact.res,2,scale),3)
    res$commCor <- C
    
    return(res)
}
#----