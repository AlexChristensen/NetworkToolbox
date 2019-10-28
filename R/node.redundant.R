#' Detects Redundant Nodes in a Network
#' @description Identifies redundant nodes in the network based on several
#' measures. Computes the weighted topological overlap between
#' each node and every other node in the network. The weighted topological
#' overlap is implemented using the method from Nowick et al. (2009; see references)
#' and the function \link[wTO]{wTO} from the wTO package. 
#' 
#' @param A Matrix or data frame.
#' An adjacency matrix of network data (if \code{type = "wTO"}).
#' Dataset (if \code{type = "pcor"})
#' 
#' @param sig Numeric.
#' \emph{p}-value for significance of overlap (defaults to \code{.05}).
#' If more than 200 connections, then \code{\link[fdrtool]{fdrtool}}
#' is used to correct for false positives. In these instances, \code{sig}
#' sets the \emph{q}-value for significance of overlap (defaults to \code{.10})
#' 
#' @param type Character.
#' Computes weighted topological overlap (\code{"wTO"})
#' or partial correlations (\code{"pcor"})
#' 
#' @param method Character.
#' Computes significance using the standard \emph{p}-value (\code{"alpha"}),
#' bonferonni corrected \emph{p}-value (\code{"bonferroni"}),
#' false-discovery rate corrected \emph{p}-value (\code{"FDR"}),
#' or adaptive alpha \emph{p}-value (\code{\link[NetworkToolbox]{adapt.a}}).
#' Defaults to \code{"alpha"}
#' 
#' @return Returns a list with vectors nested within the list corresponding
#' to redundant nodes with the name of object in the list
#' 
#' @examples
#' # normal set to FALSE for CRAN tests
#' net <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' # weighted topological overlap
#' result <- node.redundant(A = net, method = "adapt", type = "wTO")
#' 
#' # partial correlation
#' result <- node.redundant(A = neoOpen, method = "adapt", type = "pcor")
#' 
#' @references
#' #wTO
#' Nowick, K., Gernat, T., Almaas, E., & Stubbs, L. (2009).
#' Differences in human and chimpanzee gene expression patterns define an evolving network of transcription factors in brain.
#' \emph{Proceedings of the National Academy of Sciences}, \emph{106}, 22358-22363.
#' doi: \href{https://doi.org/10.1073/pnas.0911376106}{10.1073/pnas.0911376106}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @importFrom stats pgamma pnorm
#' 
#' @export
#Redundant Nodes Function
node.redundant <- function (A, sig, type = c("wTO", "pcor"),
                            method = c("alpha", "bonferroni", "FDR", "adapt"))
{
    if(missing(method))
    {method <- "adapt"
    }else{method <- match.arg(method)}
    
    if(missing("type"))
    {type <- "pcor"
    }else{method <- match.arg(type)}
    
    #number of nodes
    nodes <- ncol(A)
    
    #compute type of overlap method
    if(type == "wTO")
    {tom <- wTO::wTO(A,sign="sign")
    }else if(type == "pcor")
    {
        tom <- qgraph::cor_auto(A)
        tom <- -cov2cor(solve(tom))
    }
    
    diag(tom) <- 0
    
    #lower triangle of TOM
    lower <- tom[lower.tri(tom)]
    
    #grab names for lower triangle
    name1 <- colnames(tom)
    name2 <- colnames(tom)
    
    #initialize name matrix
    name.mat <- tom
    
    #produce name matrix
    for(i in 1:nodes)
        for(j in 1:nodes)
        {name.mat[i,j] <- paste(name1[j],name2[i],sep="--")}
    
    #grab lower triangle names
    names(lower) <- name.mat[lower.tri(name.mat)]
    
    #obtain only positive values
    pos.vals <- na.omit(ifelse(lower<=0,NA,lower))
    attr(pos.vals, "na.action") <- NULL
    
    #determine distribution
    ##distributions
    distr <- c("norm","gamma")
    aic <- vector("numeric", length = length(distr))
    names(aic) <- c("normal","gamma")
    
    for(i in 1:length(distr))
    {aic[i] <- fitdistrplus::fitdist(pos.vals,distr[i],method="mle")$aic}
    
    #obtain distribution parameters
    g.dist <- suppressWarnings(MASS::fitdistr(pos.vals, names(aic)[which.min(aic)]))
    
    #compute significance values
    pval <- switch(names(aic)[which.min(aic)],
           
                  normal = 1 - unlist(lapply(pos.vals, #positive wTo
                                      pnorm, #probability in normal distribution
                                      mean = g.dist$estimate["mean"], #mean of normal
                                      sd = g.dist$estimate["sd"]) #sd of normal
                                      ),
                  
                  gamma = 1 - unlist(lapply(pos.vals, #positive wTo
                                             pgamma, #probability in gamma distribution
                                             shape = g.dist$estimate["shape"], #shape of gamma
                                             rate = g.dist$estimate["rate"]) #rate of gamma
                                     ),
                  )
    
    #compute false-discvoery rate
    if(method == "FDR")
    {
        sig <- ifelse(missing(sig),.10,sig)
        
        pval <- suppressWarnings(fdrtool::fdrtool(pval, statistic = "pvalue", plot = FALSE,verbose = FALSE)$qval)
    }else{
        
        sig <- ifelse(missing(sig),.05,sig)
        
        if(method == "bonferroni")
        {sig <- sig / length(pos.vals)
        }else if(method == "adapt")
        {sig <- adapt.a("cor", alpha = sig, n = length(pos.vals), efxize = "medium")$adapt.a}
    }
    
    #identify q-values less than sigificance
    res <- pos.vals[which(pval<=sig)]
    
    #if no redundant, then print no redundant nodes
    if(length(res)==0)
    {
        message(paste("No redundant nodes identified. Increase 'sig' arugment to detect more nodes."))
        res.list <- NA
    }else{
        
        #create result matrix
        split.res <- unlist(strsplit(names(res), split = "--"))
        res.mat <- t(simplify2array(sapply(names(res), strsplit, split = "--")))
        
        #initialize result list
        res.list <- list()
        
        #initialize count
        count <- 0
        
        while(nrow(res.mat)!=0)
        {
            #increase count
            count <- count + 1
            
            #get variable counts
            var.counts <- sort(table(split.res), decreasing = TRUE)
            
            if(!all(var.counts==1))
            {
                #identify targets
                target <- which(res.mat == names(var.counts[1]), arr.ind = TRUE)[,"row"]
                
                #insert values into list
                res.list[[names(var.counts[1])]] <- setdiff(unlist(strsplit(names(target),split="--")),names(var.counts[1]))
                
                #remove rows from result matrix
                res.mat <- res.mat[-target,]
                
                #remove variables from split result
                split.res <- as.vector(res.mat)
                
                #force matrix
                if(is.vector(res.mat))
                {res.mat <- t(as.matrix(res.mat))}
                
            }else
            {
                for(i in 1:nrow(res.mat))
                {res.list[[res.mat[i,1]]] <- unname(res.mat[i,2])}
                
                res.mat <- res.mat[-c(1:nrow(res.mat)),]
            }
        }
    }
    
    return(res.list)
}
#----