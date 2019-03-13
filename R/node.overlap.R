#' Detects Redundant Nodes in a Network
#' @description Computes the weighted topological overlap between
#' each node and every other node in the network. The weighted topological
#' overlap is implemented using the method from Nowick et al. (2009; see references)
#' and the function \link[wTO]{wTO} from the wTO package
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param percentile Sets the number of overlapping nodes to show.
#' Defaults to \code{.99} or the nodes in the top 1% of redundancy.
#' Set lower to identify more overlapping nodes
#' 
#' @param method Sets whether edges should be consider with absolute
#' or signed weights.
#' Defaults to \code{"sign"}
#' 
#' @param bootstrap Bootstraps with resampling to estimate more
#' reliable results.
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} to bootstrap
#' 
#' @param data Data is necessary when \code{bootstrap = TRUE}.
#' Input data as a matrix or data frame
#' 
#' @param iter Number of iterations to perform in the bootstrap
#' 
#' @param network Network estimation method for bootstrap.
#' Defaults to \code{"glasso"}
#' 
#' @param ... Additional arguments to be passed to the network estimation
#' method.
#' See \link[qgraph]{EBICglasso} and \link[NetworkToolbox]{TMFG} for more details
#' 
#' @return Produces a table containing the nodes with the highest
#' topological similarity (ordered from most to least)
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' result <- node.overlap(A)
#' 
#' @references 
#' Nowick, K., Gernat, T., Almaas, E., & Stubbs, L. (2009).
#' Differences in human and chimpanzee gene expression patterns define an evolving network of transcription factors in brain.
#' \emph{Proceedings of the National Academy of Sciences}, \emph{106}, 22358-22363.
#' doi: \href{https://doi.org/10.1073/pnas.0911376106}{10.1073/pnas.0911376106}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @importFrom stats quantile
#' 
#' @export
#Node Overlap Function
node.overlap <- function(A, percentile = c(.90,.95,.99), method = c("abs","sign"),
                         bootstrap = FALSE, data = NULL, iter = NULL,
                         network = c("glasso","TMFG"), ...)
{
    #number of nodes
    if(missing(A))
    {nodes <- ncol(data)
    }else{nodes <- ncol(A)}
    
    if(missing(percentile))
    {percentile <- .99}
    
    if(missing(method))
    {method <- match.arg(method)}
    
    if(bootstrap)
    {
        if(is.null(data))
        {stop("data must be used for bootstrap")}
        
        n <- nrow(data)
        
        if(is.null(iter))
        {stop("Number of iterations must be specified for bootstrap")}
        
        if(missing(network))
        {network <- "glasso"}
        
        count <- 1
        
        pb <- txtProgressBar(min=0,max=iter+1,style=3)
        
        boot.net <- list()
        wto.mat <- matrix(NA,nrow=((nodes*nodes-nodes)/2),ncol=iter)
        
        repeat{
            rand.dat <- data[sample(n,n,replace=TRUE),]
            
            if(network=="glasso")
            {boot.net[[count]] <- invisible(qgraph::EBICglasso(qgraph::cor_auto(rand.dat),n=n, ...))
            }else if(network == "TMFG")
            {boot.net[[count]] <- TMFG(rand.dat, ...)$A}
            
            tom <- wTO::wTO(boot.net[[count]],sign=method)
            diag(tom) <- 0
            
            wto.mat[,count] <- tom[lower.tri(tom)]
            
            count <- count+1
            
            setTxtProgressBar(pb,count)
            
            if(count==iter+1)
            {break}
        }
        
        close(pb)
        
        lower <- apply(wto.mat,1,mean)
    }else{
        
        #weighted topological overlap
        tom <- wTO::wTO(A,sign=method)
        diag(tom) <- 0
        
        #lower triangle of TOM
        lower <- tom[lower.tri(tom)]
    }
    
    #grab names for lower triangle
    name1 <- colnames(tom)
    name2 <- colnames(tom)
    
    name.mat <- tom
    
    for(i in 1:nodes)
        for(j in 1:nodes)
        {name.mat[i,j] <- paste(name1[j],name2[i],sep="--")}
    
    names(lower) <- name.mat[lower.tri(name.mat)]
    
    #topological overlap above certain percentile
    above <- quantile(lower,probs=percentile)
    
    #overlap matrix
    top.items <- which(lower>=above)
    overlap.items <- matrix(NA,nrow=length(top.items),ncol=2)
    
    overlap.items[,1] <- names(top.items)
    overlap.items[,2] <- lower[top.items]
    
    #order by most overlap
    result <- overlap.items[order(overlap.items[,2],decreasing=TRUE),]
    
    #name columns
    colnames(result) <- c("Nodes", "Topo. Similarity")
    
    #convert to data frame
    result <- as.data.frame(result)
    
    return(result)
    
    if(percentile==.99)
    {message("Lower percentile for more overlapping nodes")}
    
}
#----