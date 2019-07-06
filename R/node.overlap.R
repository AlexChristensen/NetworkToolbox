#' Detects Redundant Nodes in a Network
#' @description Identifies redundant nodes in the network based on several
#' measures. Computes the weighted topological overlap between
#' each node and every other node in the network. The weighted topological
#' overlap is implemented using the method from Nowick et al. (2009; see references)
#' and the function \link[wTO]{wTO} from the wTO package. Computes the zero-order
#' correlation between highly topological similar nodes. Computes the Jaccard similarity for the
#' overlap of the connections between two nodes. In sum, these indices offer an
#' array of similarity measures that can be used to determine node redundancy
#' 
#' @param data Data is necessary when \code{bootstrap = TRUE}.
#' Data is also necessary for zero-order and n.s. correlation proportions and
#' Jaccard similarity in the output.
#' Input data as a matrix or data frame
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param percentile Sets the number of overlapping nodes to show.
#' Defaults to \code{.99} or the nodes in the top 1\% of redundancy.
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
#' @param iter Number of iterations to perform in the bootstrap
#' 
#' @param network Network estimation method for bootstrap.
#' Defaults to \code{"glasso"}
#' 
#' @param ... Arguments to be passed onto \code{\link[qgraph]{EBICglasso}}
#' or \code{\link[NetworkToolbox]{TMFG}}
#' 
#' @return Produces a list containing:
#' 
#' \item{overlap}{A table with the nodes with the highest topological similarity
#' (\code{Topo. Similarity}; ordered from most to least).
#' 
#' If \code{data} is input, then the standard deviations of each item pair (i.e., 
#' \code{Node1} and \code{Node2}) are provided. Nodes with larger standard deviations
#' are preferred. Also includes the zero-order correlation between nodes (\code{Zero-order}), and
#' Jaccard similarity (\code{Jaccard}).}
#' 
#' \item{items}{A matrix containing the names of items that have high overlap}
#'  
#' @examples
#' 
#' net <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' result <- node.overlap(A = net, normal = FALSE)
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
#' @importFrom stats quantile
#' 
#' @export
#Node Overlap Function
node.overlap <- function(data = NULL, A, percentile = c(.90,.95,.99), method = c("abs","sign"),
                         bootstrap = FALSE, iter = NULL,
                         network = c("glasso","TMFG"), ...)
{
    #jaccard function
    jaccard <- function(mat)
    {
        bin <- ifelse(mat>0,1,0)
        
        num1 <- sum(ifelse(rowSums(bin)==2,1,0))
        num2 <- sum(ifelse(rowSums(bin)==0,1,0))
        
        num <- num1 + num2
        
        jac <- num/nrow(mat)
        
        return(jac)
    }
    
    #missing argument values
    if(!is.null(data))
    {data <- as.data.frame(data)}
    
    if(missing(A))
    {
        nodes <- ncol(data)
        n <- nrow(data)
    }else{nodes <- ncol(A)}
    
    if(missing(percentile))
    {percentile <- .99
    }else{percentile <- percentile[1]}
    
    if(missing(method))
    {method <- "sign"
    }else{method <- match.arg(method)}
    
    if(missing(network))
    {network <- "glasso"
    }else{network <- match.arg(network)}
    
    #bootstrap
    if(bootstrap)
    {
        #check for data
        if(is.null(data))
        {stop("data must be used for bootstrap")}
        
        #check for iterations
        if(is.null(iter))
        {
            iter <- 100
            message("Number of bootstrap iterations set to 100")
        }
        
        #initialize count
        count <- 0
        
        #initialize progress bar
        pb <- txtProgressBar(min=0,max=iter,style=3)
        
        #initialize data stores
        boot.net <- list()
        wto.mat <- matrix(NA,nrow=((nodes*nodes-nodes)/2),ncol=iter)
        sd.mat <- matrix(NA,nrow=nodes,ncol=iter)
        
        #repeat analysis
        repeat{
            #increase count
            count <- count+1
            
            #resample data
            rand.dat <- data[sample(n,n,replace=TRUE),]
            
            #correlations
            zo.cor <- cor(rand.dat)
            
            #apply network estimation method
            if(network=="glasso")
            {boot.net[[count]] <- invisible(qgraph::EBICglasso(zo.cor,n=n))
            }else if(network == "TMFG")
            {boot.net[[count]] <- TMFG(zo.cor,...)$A}
            
            #compute wTO
            tom <- wTO::wTO(boot.net[[count]],sign=method)
            diag(tom) <- 0
            
            #put wTO in matrix
            wto.mat[,count] <- tom[lower.tri(tom)]
            
            sd.mat[,count] <- apply(rand.dat,2,sd)
            
            #increase progress bar
            setTxtProgressBar(pb,count)
            
            #break when iterations is met
            if(count==iter)
            {break}
        }
        
        #close progress bar
        close(pb)
        
        #get means of bootstrap wTO
        lower <- apply(wto.mat,1,mean)
        sd.vec <- apply(sd.mat,1,mean)
        names(sd.vec) <- colnames(data)
    }else{
        
        #apply network estimation method
        if(missing(A))
        {
            if(!is.null(data))
            {
                #correlations
                zo.cor <- qgraph::cor_auto(data)
                
                if(network=="glasso")
                {A <- invisible(qgraph::EBICglasso(zo.cor,n=n,...))
                }else if(network == "TMFG")
                {A <- TMFG(zo.cor,...)$A}
            }else{stop("Argument 'A' or 'data' must have input")}
        }
        
        #compute weighted topological overlap
        tom <- wTO::wTO(A,sign=method)
        diag(tom) <- 0
        
        #lower triangle of TOM
        lower <- tom[lower.tri(tom)]
        
        #standard deviation vector
        if(!is.null(data))
        {
            sd.vec <- apply(data,2,sd)
            names(sd.vec) <- colnames(data)
        }
    }
    
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
    
    #topological overlap above certain percentile
    above <- quantile(lower,probs=percentile)
    
    #overlap matrix
    top.items <- which(lower>=above)
    
    res.list <- strsplit(names(top.items),split="--")
    
    if(!is.null(data))
    {
        overlap.items <- matrix(NA,nrow=length(top.items),ncol=4)
        
        #input sds
        for(i in 1:length(res.list))
        {
            overlap.items[i,3] <- round(sd.vec[res.list[[i]][1]],3)
            overlap.items[i,4] <- round(sd.vec[res.list[[i]][2]],3)
        }
    }else{overlap.items <- matrix(NA,nrow=length(top.items),ncol=2)}
    
    overlap.items[,1] <- names(top.items)
    overlap.items[,2] <- lower[top.items]
    
    #order by most overlap
    result <- overlap.items[order(overlap.items[,2],decreasing=TRUE),]
    
    if(!is.null(data))
    {
        if(!bootstrap)
        {
            #compute correlations between edges
            ##number of rows for result
            len <- nrow(result)

            #non-zero comparison vector
            #nz.comp.vec <- vector("numeric",length=len)
            #goldbricker#
            
            #zero-order correlations vector
            zo.vec <- vector("numeric",length=len)
            
            #jaccard vector
            jac.vec <- vector("numeric",length=len)
            
            #sd vectors
            sd.vecs <- matrix(NA,nrow=len,ncol=2)
            
            ##compare correlations
            for(l in 1:len)
            {
                #identify target nodes
                target.nodes <- match(res.list[[l]],colnames(A))
                
                #identify correlations that are not between them
                target.cor <- zo.cor[-target.nodes,target.nodes]
                target.net <- A[-target.nodes,target.nodes]
                
                #jaccard similarity
                target.jac <- jaccard(target.net)
                
                #identify correlation between them
                shared.cor <- zo.cor[target.nodes[1],target.nodes[2]]
                
                #goldbricker
                #number of rows for comparison
                #comp.len <- nrow(target.cor)
                
                #compute corrs
                #target.vec <- vector("numeric",length=comp.len)
                
                #for(r in 1:comp.len)
                #{target.vec[r] <- suppressWarnings(cocor::cocor.dep.groups.overlap(target.cor[r,1],target.cor[r,2],shared.cor, n = n)@steiger1980$p.value)}
                
                #number of non-significant non-zero relations
                #nz.comp.vec[l] <- round(sum(ifelse(target.vec<.05,0,1))/comp.len,3)
                #goldbricker
                
                #zero-order correlations
                zo.vec[l] <- round(shared.cor,3)
                
                #jaccard
                jac.vec[l] <- round(target.jac,3)
            }
            
            #tack on nz.comp.vec to results
            result <- cbind(result,zo.vec,jac.vec)
            
            #name columns
            colnames(result) <- c("Nodes", "Topo. Similarity", "sd.Node1", "sd.Node2", "Zero-order", "Jaccard")
        }else{colnames(result) <- c("Nodes", "Topo. Similarity", "sd.Node1", "sd.Node2")}
    }else{colnames(result) <- c("Nodes", "Topo. Similarity")}
    
    #input into list
    res <- list()
    res$overlap <- as.data.frame(result)
    res$nodes <- t(simplify2array(res.list))
    
    return(res)
    
    if(percentile==.99)
    {message("Lower percentile for more overlapping nodes")}
    
}
#----