#' Network Adjusted Mean/Sum
#' @description The hybrid centrality is used to adjust the mean or sum score of participant's community scores based on each node's centrality.
#' Each participant's response values are multipled by the corresponding hybrid centrality value (uses "random" for BC argument).
#' In this way, more central nodes contribute a greater score and less central nodes contribute a lesser score
#' 
#' @param data Must be a dataset
#' 
#' @param A Adjacency matrix that has already been filtered
#' 
#' @param adjusted Should adjusted values be the mean or sum score?
#' Defaults to "mean".
#' Set to "sum" for sum scores
#' 
#' @param comm Can be a vector of community assignments or community detection algorithms
#' ("walktrap" or "louvain") can be used to determine the number of communities.
#' Defaults to 1 community.
#' Set to "walktrap" for the walktrap algortihm.
#' Set to "louvain" for louvain community detection
#' 
#' @param ... Additional arguments for community detection algorithms
#' 
#' @return Returns a list containing:
#' 
#' \item{Standardized}{The standardized network adjusted score for each participant}
#' 
#' \item{Unstandardized}{The unstandardized network adjusted score (mean or sum) for each participant}
#' 
#' \item{CommItems}{The items associated with the specified or identified communities}
#' 
#' \item{CommCor}{Correlations between the specified or identified communities}
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' #one community
#' sumadj <- nams(neoOpen, A, adjusted = "sum")
#' 
#' #theoretical communities
#' knowncomm <- nams(neoOpen, A,
#' comm = c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,8),rep(6,8)))
#' 
#' #walktrap communities
#' walkadj <- nams(neoOpen, A, adjusted = "sum", comm = "walktrap")
#' 
#' @references
#' Christensen, A. P. (2018).
#' NetworkToolbox: Methods and measures for brain, cognitive, and psychometric network analysis in R.
#' \emph{PsyArXiv}.
#' doi: \href{https://doi.org/10.31234/osf.io/6kmav}{10.31234/osf.io/6kmav}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Network Adjusted Mean/Sum----
nams <- function (data, A,
                  adjusted = c("mean","sum"),
                  comm = c("walktrap","louvain"), ...)
{
    if(missing(data))
    {stop("Data is required for analysis")}
    
    if(missing(A))
    {stop("Adjacency matrix is required for analysis")}
    
    if(missing(adjusted))
    {adjusted<-"mean"
    }else(adjusted<-match.arg(adjusted))
    
    if(missing(comm))
    {comm<-rep(1,ncol(data))}
    
    #FOR FUTURE DEVELOPMENT
    #if(missing(na.data))
    #{na.data<-"none"
    #}else{na.data<-match.arg(na.data)}
    
    if(!is.numeric(comm))
    {
        if(length(comm)!=ncol(data))
        {
            if(comm=="walktrap")
            {facts<-igraph::walktrap.community(convert2igraph(A))$membership
            }else if(comm=="louvain")
            {facts<-louvain(A,...)$community}
        }else{facts <- comm}
    }else{facts<-comm}
    
    uniq <- unique(facts)
    len <- length(uniq)
    
    if(len>1)
    {fact<-matrix(0,nrow=nrow(data),ncol=(len+1))
    }else{fact<-matrix(0,nrow=nrow(data),ncol=len)}
    
    fullhyb <- hybrid(A, BC="random")
    
    if(len>1)
    {commLCu <- comm.close(A, facts, weighted = FALSE)}
    
    for(i in 1:len)
    {
        #old algorithm
        #Ah <- A[which(facts==uniq[i]),which(facts==uniq[i])]
        #hyb <- hybrid(Ah, BC="random")
        
        fh <- fullhyb[which(facts==uniq[i])]
        
        if(len>1)
        {
            ch <- commLCu[which(names(commLCu)==uniq[i])]
            comb <- fh + ch
        }else{comb <- fh}
        
        tdata<-t(data[,which(facts==uniq[i])])
        
        mat<-matrix(0,nrow=nrow(tdata),ncol=ncol(tdata))
        
        for(j in 1:ncol(tdata))
        {mat[,j]<-comb*tdata[,j]}
        
        if(adjusted=="mean")
        {wei<-colMeans(mat)
        }else if(adjusted=="sum")
        {wei<-colSums(mat)}
        
        adj<-wei/mean(comb)
        
        fact[,i]<-adj
    }
    
    if(len>1)
    {
        tdata<-t(data)
        
        mat<-matrix(0,nrow=nrow(tdata),ncol=ncol(tdata))
        
        for(j in 1:ncol(tdata))
        {mat[,j] <- fullhyb*tdata[,j]}
        
        if(adjusted=="mean")
        {wei<-colMeans(mat)
        }else if(adjusted=="sum")
        {wei<-colSums(mat)}
        
        adj <- wei/mean(fullhyb)
        
        fact[,(len+1)] <- adj
        
        for(i in 1:nrow(fact))
            for(j in (ncol(fact)-1))
            {
                if(adjusted=="mean")
                {
                    diff <- fact[i,(j+1)] - mean(fact[i,1:j])
                    chan <- diff/j
                    fact[i,1:j] <- fact[i,1:j] + chan
                }else if(adjusted=="sum")
                {
                    diff <- fact[i,(j+1)] - sum(fact[i,1:j])
                    chan <- diff/j
                    fact[i,1:j] <- fact[i,1:j] + chan
                }
            }
    }
    
    std_fact <- scale(fact)
    fact <- as.data.frame(fact)
    std_fact <- as.data.frame(std_fact)
    
    
    for(l in 1:nrow(data))
    {row.names(fact)[l]<-paste("Part",l,sep="")}
    
    if(len>1)
    {
        if(!is.character(comm))
        {
            for(k in 1:len)
            colnames(fact)[k]<-paste("Community",uniq[k],sep="")
        }else{colnames(fact)[1:len] <- uniq}
        
        colnames(fact)[len+1]<-"overall"
        
        corr <- cor(fact)
        colnames(corr) <- colnames(fact)
        row.names(corr) <- colnames(corr)
    }else{
        colnames(fact)<-"overall"
        corr <- 1
    }
    
    colnames(std_fact) <- colnames(fact)
    row.names(std_fact) <- row.names(fact)
    
    matf<-matrix(0,nrow=len,ncol=2)
    
    for(i in 1:len)
    {matf[i,2]<-paste(colnames(data)[which(facts==uniq[i])],collapse = ", ")}
    
    matf[,1] <- uniq
    
    matf<-as.data.frame(matf)
    colnames(matf)<-c("Community","Items")
    
    return(list(Standardized=std_fact,Unstandardized=fact,CommItems=matf,CommCor=corr))
}
#----