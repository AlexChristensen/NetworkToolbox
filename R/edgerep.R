#' Edge Replication
#' @description Computes the number of edges that replicate between two cross-sectional networks
#' 
#' @param A An adjacency matrix of network A
#' 
#' @param B An adjacency matrix of network B
#' 
#' @param corr Correlation method for assessing the relationship between the replicated edge weights.
#' Defaults to \code{"pearson"}.
#' Set to \code{"spearman"} for non-linear or monotonic associations.
#' Set to \code{"kendall"} for rank-order correlations
#' 
#' @return Returns a list containing:
#' 
#' \item{replicatedEdges}{The edges that replicated and their weights}
#' 
#' \item{replicated}{Number of edges that replicated}
#' 
#' \item{meanDiff}{The average edge weight difference between the edges that replicated}
#' 
#' \item{sdDiff}{The standard deviation edge weight difference between the edges that replicated}
#' 
#' \item{cor}{The correlation between the edges that replicated}
#' 
#' Lists for each network contain:
#' 
#' \item{totalEdges}{Total possible number of edges to be replicated}
#' 
#' \item{percentage}{Percentage of edges that replicated relative to total possible}
#' 
#' \item{density}{The density of the network}
#' 
#' @examples
#' # normal set to FALSE for CRAN tests
#' tmfg <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' # normal set to FALSE for CRAN tests
#' mast <- MaST(neoOpen, normal = FALSE)
#' 
#' edges <- edgerep(tmfg, mast)
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Edge Replication----
edgerep <- function (A, B,
                     corr = c("pearson", "spearman", "kendall"))
{
    
    nameA <- deparse(substitute(A))
    nameB <- deparse(substitute(B))
    
    A <- as.matrix(A)
    B <- as.matrix(B)
    
    if(missing(corr))
    {corr<-"pearson"
    }else{corr<-match.arg(corr)}

    if(!isSymmetric(A, check.attributes = FALSE))
    {
        if(all(rowSums(A)==colSums(A)))
        {A[lower.tri(A)] <- A[upper.tri(A)]
        }else{A<-A+t(A)
        warning(paste("Adjacency matrix",nameA,"was made to be symmetric using upper triangle"))}
    }
    
    if(!isSymmetric(B, check.attributes = FALSE))
    {
        if(all(rowSums(B)==colSums(B)))
        {B[lower.tri(B)] <- B[upper.tri(B)]
        }else{B<-B+t(B)
        warning(paste("Adjacency matrix",nameB,"was made to be symmetric using upper triangle"))}
    }
    
    n<-ncol(A)
    
    diag(A)<-0
    diag(B)<-0
    
    repedges<-which(A&B!=0)
    repmat<-matrix(0,nrow=n,ncol=n)
    repmat[repedges]<-1
    count<-sum(repmat)/2
    
    repmat[lower.tri(repmat)]<-A[lower.tri(A)]*repmat[lower.tri(repmat)]
    repmat[upper.tri(repmat)]<-B[upper.tri(B)]*repmat[upper.tri(repmat)]
    
    i<-sort(c(rep(1:n,n)))
    j<-c(rep(1:n,n))
    lt<-as.vector(A)
    ut<-as.vector(B)
    
    replist<-cbind(i,j,lt,ut)
    replist<-replist[which(replist[,3]&replist[,4]!=0),]
    ord<-replist[order(replist[,3]),]
    sinlist<-ord[seq(1,nrow(ord),2),]
    sinlist<-sinlist[order(sinlist[,1],sinlist[,2]),]
    
    colnames(sinlist)<-c("nodeTo","nodeFrom",paste("weight in",nameA),paste("weight in",nameB))
    
    if(!is.null(colnames(A)))
    {
        if(!is.null(colnames(B)))
        {
            for(i in 1:nrow(sinlist))
                for(j in 1:length(colnames(A)))
                {
                    if(sinlist[i,1]==j)
                    {sinlist[i,1]<-colnames(A)[j]}
                    if(sinlist[i,2]==j)
                    {sinlist[i,2]<-colnames(A)[j]}
                }
            sinlist[,c(1,2)]<-as.character(sinlist[,c(1,2)])
            sinlist[,c(3,4)]<-round(as.numeric(sinlist[,c(3,4)]),3)
        }
    }
    
    
    colnames(repmat)<-colnames(A)
    row.names(repmat)<-colnames(A)
    
    diag(A)<-0
    diag(B)<-0
    possibleA<-sum(ifelse(A!=0,1,0)/2)
    possibleB<-sum(ifelse(B!=0,1,0)/2)
    percentA<-count/possibleA
    percentB<-count/possibleB
    densityA<-possibleA/((ncol(A)^2-ncol(A))/2)
    densityB<-possibleB/((ncol(B)^2-ncol(B))/2)
    
    mat<-matrix(0,nrow=nrow(A),ncol=ncol(A))
    wc<-0
    wveca<-0
    wvecb<-0
    for(i in 1:ncol(A))
        for(j in 1:nrow(A))
            if(A[i,j]&&B[i,j]!=0)
            {
                mat[i,j]<-abs(A[i,j]-B[i,j])
                wc<-wc+1
                wveca[wc]<-A[i,j]
                wvecb[wc]<-B[i,j]
            }
    
    corr<-cor(wveca,wvecb,method=corr)
    
    m<-0
    vec<-0
    for(i in 1:nrow(A))
        for(j in 1:ncol(A))
            if(mat[i,j]!=0)
            {m<-m+1
            vec[m]<-mat[i,j]
            mvec<-mean(vec)
            svec<-sd(vec)}else if(all(mat==0)){mvec<-0
            svec<-0}
    
    res <- list()
    
    #Between networks
    res$replicatedEdges <- as.data.frame(sinlist)
    res$replicated <- count
    res$meanDiff <- mvec
    res$sdDiff <- svec
    res$cor <- corr
    
    #For network A
    res[[nameA]]$totalEdges <- possibleA
    res[[nameA]]$percentage <- percentA
    res[[nameA]]$density <- densityA
    
    #For network B
    res[[nameB]]$totalEdges <- possibleB
    res[[nameB]]$percentage <- percentB
    res[[nameB]]$density <- densityB
    
    return(res)
}
#----