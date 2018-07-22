#' Distance
#' @description Computes distance matrix of the network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param weighted Is the network weighted?
#' Defaults to FALSE.
#' Set to TRUE for weighted measure of distance
#' 
#' @return A distance matrix of the network
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' #Unweighted
#' Du <- distance(A)
#' 
#' #Weighted
#' Dw <- distance(A, weighted = TRUE)
#' 
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}, 1059-1069.
#' doi: \href{https://doi.org/10.1016/j.neuroimage.2009.10.003}{10.1016/j.neuroimage.2009.10.003}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Distance----
distance<-function (A, weighted = FALSE)
{
    if(nrow(A)!=ncol(A))
    {stop("Input not an adjacency matrix")}
    if(!weighted)
    {B<-ifelse(A!=0,1,0)
    l<-1
    Lpath<-B
    D<-B
    Idx<-matrix(TRUE,nrow=nrow(B),ncol=ncol(B))
    while(any(Idx))
    {
        l<-l+1
        Lpath<-(as.matrix(Lpath))%*%(as.matrix(B))
        for(e in 1:nrow(Lpath))
            for(w in 1:ncol(Lpath))
                Idx[e,w]<-(Lpath[e,w]!=0&&(D[e,w]==0))
        D[Idx]<-l
    }
    
    D[!D]<-Inf
    diag(D)<-0
    }else if(weighted){
        G<-ifelse(1/A==Inf,0,1/A)
        
        if(any(G==-Inf))
        {G<-ifelse(G==-Inf,0,G)}
        
        if(any(!G==t(G)))
        {if(max(abs(G-t(G)))<1e-10)
        {G<-(G+G)/2}}
        
        n<-ncol(G)
        D<-matrix(Inf,nrow=n,ncol=n)
        diag(D)<-0
        B<-matrix(0,nrow=n,ncol=n)
        
        for(u in 1:n)
        {
            S<-matrix(TRUE,nrow=n,ncol=1)
            L1<-G
            V<-u
            while(TRUE)
            {
                S[V]<-0
                L1[,V]<-0
                for(v in V)
                {
                    W<-which(L1[v,]!=0)
                    d<-apply(rbind(D[u,W],(D[u,v]+L1[v,W])),2,min)    
                    wi<-apply(rbind(D[u,W],(D[u,v]+L1[v,W])),2,which.min)
                    D[u,W]<-d
                    ind<-W[wi==2]
                    B[u,ind]<-B[u,v]+1
                }
                
                minD<-suppressWarnings(min(D[u,S==TRUE]))
                if(length(minD)==0||is.infinite(minD)){break}
                
                V<-which(D[u,]==minD)
            }
        }
    }
    
    D<-ifelse(D==Inf,0,D)
    
    colnames(D)<-colnames(A)
    row.names(D)<-colnames(A)
    return(D)
}
#----