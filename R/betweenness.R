#' Betwenness Centrality
#' @description Computes betweenness centrlaity of each node in a network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param weighted Is the network weighted?
#' Defaults to TRUE.
#' Set to FALSE for unweighted measure of betwenness centrality
#' 
#' @return A vector of betweenness centrality values for each node in the network
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' #Weighted BC
#' BCw <- betweenness(A)
#' 
#' #Unweighted BC
#' BC <- betweenness(A, weighted = FALSE)
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
#Betweenness Centrality----
betweenness <- function (A, weighted = TRUE)
{
    if(nrow(A)!=ncol(A))
    {stop("Input not an adjacency matrix")}
    
    A <- abs(A)
    
    if(!weighted)
    {
    B <- binarize(A)
    n<-ncol(B)
    I<-diag(60)
    d<-1
    NPd<-B
    NSPd<-NPd
    NSP<-NSPd
    diag(NSP)<-1
    L<-NSPd
    diag(L)<-1
    while (!is.na(which(NSPd!=0)[1]))
    {
        d<-d+1
        NPd<-as.matrix(NPd)%*%as.matrix(B)
        NSPd<-NPd*(L==0)
        NSP<-NSP+NSPd
        L<-L+d*(NSPd!=0)
    }
    L[!L]<-Inf
    diag(L)<-0
    NSP[!NSP]<-1
    Bt<-t(B)
    DP<-matrix(0,nrow=nrow(B),ncol=ncol(B))
    diam<-d-1
    
    for(d in diam:2)
    {
        DPd1<- (as.matrix(((L==d)*(1+DP)/NSP))%*%as.matrix(Bt))*((L==(d-1))*NSP)
        DP<-DP+DPd1
    }
    BC<-round(as.matrix(colSums(DP),ncol=ncol(A)),0)
    }else{
        G<-ifelse(1/A==Inf,0,1/A)
        
        if(any(G==-Inf))
        {G<-ifelse(G==-Inf,0,G)}
        
        if(any(!G==t(G)))
        {
            if(max(abs(G-t(G)))<1e-10){G<-(G+G)/2}
        }
        
        n<-ncol(G)
        
        BC<-matrix(0,nrow=n,ncol=1)
        
        for(u in 1:n)
        {
            D<-matrix(Inf,nrow=n,ncol=1)
            D[u]<-0
            NP<-matrix(0,nrow=n,ncol=1)
            NP[u]<-1
            S<-matrix(TRUE,nrow=n,ncol=1)
            P<-matrix(FALSE,nrow=n,ncol=n)
            Q<-matrix(0,nrow=n,ncol=1)
            q<-n
            
            G1<-G
            V<-u
            
            while(TRUE)
            {
                S[V]<-0
                G1[,V]<-0
                for(v in V)
                {
                    Q[q]<-v
                    q<-q-1
                    W<-which(G1[v,]!=0)
                    
                    for(w in W)
                    {
                        Duw<-D[v]+G1[v,w]
                        if(Duw<D[w])
                        {
                            D[w]<-Duw
                            NP[w]<-NP[v]
                            P[w,]<-0
                            P[w,v]<-1
                        }else if(Duw==D[w])
                        {
                            NP[w]<-NP[w]+NP[v]
                            P[w,v]<-1
                        }
                    }
                }
                
                
                minD<-suppressWarnings(min(D[S==TRUE]))
                if(length(minD)==0){break}else if(is.infinite(minD))
                {Q[1:q]<-ifelse(length(which(is.infinite(D)))==0,break,which(is.infinite(D)))
                break}
                V<-which(D==minD)
            }
            
            DP<-matrix(0,nrow=n,ncol=1)
            for(w in Q[1:n-1])
            {BC[w]<-BC[w]+DP[w]
            for(v in which(P[w,]!=0))
                DP[v]<-(DP[v]+(1+DP[w]))*NP[v]/NP[w]}
            
        }
        BC<-round(as.matrix(BC,ncol=ncol(A)),0)
    }
    BC<-as.vector(BC)
    
    names(BC)<-colnames(A)
    
    return(BC)
}
#----