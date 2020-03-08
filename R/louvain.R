#' Louvain Community Detection Algorithm
#' @description Computes a vector of communities (community) and a global modularity measure (Q)
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param gamma Defaults to \code{1}.
#' Set to \code{gamma} > 1 to detect smaller modules and \code{gamma} < 1 for larger modules
#' 
#' @param M0 Input can be an initial community vector.
#' Defaults to \code{NULL}
#' 
#' @return Returns a list containing:
#' 
#' \item{community}{A community vector corresponding to each node's community}
#' 
#' \item{Q}{Modularity statistic. A measure of how well the communities are compartmentalized}
#' 
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' modularity <- louvain(A)
#' 
#' @references
#' Blondel, V. D., Guillaume, J. L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks. 
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}, P10008.
#'  
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{NeuroImage}, \emph{52}, 1059-1069.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Louvain Community Detection----
louvain <- function (A, gamma, M0)
{
    if(missing(gamma))
    {gamma<-1
    }else{gamma<-gamma}
    
    if(missing(M0))
    {M0<-1:ncol(A)
    }else(M0<-M0)
    
    diag(A) <- 0
    
    n<-ncol(A)
    s<-sum(A)
    
    if(min(A)<0)
    {warning("Matrix contains negative weights: absolute weights were used")
        A<-abs(A)}
    
    Mb<-M0
    M<-M0
    
    mat<-matrix(0,nrow=n,ncol=n)
    for(i in 1:n)
        for(j in 1:n)
        {mat[i,j]<-(colSums(A)[i]*rowSums(A)[j])/s}
    
    B<-A-(gamma*(mat))
    
    B<-(B+t(B))/2
    
    Hnm<-matrix(0,nrow=nrow(A),ncol=(length(unique(Mb))))
    
    
    for(m in 1:max(Mb))
    {
        if(!is.null(nrow(B[,which(Mb==m)])))
        {Hnm[,m]<-rowSums(B[,which(Mb==m)])
        }else{Hnm[,m]<-B[,which(Mb==m)]}
    }
    
    H<-colSums(Hnm)
    Hm<-rowSums(Hnm)
    
    Q0<-(-Inf)
    bsxfun<-matrix(0,nrow=n,ncol=n)
    diag(bsxfun)<-1
    Q<-sum(diag(as.matrix(B)*bsxfun))/s
    
    
    first_iter<-TRUE
    while(Q-Q0>0)
    {
        flag<-TRUE
        while(flag)
        {
            set.seed(0)
            flag<-FALSE
            for(u in sample(n))
            {
                ma<-Mb[u]
                dQ<-Hnm[u,] - Hnm[u,ma] + B[u,u]
                dQ[ma]<-0
                
                max_dQ<-max(dQ)
                mb<-which.max(dQ)
                
                if(max_dQ>0)
                {
                    flag<-TRUE
                    Mb[u]<-mb
                    
                    Hnm[,mb]<-Hnm[,mb]+B[,u]
                    Hnm[,ma]<-Hnm[,ma]-B[,u]
                    Hm[mb]<-Hm[mb]+H[u]
                    Hm[ma]<-Hm[ma]-H[u]
                }
            }
        }
        Mb<-match(Mb,unique(Mb))
        
        M0<-M
        if(first_iter)
        {
            M<-Mb
            first_iter<-FALSE
        }else{
            for(u in 1:n)
            {
                M[M0==u]<-Mb[u]
            }
        }
        
        n<-max(Mb)
        B1<-matrix(0,nrow=n,ncol=n)
        for(u in 1:n)
            for(v in u:n)
            {
                bm<-sum(sum(B[Mb==u,Mb==v]))
                B1[u,v]<-bm
                B1[v,u]<-bm
            }
        B<-B1
        
        Mb<-1:n
        Hnm<-B
        H<-colSums(B)
        
        Q0<-Q
        
        Q<-sum(diag(B))/s
        
    }
    return(list(community=M,Q=Q))
}
#----