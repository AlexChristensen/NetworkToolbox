#' Randomized Shortest Paths Betweenness Centrality
#' @description Computes betweenness centrality based on randomized shortest paths
#' of each node in a network
#' (\strong{Please see and cite Kivimaki et al., 2016})
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param beta Sets the beta parameter.
#' Defaults to \code{0.01} (recommended).
#' Beta > 0.01 measure gets closer to weighted
#' betweenness centrality (10) and beta < 0.01
#' measure gets closer to degree (.0001)
#' 
#' @param comm Vector.
#' Community vector containing a value for each node.
#' Computes "bridge" RSPBC, where the number of times
#' a node is used on a random path between to another community
#' 
#' @return A vector of randomized shortest paths betweenness
#' centrality values for each node in the network
#' 
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#' 
#' rspbc <- rspbc(A, beta=0.01)
#' 
#' @references 
#' Kivimaki, I., Lebichot, B., Saramaki, J., & Saerens, M. (2016).
#' Two betweenness centrality measures based on Randomized Shortest Paths.
#' \emph{Scientific Reports}, \emph{6}, 19668.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Randomized Shortest Paths Betweennesss Centrality----
rspbc <- function (A, beta = 0.01, comm = NULL)
{
    if(nrow(A)!=ncol(A))
    {stop("Input not an adjacency matrix")}
    
    if(is.null(comm))
    {A <- abs(A)
    }else{A[A<0] <- 0}
    
    A <- as.matrix(A)
    
    n<-ncol(A)
    e<-matrix(1,nrow=n,ncol=1)
    I<-diag(1,nrow=nrow(A),ncol=ncol(A))
    degs<-as.matrix(A)%*%as.matrix(e)
    
    if(any(degs==0))
    {stop("Graph contains unconnected nodes")}
    
    D1<-matrix(0,nrow=nrow(I),ncol=ncol(I))
    for(i in 1:nrow(I))
        for(j in 1:ncol(I))
            if(I[i,j]==1)
            {D1[i,j]<-I[i,j]/degs[i]}
    
    Pref<-as.matrix(D1)%*%as.matrix(A)
    
    bets<-matrix(0,nrow=n,ncol=1)
    C<-1/A
    C<-as.matrix(C)
    C[is.infinite(C)]<-0
    W<-Pref*exp(-(beta)*C)
    
    #Bridge RSPBC
    if(!is.null(comm))
    {
        for(i in 1:length(unique(comm)))
        {
            combos <- t(combn(which(comm==unique(comm)[i]),2))
            
            for(j in 1:nrow(combos))
            {W[c(combos[j,]),c(combos[j,])] <- 0}
        }
    }
    
    Y<-I-W
    Z<-solve(Y,I)
    Zdiv<-1/Z
    Zdiv[Zdiv==Inf]<-0
    DZdiv<-matrix(0,nrow=nrow(Zdiv),ncol=ncol(Zdiv))
    diag(DZdiv)<-diag(Zdiv)
    
    bet<-diag(as.matrix(Z)%*%as.matrix(t(Zdiv-n*DZdiv))%*%as.matrix(Z))
    bet<-round(as.data.frame(bet),0)
    minimum <- min(bet) - 1
    bet <- bet - minimum
    bet<-as.vector(as.matrix(bet))
    names(bet)<-colnames(A)
    return(bet)
}
#----