#' Randomized Shortest Paths Betweenness Centrality
#' @description Computes betweenness centrlaity based on randomized shortest paths
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
#' @return A vector of randomized shortest paths betweenness
#' centrality values for each node in the network
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' rspbc <- rspbc(A, beta=0.01)
#' 
#' @references 
#' Kivimaki, I., Lebichot, B., Saramaki, J., & Saerens, M. (2016).
#' Two betweenness centrality measures based on Randomized Shortest Paths.
#' \emph{Scientific Reports}, \emph{6}, 19668.
#' doi: \href{https://doi.org/10.1038/srep19668}{10.1038/srep19668}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Randomized Shortest Paths Betweennesss Centrality----
rspbc <- function (A, beta = 0.01)
{
    if(nrow(A)!=ncol(A))
    {stop("Input not an adjacency matrix")}
    
    A <- abs(A)
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
    rsums<-rowSums(W)
    
    
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
    rownames(bet)<-colnames(A)
    bet<-as.matrix(bet)
    return(bet)
}
#----