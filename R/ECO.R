#' ECO Neural Network Filter
#' @description Applies the ECO neural network filtering method
#' 
#' @param data Can be a dataset or a correlation matrix
#' 
#' @param directed Is the network directed?
#' Defaults to FALSE.
#' Set TRUE if the network is directed
#' 
#' @return A sparse association matrix
#' 
#' @examples
#' eco.net <- ECO(neoOpen)
#' 
#' @references 
#' Fallani, F. D. V., Latora, V., & Chavez, M. (2017).
#' A topological criterion for filtering information in complex brain networks.
#' \emph{PLoS Computational Biology}, \emph{13}, e1005305.
#' doi: \href{https://doi.org/10.1371/journal.pcbi.1005305}{10.1371/journal.pcbi.1005305}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#ECO Neural Network Filter----
ECO <- function (data, directed = FALSE)
{
    #corrlation matrix
    if(nrow(data)==ncol(data)){cormat<-data
    }else{cormat<-cor(data)}
    
    C<-cormat
    
    n<-ncol(C)
    S<-C
    if(directed)
    {
        numcon<-3*n
        ind<-which(C!=0)
    }else{C<-upper.tri(C,diag=TRUE)
    numcon<-1.5*n
    ind<-which(upper.tri(C,diag=TRUE)!=0)}
    
    S<-ifelse(C==1,S,0)
    
    if(numcon>length(ind))
    {
        stop("Input matrix is too sparse")
    }
    
    sorind<-matrix(0,nrow=length(ind),ncol=2)
    
    G<-S
    S<-abs(S)
    
    x<-S[ind]
    y<-ind
    h<-cbind(ind,S[ind])
    sorind<-h[order(-h[,2]),]
    C[sorind[(numcon+1):nrow(sorind),1]]<-0
    
    if(directed)
    {W<-C}else{W<-C+t(C)
    diag(W)<-1}
    J<-G+t(G)
    diag(J)<-1
    W<-ifelse(W!=0,J,0)
    W<-as.data.frame(W)
    colnames(W)<-colnames(data)
    row.names(W)<-colnames(data)
    W<-as.matrix(W)
    return(W)
}
#----