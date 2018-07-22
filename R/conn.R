#' Network Connectivity
#' @description Computes the average and standard deviation of the weights in the network
#' 
#' @param A An adjacency matrix of network A
#' 
#' @return Returns a list containing:
#' 
#' \item{weights}{Each edge weight in the network}
#' 
#' \item{mean}{The mean of the edge weights in the network}
#' 
#' \item{sd}{The standard deviation of the edge weights in the network}
#' 
#' \item{total}{The sum total of the edge weights in the network}
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' connectivity <- conn(A)
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Network Connectivity----
conn <- function (A)
{
    diag(A)<-0
    
    weights<-0
    wc<-0
    B<-A[lower.tri(A)]
    for(i in 1:length(B))
        if (B[i]!=0)
        {
            wc <- wc+1
            weights[wc] <- B[i]
        }
    tot<-sum(weights)
    mea<-mean(weights)
    s<-sd(weights)
    
    possible<-sum(ifelse(A!=0,1,0)/2)
    den<-possible/((ncol(A)^2-ncol(A))/2)
    
    return(list(weights=weights,mean=mea,sd=s,total=tot,density=den))
}
#----