#' Generates a Random Network
#' @description Generates a random network
#' 
#' @param nodes Number of nodes in random network
#' 
#' @param edges Number of edges in random network
#' 
#' @return Returns an adjacency matrix of a random network
#' 
#' @examples
#' rand <- randnet(10, 27)
#' 
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}, 1059-1069.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Random Network----
randnet <- function (nodes, edges)
{
    mat<-matrix(1,nrow=nodes,ncol=nodes)
    diag(mat)<-0
    ind<-ifelse(upper.tri(mat)==TRUE,1,0)
    i<-which(ind==1)
    rp<-sample(length(i))
    irp<-i[rp]
    
    rand<-matrix(0,nrow=nodes,ncol=nodes)
    rand[irp[1:edges]]<-1
    rand<-rand+t(rand)
    
    return(rand)
}
#----