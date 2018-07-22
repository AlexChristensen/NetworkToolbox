#' Generates a Lattice Network
#' @description Generates a lattice network
#' 
#' @param nodes Number of nodes in lattice network
#' 
#' @param edges Number of edges in lattice network
#' 
#' @return Returns an adjacency matrix of a lattice network
#' 
#' @examples
#' latt <- lattnet(10, 27)
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
#Lattice Network----
lattnet <- function (nodes, edges)
{
    dlat<-matrix(0,nrow=nodes,ncol=nodes)
    lat<-matrix(0,nrow=nodes,ncol=nodes)
    
    for(i in 1:nodes)
    {
        if(i!=nodes)
        {dlat[i,i+1]<-1}
        if(i<nodes-1)
        {dlat[i,i+2]<-1}
    }
    lat<-dlat+t(dlat)
    
    over<-sum(lat)-edges
    if(over!=0)
    {rp<-sample(which(dlat==1))
    for(i in 1:over)
    {lat[rp[i]]<-0}}
    
    return(lat)   
}
#----