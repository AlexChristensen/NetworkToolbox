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
#' \emph{NeuroImage}, \emph{52}, 1059-1069.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Lattice Network----
#Updated 12.05.2021
lattnet <- function (nodes, edges)
{
    dlat<-matrix(0,nrow=nodes,ncol=nodes)
    lat<-matrix(0,nrow=nodes,ncol=nodes)
    
    balance <- sum(lat) - edges
    
    count <- 0
    
    while(sign(balance) == -1){
        
        if(count == 0){
            
            for(i in 1:nodes){
                
                if(i != nodes){
                    dlat[i, (i + 1)] <- 1
                }
            }
            
        }else{
            
            for(i in 1:nodes){
                
                if(i < (nodes - count)){
                    dlat[i, (i + (count + 1))] <- 1
                }
                
            }
            
        }
        
        count <- count + 1
        
        balance <- sum(dlat) - edges
    
    }
    
    over <- sum(dlat) - edges
    
    if(over != 0){
        
        rp <- sample(which(dlat==1), over, replace = FALSE)
        
        dlat[rp] <- 0
        
    }
    
    lat <- dlat + t(dlat)
    
    return(lat)   
}
#----