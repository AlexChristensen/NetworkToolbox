#' ECO+MaST Network Filter
#' @description Applies the \link[NetworkToolbox]{ECO} neural network filtering method
#' combined with the \link[NetworkToolbox]{MaST} filtering method
#' 
#' @param data Can be a dataset or a correlation matrix
#' 
#' @return A sparse association matrix
#' 
#' @examples
#' ECOplusMaST.net <- ECOplusMaST(neoOpen)
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
#ECO Filter + MaST----
ECOplusMaST <- function (data)
{
    #corrlation matrix
    if(nrow(data)==ncol(data)){cormat<-data
    }else{cormat<-cor(data)}
    
        a<-MaST(data)
        b<-ECO(data)
        k<-matrix(NA,nrow=nrow(a),ncol=ncol(a))
        for(i in 1:nrow(a))
            for(j in 1:ncol(a))
                if(a[i,j]==b[i,j]){k[i,j]<-a[i,j]}else k[i,j]<-a[i,j]+b[i,j]
        
    k<-as.data.frame(k)
    colnames(k)<-colnames(data)
    row.names(k)<-colnames(data)
    k<-as.matrix(k)
    return(k)
}
#----