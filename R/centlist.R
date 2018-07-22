#' List of Centrality Measures
#' @description Computes centrality measures of the network
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param weighted Is the network weighted?
#' Defaults to TRUE.
#' Set to FALSE for unweighted list of centrality measures
#' 
#' @return Returns a list containing:
#' 
#' \item{BC}{Betweenness centrality}
#' 
#' \item{LC}{Closeness centrality}
#' 
#' \item{k}{Degree (weighted = FALSE)}
#' 
#' \item{Str}{Strength (weighted = TRUE)}
#' 
#' \item{EC}{Eigenvector centrality}
#' 
#' \item{lev}{Leverage centrality}
#' 
#' \item{HC}{Hybrid centrality (weighted = TRUE)}
#' 
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' #Weighted
#' centW <- centlist(A, weighted = TRUE)
#' 
#' #Unweighted
#' centU <- centlist(A, weighted = FALSE)
#' 
#' @references 
#' Pozzi, F., Di Matteo, T., & Aste, T. (2013).
#' Spread of risk across financial markets: Better to invest in the peripheries. 
#' \emph{Scientific Reports}, \emph{3}(1655), 1-7.
#' doi: \href{https://doi.org/10.1038/srep01665}{10.1038/srep01665}
#' 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}, 1059-1069.
#' doi: \href{https://doi.org/10.1016/j.neuroimage.2009.10.003}{10.1016/j.neuroimage.2009.10.003}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Centrality List----
centlist <- function (A, weighted = TRUE)
{
    if(nrow(A)!=ncol(A))
    {stop("Input not an adjacency matrix")}
    if(!weighted)
    {
        BC<-betweenness(A,weighted=FALSE)
        CC<-closeness(A,weighted=FALSE)
        Deg<-degree(A)
        EC<-eigenvector(A,weighted=FALSE)
        lev<-leverage(A,weighted=FALSE)
    
        return(list(BC=BC,LC=CC,k=Deg,EC=EC,lev=lev))
    }else{
        BC<-betweenness(A)
        CC<-closeness(A)
        Str<-strength(A)
        EC<-eigenvector(A)
        lev<-leverage(A)
        hyb<-hybrid(A)
        
        return(list(BC=BC,LC=CC,Str=Str,EC=EC,lev=lev,HC=hyb))}
}
#----