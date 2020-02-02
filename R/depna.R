#' Dependency Neural Networks
#' @description Applies the dependency network approach to neural network array
#' 
#' @param neuralarray Array from \code{\link{convertConnBrainMat}} function
#' 
#' @param pB Should progress bar be displayed?
#' Defaults to \code{TRUE}.
#' Set \code{FALSE} for no progress bar
#' 
#' @param ... Additional arguments from \code{\link{depend}} function
#' 
#' @return Returns an array of n x n x m dependency matrices
#' 
#' @examples
#' \dontrun{
#' neuralarray <- convertConnBrainMat()
#' 
#' dependencyneuralarray <- depna(neuralarray)
#' }
#' 
#' @references
#' Jacob, Y., Winetraub, Y., Raz, G., Ben-Simon, E., Okon-Singer, H., Rosenberg-Katz, K., ... & Ben-Jacob, E. (2016).
#' Dependency Network Analysis (DEPNA) reveals context related influence of brain network nodes.
#' \emph{Scientific Reports}, \emph{6}, 27444.
#' 
#' Kenett, D. Y., Tumminello, M., Madi, A., Gur-Gershgoren, G., Mantegna, R. N., & Ben-Jacob, E. (2010).
#' Dominating clasp of the financial sector revealed by partial correlation analysis of the stock market.
#' \emph{PLoS one}, \emph{5}, e15032.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Dependency Network Analysis----
depna <- function (neuralarray, pB = TRUE, ...)
{
    n<-length(neuralarray)/nrow(neuralarray)/ncol(neuralarray)
    
    for(i in 1:n)    
        if(nrow(neuralarray)!=ncol(neuralarray))
        {stop(paste("Participant ",i,"'s matrix is not symmetric",sep=""))}
    
    deparray<-neuralarray
    
    if(pB)
    {pb <- txtProgressBar(max=n, style = 3)}
    
    for(i in 1:n)
    {deparray[,,i]<-depend(neuralarray[,,i],progBar=FALSE,...)
    if(pB){setTxtProgressBar(pb, i)}}
    
    if(pB){close(pb)}
    
    return(deparray)
}
#----