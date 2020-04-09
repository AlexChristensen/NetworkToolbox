#' Dependency Neural Networks
#' @description Applies the dependency network approach to neural network array
#' 
#' @param neuralarray Array from \code{\link{convertConnBrainMat}} function
#' 
#' @param cores Numeric.
#' Number of cores to use in computing results.
#' Set to \code{1} to not use parallel computing.
#' Recommended to use maximum number of cores minus one
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
# Updated 09.04.2020
depna <- function (neuralarray, cores, ...)
{
    # Convert to list for parallelization
    deplist <- list()
    
    # Loop through neuralarray
    for(i in 1:dim(neuralarray)[3])
    {deplist[[i]] <- neuralarray[,,i]}
    
    # Initialize depna array
    deparray <- neuralarray
    
    # Let user know data generation has started
    message("\nComputing dependency matrices...\n", appendLF = FALSE)
    
    # Parallel processing
    cl <- parallel::makeCluster(cores)
    
    # Export variables
    parallel::clusterExport(cl = cl,
                            varlist = c("deplist"),
                            envir=environment())
    
    # Compute dependency matrices
    depnalist <- pbapply::pblapply(X = deplist, cl = cl, FUN = depend,
                                   na.data = "pairwise", progBar = FALSE)
    
    # Stop cluster
    parallel::stopCluster(cl)
    
    # Convert back to array
    for(i in 1:length(depnalist))
    {deparray[,,i] <- depnalist[[i]]}
    
    return(deparray)
}
#----