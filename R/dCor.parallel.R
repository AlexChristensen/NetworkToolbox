#' Parallelization of Distance Correlation for ROI Time Series
#' 
#' @description Parallelizes the \code{\link[NetworkToolbox]{dCor}} function
#' for faster computation times
#' 
#' @param neurallist List of lists.
#' A list containing the time series list from all participants imported from the
#' \code{\link[NetworkToolbox]{convertConnBrainMat}} function
#' 
#' @param cores Number of computer processing cores to use when performing covariate analyses.
#' Defaults to \emph{n} - 1 total number of cores.
#' Set to any number between 1 and maximum amount of cores on your computer
#' 
#' @return Returns a \emph{m} x \emph{m} x \emph{n} array corresponding to distance correlations
#' between ROIs (\emph{m} x \emph{m} matrix) for \emph{n} participants
#' 
#' @examples
#' \dontrun{
#' # Import time series data 
#' for(i in 1:5)
# {mat.list[[i]] <- convertConnBrainMat()}
#' 
#' # Run distance correlation
#' dCor.parallel(mat.list, cores = 2)
#' 
#' }
#' 
#' @references
#' Yoo, K., Rosenberg, M. D., Noble, S., Scheinost, D., Constable, R. T., & Chun, M. M. (2019).
#' Multivariate approaches improve the reliability and validity of functional connectivity and prediction of individual behaviors.
#' \emph{NeuroImage}, \emph{197}, 212-223.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Parallelization of distance correlation----
#Updated 03.05.2020
dCor.parallel <- function(neurallist, cores)
{
  ###########################
  #### MISSING ARGUMENTS ####
  ###########################
  
  if(missing(cores))
  {cores <- parallel::detectCores() - 1
  }else{cores <- cores}
  
  # Check data format
  if(!is.matrix(neurallist[[1]][[1]]))
  {
    neurallist <- lapply(neurallist, as.list)
    neurallist <- lapply(neurallist, function(x){lapply(x, as.matrix)})
  }
  
  #######################
  #### MAIN FUNCTION ####
  #######################
  
  # Initialize list
  dCor.list <- list()
  
  # Set up parallelization
  cl <- parallel::makeCluster(cores)
  
  # Export variables
  parallel::clusterExport(cl = cl,
                          varlist = c("dCor.list"),
                          envir = environment())
  
  # Compute distance correlation matrices
  dCor.list <- pbapply::pblapply(X = neurallist, cl = cl, FUN = NetworkToolbox::dCor)
  
  # Stop cluster
  parallel::stopCluster(cl)
  
  # Convert to array
  neuralarray <- simplify2array(dCor.list)
  
  # Change matrix names
  colnames(neuralarray) <- names(neurallist[[1]])
  row.names(neuralarray) <- colnames(neuralarray)
  
  return(neuralarray)
}
#----