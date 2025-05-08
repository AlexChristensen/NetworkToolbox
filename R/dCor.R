#' Distance Correlation for ROI Time Series
#'
#' @description Computes the distance correlation (Yoo et al., 2019) for
#' ROI time series data. This function is mainly a subroutine for the
#' \code{\link[NetworkToolbox]{dCor.parallel}} function
#'
#' @param neurallist List.
#' A time series list from \code{\link[NetworkToolbox]{convertConnBrainMat}} function
#'
#' @param centering Character.
#' Options for centering the Euclidean distances.
#'
#' \itemize{
#' \item \code{"U"} --- Uses number of time points minus 2 in the computation of the mean
#'
#' \item \code{"double"} --- Uses the mean
#'
#' }
#'
#' @return Returns a \emph{m} x \emph{m} matrix corresponding to distance correlations
#' between ROIs
#'
#' @examples
#' \dontrun{
#' # Import time series data
#' neurallist <- convertConnBrainMat()
#'
#' # Run distance correlation
#' dCor(neurallist)
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
#' @importFrom stats dist
#'
#' @export
#Distance correlation----
#Updated 07.03.2020
dCor <- function(neurallist, centering = c("U", "double"))
{
  ###########################
  #### MISSING ARGUMENTS ####
  ###########################

  if(missing(centering))
  {centering <- "U"
  }else{centering <- match.arg(centering)}

  # Number of ROIs
  ROIs <- length(neurallist)

  # Length of time-series
  nTime <- nrow(neurallist[[1]])

  # Initialize array
  dist.arr <- array(0, dim = c(nTime, nTime, ROIs))

  ######################
  #### COMPUTE dCOR ####
  ######################

  for(i in 1:ROIs)
  {
    # Check for ROIs that are zeros
    if(any(apply(neurallist[[i]],2,sum) == 0))
    {neurallist[[i]] <- neurallist[[i]][,-which(colSums(neurallist[[i]]) == 0)]}

    # Compute Euclidean distance
    euclid <- dist(neurallist[[i]])

    # 1D distance array to 2D matrix
    n2 <- ceiling(sqrt(length(euclid)*2))
    eDist <- matrix(0, nrow = n2, ncol = n2)
    eDist[lower.tri(eDist)] <- euclid
    eDist <- eDist + t(eDist)

    # Centering
    if(centering == "U")
    {
      eDistCent <- eDist - matrix(colSums(eDist) / (nTime - 2),
                                  nrow = nrow(eDist),
                                  ncol = nTime,
                                  byrow = FALSE) - matrix(colSums(eDist) / (nTime - 2),
                                                          nrow = nrow(eDist),
                                                          ncol = nTime,
                                                          byrow = TRUE) + matrix((sum(colSums(eDist))) / ((nTime-1) * (nTime-2)),
                                                                                 nrow = nTime,
                                                                                 ncol = nTime)
      diag(eDistCent) <- 0
    }else if(centering == "double")
    {
      eDistCent <- eDist - matrix(colMeans(eDist),
                                  nrow = nrow(eDist),
                                  ncol = nTime,
                                  byrow = FALSE) - matrix(colMeans(eDist),
                                                          nrow = nrow(eDist),
                                                          ncol = nTime,
                                                          byrow = TRUE) + matrix((mean(colMeans(eDist))),
                                                                                 nrow = nTime,
                                                                                 ncol = nTime)
    }

    # Input into array
    dist.arr[,,i] <- eDistCent
  }

  # Compute distance variance
  K <- switch(centering,
              U = nTime * (nTime - 3),
              double = nTime^2)

  dVar <- colSums(colSums(dist.arr^2)) / K

  # Compute distnace covariance
  dCov <- matrix(0, nrow = ROIs, ncol = ROIs)

  for(i in 1:(ROIs-1))
    for(j in (i+1):ROIs)
    {dCov[i,j] <- sum(colSums(dist.arr[,,i] * dist.arr[,,j])) / K}

  dCov <- dCov + t(dCov)
  diag(dCov) <- dVar

  # Compute distance correlation
  dVarSqrt <- sqrt(dVar %*% t(dVar))
  divCor <- dCov / dVarSqrt
  signs <- sign(divCor)
  dCor <- sqrt(abs(divCor))
  dCor <- signs * dCor
  dCor[dCov <= 0] <- 0
  diag(dCor) <- 0

  return(dCor)

}
#----