#' Neural Network Filter
#' @description Applies a network filtering methodology to neural network array.
#' Removes edges from the neural network output from \code{\link{convertConnBrainMat}} 
#' using a network filtering approach
#' 
#' @param neuralarray Array from \code{\link{convertConnBrainMat}} function
#' 
#' @param method Filtering method to be applied
#' 
#' @param progBar Should progress bar be displayed?
#' Defaults to \code{TRUE}.
#' Set \code{FALSE} for no progress bar
#' 
#' @param ... Additional arguments from network filtering methods
#' 
#' @return Returns an array of n x n x m filtered matrices
#' 
#' @examples
#' \dontrun{neuralarray <- convertConnBrainMat()
#' 
#' filteredneuralarray <- neuralnetfilter(neuralarray, method = "threshold", thresh = .50)
#' 
#' dependencyarray <- depna(neuralarray)
#' 
#' filtereddependencyarray <- neuralnetfilter(dependencyarray, method = "TMFG", depend = TRUE)
#' }
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Neural Network Filter----
neuralnetfilter <- function (neuralarray, method = c("TMFG","MaST","ECOplusMaST","ECO","threshold"),progBar = TRUE, ...)
{
    if(missing(method))
    {method<-"TMFG"
    }else{method<-match.arg(method)}
    
    n<-length(neuralarray)/nrow(neuralarray)/ncol(neuralarray)  
    
    for(i in 1:n)    
        if(nrow(neuralarray)!=ncol(neuralarray))
        {stop(paste("Participant ",i,"'s matrix is not symmetric",sep=""))}
    
    filarray<-neuralarray
    
    if(progBar)
    {pb <- txtProgressBar(max=n, style = 3)}
    
    for(i in 1:n)
    {
        if(method=="TMFG")
        {filarray[,,i]<-TMFG(neuralarray[,,i],...)$A
        }else if(method=="MaST")
        {filarray[,,i]<-MaST(neuralarray[,,i],...)
        }else if(method=="ECO")
        {filarray[,,i]<-ECO(neuralarray[,,i],...)
        }else if(method=="ECOplusMaST")
        {filarray[,,i]<-ECOplusMaST(neuralarray[,,i],...)
        }else if(method=="threshold")
        {filarray[,,i]<-threshold(neuralarray[,,i],...)$A
        }else{stop("Method not available")}
        if(progBar){setTxtProgressBar(pb, i)}
    }
    if(progBar){close(pb)}
    
    return(filarray)
}
#----