#' Dataset Descriptive Statisitcs
#' @description Computes \code{mean}, standard deviaiton (\code{sd}), minimum value (\code{min}),
#' maximum value (\code{max}),
#' and univariate normal statistics (\code{normal?}) for the entire dataset
#' 
#' @param data A matrix or data frame
#' 
#' @return A data frame containing values for \code{n} (number of cases),
#' \code{missing} (number of missing cases), \code{mean}, \code{sd}, \code{min}, and \code{max}. \code{normal?}
#' will contain yes/no for whether the variable is normally distributed based
#' on the \code{\link{shapiro.test}} for the entire dataset
#' 
#' @examples
#' 
#' desc.all(neoOpen)
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Variable Descriptive Statistics----
desc.all <- function(data)
{
    n <- ncol(data)
    
    desc.list <- list()
    
    for(i in 1:n)
    {
        if(all(is.numeric(data[,i])))
        {
                desc.list[[colnames(data)[i]]] <- desc(data,i,histplot=FALSE)
        }
    }
    
    len <- length(desc.list)
    
    desc.mat <- matrix(NA,nrow=len,ncol=7)
    name <- vector("character",length=len)
    
    for(i in 1:len)
    {
        name[i] <- names(desc.list[i])
        desc.mat[i,] <- t(data.frame(unlist(desc.list[[i]]),stringsAsFactors = FALSE))
    }
    
    desc.df <- as.data.frame(desc.mat)
    
    colnames(desc.df) <- c("n","missing","mean","sd","min","max","normal?")
    row.names(desc.df) <- name
    
    return(desc.df)
}
#----