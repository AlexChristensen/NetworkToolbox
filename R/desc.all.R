#' Dataset Descriptive Statisitcs
#' @description Computes mean, standard deviaiton (sd), minimum value (min),
#' maximum value (max),
#' and univariate normal statistics (normal?) for the entire dataset
#' 
#' @param data A matrix or data frame
#' 
#' @return A data frame containing values for n (number of cases),
#' missing (number of missing cases), mean, sd, min, and max. normal?
#' will contain yes/no for whether the variable is normally distributed based
#' on the (\code{\link{shapiro.test}}) for the entire dataset
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