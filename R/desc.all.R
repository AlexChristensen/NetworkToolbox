#' Dataset Descriptive Statistics
#' @description Computes \code{mean}, standard deviation (\code{sd}), minimum value (\code{min}),
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
    # Number of variables
    n <- ncol(data)
    
    # Descriptives list
    desc.list <- list()
    
    # Loop for descriptives
    for(i in 1:n)
    {
        # Check for variables that are factors
        # that can be converted to numeric
        if(is.factor(data[,i]))
        {data[,i] <- suppressWarnings(as.numeric(levels(data[,i]))[data[,i]])}
        
        # If variables are numeric,
        # then get descriptives
        if(all(is.numeric(data[,i])))
        {desc.list[[colnames(data)[i]]] <- desc(data,i,histplot=FALSE)}
    }
    
    # Get number of descriptive variables
    len <- length(desc.list)
    
    # Initialize descriptive matrix
    desc.mat <- matrix(NA,nrow=len,ncol=7)
    # Initialize variable name vector
    name <- vector("character",length=len)
    
    # Loop for descriptive matrix
    for(i in 1:len)
    {
        # Grab names
        name[i] <- names(desc.list[i])
        
        # Grab descriptives
        desc.mat[i,] <- t(data.frame(unlist(desc.list[[i]]),stringsAsFactors = FALSE))
    }
    
    # Convert descriptive matrix to a data frame
    desc.df <- as.data.frame(desc.mat)
    
    # Change column names to descriptives
    colnames(desc.df) <- c("n","missing","mean","sd","min","max","normal?")
    # Change row names to variable names
    row.names(desc.df) <- name
    
    return(desc.df)
}
#----