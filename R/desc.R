#' Variable Descriptive Statistics
#' @description Computes \code{mean}, standard deviation (\code{sd}), minimum value (\code{min}),
#' maximum value (\code{max}),
#' and univariate normal statistics (\code{normal?}) for a variable
#' 
#' @param data A matrix or data frame
#' 
#' @param column Column name or number in \code{data}
#' 
#' @param histplot A histogram plot of the variable
#' 
#' @return A data frame containing values for \code{n} (number of cases),
#' \code{missing} (number of missing cases), \code{mean}, \code{sd}, \code{min}, and \code{max}. \code{normal?}
#' will contain yes/no for whether the variable is normally distributed based
#' on the \code{\link{shapiro.test}} for a variable
#' 
#' @examples
#' 
#' desc(neoOpen,1)
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @importFrom stats shapiro.test
#' 
#' @export
#Variable Descriptive Statistics----
desc <- function(data, column, histplot = TRUE)
{
    if(missing(column))
    {stop("Column name or number must be input")}
    
    if(is.character(column))
    {num <- which(colnames(data)==column)
    }else{num <- column}
    
    vec <- data[,num]
    
    len <- length(which(is.na(vec)))
    
    nas <- ifelse(len==0,0,len)
    
    # shapiro.test requires between 3 and 5000 (non-missing) values
    n.valid <- length(vec) - nas
    
    if(sd(vec,na.rm=TRUE)==0)
    {normal <- "uniform"
    }else if(n.valid >= 3 && n.valid <= 5000)
    {normal <- ifelse(shapiro.test(vec[!is.na(vec)])$p.value<.05,"no","yes")
    }else{normal <- NA_character_}
    
    desc.tab <- data.frame(n = n.valid,
                           missing = nas,
                           mean = round(mean(vec,na.rm=TRUE),3),
                           sd = round(sd(vec,na.rm=TRUE),3),
                           min = round(min(vec,na.rm=TRUE),3),
                           max = round(max(vec,na.rm=TRUE),3),
                           normal = normal,
                           stringsAsFactors = FALSE)
    
    colnames(desc.tab) <- c("n","missing","mean","sd","min","max","normal?")
    row.names(desc.tab) <- column
    
    if(histplot)
    {
        hist(vec, main = paste("Histogram of ",column,sep=""),
             xlab = paste(column),breaks=20)
    }
    
    return(desc.tab)
    
}
#----