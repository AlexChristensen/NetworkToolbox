#' Variable Descriptive Statisitcs
#' @description Computes mean, standard deviaiton (sd), minimum value (min),
#' maximum value (max),
#' and univariate normal statistics (normal?) for a variable
#' 
#' @param data A matrix or data frame
#' 
#' @param column Column name or number in \code{data}
#' 
#' @param histplot A histogram plot of the variable
#' 
#' @return A data frame containing values for n (number of cases),
#' missing (number of missing cases), mean, sd, min, and max. normal?
#' will contain yes/no for whether the variable is normally distributed based
#' on the (\code{\link{shapiro.test}})
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
    
    desc.tab <- matrix(NA,nrow=1,ncol=7)
    
    desc.tab[1,1] <- length(vec) - nas
    desc.tab[1,2] <- nas
    desc.tab[1,3] <- round(mean(vec,na.rm=TRUE),3)
    desc.tab[1,4] <- round(sd(vec,na.rm=TRUE),3)
    desc.tab[1,c(5:6)] <- round(range(vec,na.rm=TRUE),3)
    desc.tab[1,7] <- ifelse(sd(vec,na.rm=TRUE)!=0,
                            ifelse(shapiro.test(vec)$p.value<.05,"no","yes"),
                            "uniform")
    
    desc.tab <- as.data.frame(desc.tab)
    
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