#' Regression Matrix
#' @description Computes regression such that one variable is regressed over all other variables
#' 
#' @param data A dataset
#' 
#' @param family Error distribution to be used in the regression model.
#' Defaults to "logistic".
#' Set to any family used in function \link[stats]{family}
#' 
#' @param symmetric Should matrix be symmetric?
#' Defaults to TRUE, taking the mean of the two edge weights
#' (i.e., [\emph{i},\emph{j}] and [\emph{j},\emph{i}])
#' Set to FALSE for asymmwtric weights
#' (i.e., [\emph{i},\emph{j}] does not equal [\emph{j},\emph{i}])
#' 
#' @return A matrix of fully regressed coefficients where one variable is regressed over all others
#' 
#' @examples
#' #binarize responses
#' psyb <- ifelse(neoOpen>=4, 1, 0)
#' 
#' #perform logistic regression
#' mat <- reg(psyb)
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @importFrom stats glm
#' 
#' @export 
#Regression----
reg <- function (data,
                 family = c("binomial" ,"gaussian", "Gamma", "poisson"),
                 symmetric = TRUE)
{
    if(missing(family))
    {family<-"binomial"
    }else(family<-match.arg(family))
    
    n<-ncol(data)
    
    data<-as.data.frame(data)
    
    mat<-matrix(0,nrow=(n-1),ncol=n)
    
    pb <- txtProgressBar(max=n, style = 3)
    for(i in 1:ncol(data))
    {
        res<-cbind(data[,i],data[,-i])
        mat[,i]<-glm(res,family=family)$coefficients[2:(ncol(data))]
        
        setTxtProgressBar(pb, i)
    }
    close(pb)
    
    nmat<-matrix(0,nrow=n,ncol=n)
    
    for(i in 1:n)
    {
        if(i==1)
        {nmat[,i]<-c(0,mat[,i])
        }else if(i!=n)
        {nmat[,i]<-c(mat[1:(i-1),i],0,mat[i:nrow(mat),i])
        }else if(i==n)
        {nmat[,i]<-c(mat[,i],0)}
    }
    
    if(symmetric)
    {nmat<-(nmat+t(nmat))/2}
    
    row.names(nmat)<-colnames(data)
    colnames(nmat)<-colnames(data)
    
    return(nmat)
}
#----