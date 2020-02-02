#' Dependency Network Approach
#' @description Generates a dependency matrix of the data (index argument is still in testing phase)
#' 
#' @param data A set of data
#' 
#' @param normal Should data be transformed to a normal distribution?
#' Defaults to \code{FALSE}. Data is not transformed to be normal.
#' Set to \code{TRUE} if data should be transformed to be normal
#' (computes correlations using the \link[qgraph]{cor_auto} function)
#' 
#' @param na.data How should missing data be handled?
#' For \code{"listwise"} deletion the \code{\link{na.omit}} function is applied.
#' Set to \code{"fiml"} for Full Information Maximum Likelihood (\link[psych]{corFiml}).
#' Full Information Maximum Likelihood is \strong{recommended} but time consuming
#' 
#' @param index Should correlation with the latent variable
#' (i.e., weighted average of all variables) be removed?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} to remove common latent factor
#' 
#' @param fisher Should Fisher's Z-test be used to keep significantly higher influences (index only)?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} to remove non-significant influences
#' 
#' @param progBar Should progress bar be displayed?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} for no progress bar
#' 
#' @return Returns an adjacency matrix of dependencies
#' 
#' @examples
#' \dontrun{
#' D <- depend(neoOpen)
#' 
#' Dindex <- depend(neoOpen, index = TRUE)
#' }
#' @references
#' Kenett, D. Y., Tumminello, M., Madi, A., Gur-Gershgoren, G., Mantegna, R. N., & Ben-Jacob, E. (2010).
#' Dominating clasp of the financial sector revealed by partial correlation analysis of the stock market.
#' \emph{PLoS one}, \emph{5}, e15032.
#' 
#' Kenett, D. Y., Huang, X., Vodenska, I., Havlin, S., & Stanley, H. E. (2015).
#' Partial correlation analysis: Applications for financial markets.
#' \emph{Quantitative Finance}, \emph{15}, 569-578.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @importFrom utils txtProgressBar setTxtProgressBar
#' 
#' @export
#Dependency----
depend <- function (data, normal = FALSE,
                    na.data = c("pairwise","listwise","fiml", "none"),
                    index = FALSE, fisher = FALSE, progBar = TRUE)
{
    #missing data handling
    if(missing(na.data))
    {
        if(any(is.na(data)))
        {stop("Missing values were detected! Set 'na.data' argument")
        }else{na.data<-"none"}
    }else{na.data<-match.arg(na.data)}
    
    if(na.data=="pairwise")
    {
        if(normal)
        {cormat<-qgraph::cor_auto(data,missing=na.data)
        }else{cormat<-psych::cor2(data,use=na.data)}
    }else if(na.data=="listwise")
    {
        if(normal)
        {cormat<-qgraph::cor_auto(data,missing=na.data)
        }else{
            rem<-na.action(na.omit(data))
            warning(paste(length(na.action(na.omit(data)))),
                    " rows were removed for missing data\nrow(s): ",
                    paste(na.action(na.omit(data)),collapse = ", "))
            data<-na.omit(data)
        }
    }else if(na.data=="fiml")
    {
        if(normal)
        {cormat<-qgraph::cor_auto(data,missing=na.data)
        }else{cormat<-psych::corFiml(data)}
    }else if(na.data=="none")
    {
        if(nrow(data)==ncol(data)){cormat<-data
        }else if(normal){cormat<-qgraph::cor_auto(data)
        }else{cormat<-cor(data)}
    }
    
    inter<-((ncol(cormat)*(ncol(cormat)-1)*(ncol(cormat)-2)))
    
    if(index)
    {
        m<-rowMeans(data)
        dat<-cbind(data,m)
        
        #missing data handling
        
        if(na.data=="pairwise")
        {
            if(normal)
            {cordat<-qgraph::cor_auto(dat,missing=na.data)
            }else{cordat<-cor(dat,use="pairwise.complete.obs")}
        }else if(na.data=="listwise")
        {
            if(normal)
            {cordat<-qgraph::cor_auto(dat,missing=na.data)
            }else{
                rem<-na.action(na.omit(dat))
                warning(paste(length(na.action(na.omit(dat)))),
                        " rows were removed for missing data\nrow(s): ",
                        paste(na.action(na.omit(dat)),collapse = ", "))
                dat<-na.omit(dat)
            }
        }else if(na.data=="fiml")
        {
            if(normal)
            {cordat<-qgraph::cor_auto(dat,missing=na.data)
            }else{cormat<-psych::corFiml(dat)}
        }else if(na.data=="none")
        {
            if(nrow(dat)==ncol(dat)){cordat<-dat
            }else if(normal){cordat<-qgraph::cor_auto(dat)
            }else{cordat<-cor(dat)}
        }
        
        indpartial <- function (data,i,k,m=ncol(cordat))
        {(data[i,k]-(data[i,m]*data[k,m]))/sqrt((1-(data[i,m]^2))*(1-(data[k,m]^2)))}
        indmat<-matrix(0,nrow=nrow(cordat)-1,ncol=ncol(cordat)-1)
        for(i in 1:ncol(cordat)-1)
            for(k in 1:ncol(cordat)-1)
                if(i!=k)
                {indmat[i,k]<-cordat[i,k]-indpartial(cordat,i,k)}
        
        if(progBar){pb <- txtProgressBar(max=inter, style = 3)}
        count<-0
        
        partial <- function (data,i,k,j)
        {(data[i,k]-(data[i,j]*data[k,j]))/sqrt((1-(data[i,j]^2))*(1-(data[k,j]^2)))}
        
        z <- function (r)
        {.5*log((1+r)/(1-r))}
        
        parmat<-array(0,dim=c(nrow=ncol(indmat),ncol=ncol(indmat),ncol(indmat)))
        for(i in 1:ncol(indmat))
            for(k in 1:ncol(indmat))
                for(j in 1:ncol(indmat))
                    if(i!=j&&k!=j&&i!=k)
                    {count<-count+1
                    parmat[i,k,j]<-(z(indmat[i,k])-z(partial(indmat,i,k,j)))
                    if(progBar){setTxtProgressBar(pb, count)}}
        if(progBar){close(pb)}
    }
    
    if(!index)
    {
        if(progBar){pb <- txtProgressBar(max=inter, style = 3)}
        count<-0
        
        partial <- function (data,i,k,j)
        {(data[i,k]-(data[i,j]*data[k,j]))/sqrt((1-(data[i,j]^2))*(1-(data[k,j]^2)))}
        
        parmat<-array(0,dim=c(nrow=ncol(cormat),ncol=ncol(cormat),ncol(cormat)))
        for(i in 1:ncol(cormat))
            for(k in 1:ncol(cormat))
                for(j in 1:ncol(cormat))
                    if(i!=j&&k!=j&&i!=k)
                    {count<-count+1
                    parmat[i,k,j]<-(cormat[i,k]-partial(cormat,i,k,j))
                    if(progBar){setTxtProgressBar(pb, count)}}
        if(progBar){close(pb)}
    }
    
    for(h in 1:j)
        diag(parmat[,,h])<-1
    
    depmat<-matrix(0,nrow=nrow(parmat),ncol=ncol(parmat))
    for(i in 1:ncol(parmat))
        for(j in 1:ncol(parmat))
        {depmat[j,i]<-mean(parmat[i,-j,j])}
    
    if(fisher)
    {
        fish <- function (r)
        {z<-.5*log((1+abs(r))/(1-abs(r)))
        return(z)}
        
        zsd <- function (n)
        {zsd<-1/(n-3)
        return(zsd)}
        
        fishtest <- function (r1,r2,n1,n2)
        {test<-(fish(r1)-fish(r2))/sqrt(zsd(n1)+zsd(n2))
        return(test)}
        
        sig<-matrix(0,nrow=nrow(indmat),ncol=ncol(indmat))
        for(i in 1:nrow(indmat))
            for(j in 1:ncol(indmat))
                if(i!=j)
                {sig[i,j]<-fishtest(indmat[i,j],depmat[i,j],nrow(data),nrow(data))}
        sig<-ifelse(sig>=1.96,1,0)
    }
    
    colnames(depmat)<-colnames(data)
    
    if(!fisher)
    {return(depmat)
    }else if(fisher)
    {return(list(depmat=depmat,sigmat=sig))}
}
#----