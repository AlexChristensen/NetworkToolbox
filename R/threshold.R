#' Threshold Network Estimation Methods
#' @description Filters the network based on an r-value, alpha, adaptive alpha,
#' bonferroni, false-discovery rate (FDR), or proportional density (fixed number of edges) value
#' 
#' @param data Can be a dataset or a correlation matrix
#' 
#' @param normal Should data be transformed to a normal distribution?
#' Defaults to FALSE.
#' Data is not transformed to be normal.
#' Set to TRUE if data should be transformed to be normal
#' (computes correlations using the \link[qgraph]{cor_auto} function)
#' 
#' @param a When \code{thresh = "alpha"}, \code{"adaptive"}, and \code{"bonferroni"}
#' an \eqn{\alpha} threshold is applied (defaults to \code{.05}).
#' For \code{"adaptive"}, beta (Type II error) is set to \eqn{\alpha*5} for a medium effect size (\emph{r} = .3).
#' When \code{thresh = "FDR"}, a q-value threshold is applied (defaults to \code{.10}).
#' When \code{thresh = "proportional"}, a density threshold is applied (defaults to \code{.15})
#' 
#' @param thresh Sets threshold. Defaults to \code{"alpha"}.
#' Set to any value 0> \emph{r} >1 to retain values greater than set value,
#' \code{"adaptive"} for an \code{\link{adapt.a}} based on sample size (Perez & Pericchi, 2014),
#' \code{"bonferroni"} for the bonferroni correction on alpha,
#' \code{"FDR"} for local false discovery rate,
#' and \code{"proportional"} for a fixed density of edges (keeps strongest correlations within density)
#' 
#' @param na.data How should missing data be handled?
#' For \code{"listwise"} deletion the \code{\link{na.omit}} function is applied.
#' Set to \code{"fiml"} for Full Information Maxmimum Likelihood (\code{\link[psych]{corFiml}}).
#' Full Information Maxmimum Likelihood is \strong{recommended} but time consuming
#' 
#' @param ... Additional arguments for \code{\link[fdrtool]{fdrtool}} and \code{\link[NetworkToolbox]{adapt.a}}
#' 
#' @return Returns a list containing:
#' 
#' \item{A}{The filtered adjacency matrix}
#' 
#' \item{r.cv}{The critical correlation value used to filter the network}
#' 
#' @examples
#' threshnet<-threshold(neoOpen)
#' 
#' alphanet<-threshold(neoOpen, thresh = "alpha", a = .05)
#' 
#' bonnet<-threshold(neoOpen, thresh = "bonferroni", a = .05)
#' 
#' FDRnet<-threshold(neoOpen, thresh = "FDR", a = .10)
#' 
#' propnet<-threshold(neoOpen, thresh = "proportional", a = .15)
#' @references 
#' Strimmer, K. (2008).
#' fdrtool: A versatile R package for estimating local and tail area-based false discovery rates.
#' \emph{Bioinformatics}, \emph{24}, 1461-1462.
#' doi: \href{https://doi.org/10.1093/bioinformatics/btn209}{10.1093/bioinformatics/btn209}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Threshold filtering----
threshold <- function (data, a,
                       thresh = c("alpha","adaptive","bonferroni","FDR","proportional"),
                       normal = FALSE,
                       na.data = c("pairwise","listwise","fiml","none"), ...)
{
    if(missing(thresh))
    {thresh<-"alpha"
    }else{thresh<-thresh}
    
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
        }else{cormat<-cor(data,use="pairwise.complete.obs")}
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
    }
    
    #corrlation matrix
    if(nrow(data)==ncol(data)){cormat<-data
    }else if(normal){cormat<-qgraph::cor_auto(data)
    }else{cormat<-cor(data)}
    
    critical.r <- function(nrow, a){
        df <- nrow - 2
        critical.t <- qt( a/2, df, lower.tail = F )
        cvr <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
        return(cvr)}
    
    if(thresh=="alpha")
    {
        if(missing(a))
        {a<-.05}
        
        thr<-critical.r(nrow(data),a)
        cormat<-ifelse(abs(cormat)>=thr,cormat,0)
    }else if(thresh=="adaptive")
    {
        if(missing(a))
        {a<-.05}
        
        a <- adapt.a(test = "cor", a, nrow(data), ...)
        
        thr<-critical.r(nrow(data),a)
        cormat<-ifelse(abs(cormat)>=thr,cormat,0)
    }else if(thresh=="bonferroni")
    {
        if(missing(a))
        {a<-.05}
        
        thr<-critical.r(nrow(data),(a/((ncol(cormat)^2)-(ncol(cormat))/2)))
        cormat<-ifelse(abs(cormat)>=thr,cormat,0)
    }else if(thresh=="FDR")
    {
        if(missing(a))
        {a<-.10}
        
        corsig<-cormat
        for(i in 1:ncol(data))
            for(j in 1:ncol(data))
                corsig[i,j]<-cor.test(data[,i],data[,j])$p.value
        
        fdrmat<-matrix(0,nrow=((ncol(cormat)^2)-(ncol(cormat))),ncol=3)
        fdrvec<-as.vector(corsig)
        
        fdrvec<-fdrtool::fdrtool(fdrvec,plot=FALSE,verbose=FALSE,statistic = "pvalue",...)$qval
        fdrvec<-ifelse(fdrvec<=a,fdrvec,0)
        corsig<-matrix(fdrvec,nrow=nrow(cormat),ncol(cormat))
        cormat<-ifelse(corsig!=0,cormat,0)
        thr<-min(cormat[cormat!=0])
    }else if(thresh=="proportional")
    {
        if(missing(a))
        {a<-.15}
        
        ltri<-cormat[lower.tri(cormat)]
        altri<-abs(ltri)
        
        poss<-(ncol(cormat)^2-ncol(cormat))/2
        
        den<-round(poss*a,0)
        
        thr<-min(sort(altri,decreasing=TRUE)[1:den])
        
        ftri<-ifelse(altri>=thr,ltri,0)
        
        cormat[lower.tri(cormat)]<-ifelse(ftri!=0,ltri,0)
        cormat[upper.tri(cormat)]<-0
        cormat<-cormat+t(cormat)
        
    }else if(is.numeric(thresh))
    {
        thr<-thresh
        cormat<-ifelse(abs(cormat)>=thr,cormat,0)
    }
    
    diag(cormat)<-1
    
    return(list(A=cormat, r.cv=thr))
}
#----