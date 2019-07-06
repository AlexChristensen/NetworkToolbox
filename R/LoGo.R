#' Local/Global Inversion Method
#' 
#' @description Applies the Local/Global method to estimate
#' a Gaussian Graphical Model (GGM) using a \code{\link[NetworkToolbox]{TMFG}}-filtered network
#' (\strong{see and cite Barfuss et al., 2016}). Also used to 
#' convert clique and separator structure from
#' \code{\link[NetworkToolbox]{MFCF}} into partial correlation
#' and precision matrices
#' 
#' @param data Must be a dataset
#' 
#' @param cliques Cliques defined in the network.
#' Input can be a list or matrix
#' 
#' @param separators Separators defined in the network.
#' Input can be a list or matrix
#' 
#' @param normal Should data be transformed to a normal distribution?
#' Defaults to \code{TRUE} (computes correlations using the \code{\link[qgraph]{cor_auto}} function).
#' Set to \code{FALSE} for Pearson's correlations
#' 
#' @param na.data How should missing data be handled?
#' For \code{"listwise"} deletion the \code{\link{na.omit}} function is applied.
#' Set to \code{"fiml"} for Full Information Maxmimum Likelihood (\code{\link[psych]{corFiml}}).
#' Full Information Maxmimum Likelihood is \strong{recommended} but time consuming
#' 
#' @param partial Should the output network's connections be the partial correlation between two nodes given all other nodes?
#' Defaults to \code{TRUE}, which returns a partial correlation matrix.
#' Set to \code{FALSE} for a sparse inverse covariance matrix
#' 
#' @param ... Additional arguments (deprecated arguments)
#' 
#' @return Returns the sparse LoGo-filtered inverse covariance matrix (\code{partial = FALSE})
#' or LoGo-filtered partial correlation matrix (\code{partial = TRUE})
#' 
#' @examples
#' \dontrun{
#' LoGonet <- LoGo(neoOpen, partial = TRUE)
#' }
#' 
#' @references
#' Barfuss, W., Massara, G. P., Di Matteo, T., & Aste, T. (2016).
#' Parsimonious modeling with information filtering networks.
#' \emph{Physical Review E}, \emph{94}, 062306.
#' doi: \href{https://doi.org/10.1103/PhysRevE.94.062306}{10.1103/PhysRevE.94.062306}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @importFrom stats cov
#' @importFrom psych corFiml
#' 
#' @export
#LoGo Sparse Inverse Covariance Matrix----
LoGo <- function (data, cliques, separators,
                  normal = TRUE, 
                  na.data = c("pairwise","listwise","fiml","none"),
                  partial = TRUE, ...)
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
    }else if(na.data=="none")
    {
        if(nrow(data)==ncol(data)){cormat<-data
        }else if(normal){cormat<-qgraph::cor_auto(data)
        }else{cormat<-cor(data)}
    }
    
    #covariance matrix
    standardize <- TRUE
    
    
    S <- cormat
    
    if(missing(separators))
    {separators<-NULL}
    
    if(missing(cliques))
    {cliques<-NULL}
    
    
    if(is.null(separators)&is.null(cliques))
    {
        tmfg<-TMFG(cormat)
        separators<-tmfg$separators
        cliques<-tmfg$cliques
    }
    
    n<-ncol(S)
    Jlogo<-matrix(0,nrow=n,ncol=n)
    
    if(!is.list(cliques)&!is.list(separators))
    {
        for(i in 1:nrow(cliques))
        {
            v<-cliques[i,]
            Jlogo[v,v]<-Jlogo[v,v]+solve(S[v,v])
        }
        
        for(i in 1:nrow(separators))
        {
            v<-separators[i,]
            Jlogo[v,v]<-Jlogo[v,v]-solve(S[v,v])
        }
    }else{
        for(i in 1:length(cliques))
        {
            v<-cliques[[i]]
            Jlogo[v,v]<-Jlogo[v,v]+solve(S[v,v])
        }
        
        for(i in 1:length(separators))
        {
            v<-separators[[i]]
            Jlogo[v,v]<-Jlogo[v,v]-solve(S[v,v])
        }
    }
    
    if(partial)
    {
        Jlogo<-(-cov2cor(Jlogo))
        if(any(is.na(Jlogo)))
        {Jlogo <- ifelse(is.na(Jlogo),0,Jlogo)}
        diag(Jlogo)<-0
    }
    
    colnames(Jlogo)<-colnames(data)
    row.names(Jlogo)<-colnames(data)
    
    if(!isSymmetric(Jlogo))
    {Jlogo<-as.matrix(Matrix::forceSymmetric(Jlogo))}
    
    return(logo=Jlogo)
}
#----