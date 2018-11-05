#' Bootstrapped Network Generalization
#' @description Bootstraps the sample to identify the most stable correlations.
#' Also produces a network that is penalizes low reliability edges.
#' This function is useful for overcoming the structural constraint of the IFN approach.
#' STILL BEING DEVELOPED
#' 
#' @param data A set of data
#' 
#' @param method A network filtering method.
#' Defaults to "TMFG"
#' 
#' @param n Number of people to use in the bootstrap.
#' Defaults to full sample size
#' 
#' @param iter Number of bootstrap iterations.
#' Defaults to 1000 iterations
#' 
#' @param normal Should data be transformed to a normal distribution?
#' Defaults to FALSE.
#' Data is not transformed to be normal.
#' Set to TRUE if data should be transformed to be normal
#' (computes correlations using the \link[qgraph]{cor_auto} function)
#' 
#' @param na.data How should missing data be handled?
#' For "listwise" deletion the \code{\link{na.omit}} function is applied.
#' Set to "fiml" for Full Information Maxmimum Likelihood (\link[psych]{corFiml}).
#' Full Information Maxmimum Likelihood is \strong{recommended} but time consuming
#' 
#' @param cores Number of computer processing cores to use for bootstrapping samples.
#' Defaults to \emph{n} - 1 total number of cores.
#' Set to any number between 1 and maxmimum amount of cores on your computer
#' 
#' @param ... Additional arguments for filtering methods
#' 
#' @return Returns a list that includes:
#' 
#' \item{orignet}{the original filtered network}
#' 
#' \item{bootmat}{correlation matrix of the mean bootstrapped network}
#' 
#' \item{netrel}{unfiltered reliabilities of all of the connections}
#' 
#' @examples
#' \dontrun{
#' bootTMFG<-bootgen(neoOpen)
#' 
#' bootLoGo<-bootgen(neoOpen,method="LoGo")
#' 
#' bootMaST<-bootgen(neoOpen,method="MaST")
#' 
#' bootThreshold<-bootgen(neoOpen,method="threshold")
#' }
#' 
#' @references
#' Musciotto, F., Marotta, L., Micciche, S., & Mantegna, R. N. (2018).
#' Bootstrap validation of links of a minimum spanning tree.
#' \emph{arXiv}, 1802.03395.
#' doi: \href{https://arxiv.org/pdf/1802.03395.pdf}{1802.03395}
#' 
#' Tumminello, M., Coronnello, C., Lillo, F., Micciche, S., & Mantegna, R. N. (2007).
#' Spanning trees and bootstrap reliability estimation in correlation-based networks.
#' \emph{International Journal of Bifurcation and Chaos}, \emph{17}, 2319-2329.
#' doi: \href{https://doi.org/10.1142/S0218127407018415}{10.1142/S0218127407018415}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @importFrom stats cov2cor pt
#' 
#' @export
#Bootstrap Network Generalization----
bootgen <- function (data, method = c("MaST", "TMFG", "LoGo", "threshold"),
                     n = nrow(data), iter = 1000, normal = FALSE,
                     na.data = c("pairwise", "listwise", "fiml","none"),
                     cores, ...)
{
    
    #arguments
    ##########################################################
    if(missing(method))
    {method<-"TMFG"
    }else{method<-match.arg(method)}
    
    ##########################################################
    #fisher z and alpha functions
    ##########################################################
    #fisher's z
    fish <- function (r)
    {z<-.5*log((1+abs(r))/(1-abs(r)))
    if(nrow(r)>1&&ncol(r)>1&&length(r)>1)
    {diag(z)<-0}
    for(i in 1:ncol(r))
        for(j in 1:nrow(r))
        {
            if(r[i,j]<0)
            {z[i,j]<--z[i,j]}
        }
    return(z)}
    
    #mean fisher's z
    zw <- function (z,iter)
    {sum((iter-3)*(z))/((iter-3)*iter)}
    
    ##########################################################
    
    #missing data handling
    ##########################################################
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
    ##########################################################
    
    if(method=="LoGo")
    {
        invmat<-solve(cov(data))
        realmat<--cov2cor(invmat)
        diag(realmat)<-0
    }else{realmat<-cormat}
    
    #general bootstrapping    
    ##########################################################
    sampslist<-list() #initialize sample list
    
    #Parallel processing
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    
    sampslist<-foreach::foreach(i=1:iter,
                                .packages = c("NetworkToolbox","psych","qgraph"))%dopar%
                                {
                                    mat<-data[sample(1:n,replace=TRUE),]
                                    
                                    #missing data handling
                                    if(any(is.na(mat)))
                                    {
                                        if(na.data=="pairwise")
                                        {
                                            if(normal)
                                            {cormat<-qgraph::cor_auto(mat,missing=na.data)
                                            }else{cormat<-cor(mat,use="pairwise.complete.obs")}
                                        }else if(na.data=="listwise")
                                        {
                                            if(normal)
                                            {cormat<-qgraph::cor_auto(mat,missing=na.data)
                                            }else{mat<-na.omit(mat)
                                            }
                                        }else if(na.data=="fiml")
                                        {
                                            if(normal)
                                            {cormat<-qgraph::cor_auto(mat,missing=na.data)
                                            }else{cormat<-psych::corFiml(mat)}
                                        }
                                    }else{
                                        if(nrow(mat)==ncol(mat)){cormat<-mat
                                        }else if(normal){cormat<-qgraph::cor_auto(mat)
                                        }else{cormat<-cor(mat)}
                                    }
                                    
                                    if(method=="TMFG")
                                    {samps<-fish(TMFG(cormat)$A)
                                    }else if(method=="LoGo")
                                    {samps<-fish(suppressWarnings(LoGo(mat,partial=TRUE)))
                                    }else if(method=="MaST")
                                    {samps<-fish(MaST(cormat,...))
                                    }else if(method=="threshold")
                                    {samps<-fish(threshold(cormat,...)$A)
                                    }else{stop("Method not available")}
                                    
                                    return(samps)
                                }
    parallel::stopCluster(cl)
    ##########################################################
    
    #Original Networks
    if(method=="TMFG")
    {tru<-TMFG(data)$A
    }else if(method=="LoGo")
    {tru<-LoGo(data,partial=TRUE)
    diag(tru)<-0
    }else if(method=="MaST")
    {tru<-MaST(data,...)
    }else if(method=="MaST")
    {tru<-MaST(data,...)
    }else if(method=="threshold")
    {tru<-threshold(data,...)$A
    }else stop("Method not available")
    
    diag(tru)<-0
    
    #convert samples list from foreach to array
    samps<-array(0,dim=c(nrow=nrow(realmat),ncol=ncol(realmat),iter))
    
    for(i in 1:iter) #populate array
    {samps[,,i]<-sampslist[[i]]}
    
    #Mean fisher matrix
    meanmat<-matrix(0,nrow=nrow(realmat),ncol=ncol(realmat)) #Initialize Mean matrix
    for(j in 1:nrow(realmat))
        for(k in 1:ncol(realmat))
        {meanmat[j,k]<-zw(samps[j,k,],iter)}
    
    #convert fisher z to r
    meanmat<-psych::fisherz2r(meanmat)
    
    #ensure model is graphical
    if(method=="LoGo")
    {
        if(!is.graphical(bootmat))
        {
            boot<-abs(bootmat[lower.tri(bootmat)])
            real<-abs(meanmat[lower.tri(meanmat)])
            comb<-as.data.frame(cbind(real,boot))
            ordered<-as.data.frame(comb[order(comb[,1],decreasing = TRUE),])
            
            zeros<-which(ordered$real!=0&ordered$boot==0)
            
            for(i in 1:length(zeros))
            {
                bootmat[which(round(abs(meanmat),5)==round(ordered[zeros[i],]$real,5))]<-meanmat[which(round(abs(meanmat),5)==round(ordered[zeros[i],]$real,5))]
                
                if(is.graphical(bootmat)){break}
            }
        }
    }else{bootmat <- meanmat}
    
    #return bootmat
    bootmat<-ifelse(bootmat!=0,realmat,0)
    colnames(bootmat)<-colnames(realmat)
    
    diag(bootmat)<-0
    
    if(method=="LoGo")
    {
        diag(bootmat)<-1
        
        invmat<-ifelse(bootmat!=0,invmat,0)
    }
    
    if(!isSymmetric(bootmat))
    {bootmat<-as.matrix(Matrix::forceSymmetric(bootmat))}
    
    #Reliability matrix
    samp<-array(0,c(nrow=nrow(realmat),ncol=ncol(realmat),iter))
    
    rel<-matrix(0,nrow=nrow(realmat),ncol=ncol(realmat))
    for(j in 1:nrow(realmat))
        for(k in 1:ncol(realmat))
            for(l in 1:iter)
                if(samps[j,k,l]!=0)
                {samp[j,k,l]<-1}
    
    #reliablity plot
    for(j in 1:nrow(realmat))
        for(k in 1:ncol(realmat))
            rel[j,k]<-sum(samp[j,k,])/iter
    colnames(rel)<-colnames(data)
    
    if(method=="LoGo")
    {return(list(orignet=tru,bootmat=bootmat,netrel=rel,invcov=invmat))
    }else{return(list(orignet=tru,bootmat=bootmat,netrel=rel))}
}
#----