#' Triangulated Maximally Filtered Graph
#' @description Applies the Triangulated Maximally Filtered Graph (TMFG) filtering method
#' (\strong{Please see and cite Massara et al., 2016}). The TMFG method uses a structural
#' constraint that limits the number of zero-order correlations included in the network
#' (3n - 6; where \emph{n} is the number of variables). The TMFG algorithm begins by
#' identifying four variables which have the largest sum of correlations to all other
#' variables. Then, it iteratively adds each variable with the largest sum of three
#' correlations to nodes already in the network until all variables have been added to
#' the network. This structure can be associated with the inverse correlation matrix
#' (i.e., precision matrix) to be turned into a GGM (i.e., partial correlation network)
#' by using \code{\link[NetworkToolbox]{LoGo}}.
#' 
#' @param data Can be a dataset or a correlation matrix
#' 
#' @param normal Should data be transformed to a normal distribution?
#' Input must be a dataset.
#' Defaults to \code{FALSE}.
#' Data is not transformed to be normal.
#' Set to \code{TRUE} if data should be transformed to be normal
#' (computes correlations using the \code{\link[qgraph]{cor_auto}} function)
#' 
#' @param na.data How should missing data be handled?
#' For \code{"listwise"} deletion the \code{\link{na.omit}} function is applied.
#' Set to \code{"fiml"} for Full Information Maxmimum Likelihood (\code{\link[psych]{corFiml}}).
#' Full Information Maxmimum Likelihood is \strong{recommended} but time consuming
#' 
#' @param depend Is network a dependency (or directed) network?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} to generate a TMFG-filtered dependency network
#' (output obtained from the \code{\link[NetworkToolbox]{depend}} function)
#' 
#' @return Returns a list containing:
#' 
#' \item{A}{The filtered adjacency matrix}
#' 
#' \item{separators}{The separators (3-cliques) in the network
#' (wrapper output for \code{\link[NetworkToolbox]{LoGo}})}
#' 
#' \item{cliques}{The cliques (4-cliques) in the network
#' (wrapper output for \code{\link[NetworkToolbox]{LoGo}})}
#' 
#' @examples
#' TMFG.net <- TMFG(neoOpen)
#' 
#' @references
#' Christensen, A. P., Kenett, Y. N., Aste, T., Silvia, P. J., & Kwapil, T. R. (2018).
#' Network structure of the Wisconsin Schizotypy Scales-Short Forms: Examining psychometric network filtering approaches.
#' \emph{Behavior Research Methods}, 1-20.
#' doi: \href{https://doi.org/10.3758/s13428-018-1032-9}{10.3758/s13428-018-1032-9}
#' 
#' Massara, G. P., Di Matteo, T., & Aste, T. (2016).
#' Network filtering for big data: Triangulated maximally filtered graph.
#' \emph{Journal of Complex Networks}, \emph{5}, 161-178.
#' doi: \href{https://doi.org/10.1093/comnet/cnw015}{10.1093/comnet/cnw015}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @importFrom stats cor sd runif qt na.action qchisq
#' @importFrom utils capture.output
#' 
#' @export
#TMFG Filtering Method----
TMFG <-function (data, normal = FALSE,
                 na.data = c("pairwise","listwise","fiml","none"),
                 depend = FALSE)
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
    
    n<-ncol(cormat)
    
    tcormat<-cormat
    #cormat<-abs(cormat)
    
    if(n<9){print("Matrix is too small")}
    #nodeTO<-array()
    #nodeFROM<-array()
    #nodeWEIGHT<-array()
    #count<-0
    #for(i in 1:nrow(cormat))
    #    for(j in 1:ncol(cormat))
    #        if(cormat[i,j] != 0)
    #        {
    #            count<-count+1
    #            nodeTO[count]<-i
    #            nodeFROM[count]<-j
    #            nodeWEIGHT[count]<-cormat[i,j]
    #        }
    ###vectorized improvement
    nodeTO<-sort(c(rep(1:n,n)))
    nodeFROM<-c(rep(1:n,n))
    nodeWEIGHT<-as.vector(cormat)
    
    M<-cbind(nodeTO,nodeFROM,nodeWEIGHT) #create node-weight matrix
    in_v<-matrix(nrow=nrow(cormat),ncol=1) #initialize inserted vertices
    ou_v<-matrix(nrow=nrow(cormat),ncol=1) #initialize not yet inserted vertices
    tri<-matrix(nrow=((2*n)-4),ncol=3) #initializaes triangles
    separators<-matrix(nrow=n-4,ncol=3)#initialize list of 3-cliques (non-face triangles)
    #find 3 vertices with largest strength
    #s<-colSums(cormat*(cormat>mean(matrix(cormat,nrow=1)))*1) ##old s
    s<-rowSums(cormat*(cormat>mean(matrix(unlist(cormat),nrow=1)))*1)
    in_v[1:4]<-order(s,decreasing=TRUE)[1:4]
    ou_v<-setdiff(1:nrow(in_v),in_v)
    #build tetrahedron with the largest strength
    tri[1,]<-in_v[1:3,]
    tri[2,]<-in_v[2:4,]
    tri[3,]<-in_v[c(1,2,4),]
    tri[4,]<-in_v[c(1,3,4),]
    S<-matrix(nrow=(3*nrow(cormat)-6),ncol=3) #initialize sparse matrix
    if(!depend){S[1,]<-c(in_v[1],in_v[2],1)
    S[2,]<-c(in_v[1],in_v[3],1)
    S[3,]<-c(in_v[1],in_v[4],1)
    S[4,]<-c(in_v[2],in_v[3],1)
    S[5,]<-c(in_v[2],in_v[4],1)
    S[6,]<-c(in_v[3],in_v[4],1)
    }else{if(cormat[in_v[1],in_v[2]]>cormat[in_v[2],in_v[1]])
    {S[1,]<-c(in_v[1],in_v[2],1)
    }else{S[1,]<-c(in_v[2],in_v[1],1)}
        if(cormat[in_v[1],in_v[3]]>cormat[in_v[3],in_v[1]])
        {S[2,]<-c(in_v[1],in_v[3],1)
        }else{S[2,]<-c(in_v[3],in_v[1],1)}
        if(cormat[in_v[1],in_v[4]]>cormat[in_v[4],in_v[1]])
        {S[3,]<-c(in_v[1],in_v[4],1)
        }else{S[3,]<-c(in_v[4],in_v[1],1)}
        if(cormat[in_v[2],in_v[3]]>cormat[in_v[3],in_v[2]])
        {S[4,]<-c(in_v[2],in_v[3],1)
        }else{S[4,]<-c(in_v[3],in_v[2],1)}
        if(cormat[in_v[2],in_v[4]]>cormat[in_v[4],in_v[2]])
        {S[5,]<-c(in_v[2],in_v[4],1)
        }else{S[5,]<-c(in_v[4],in_v[2],1)}
        if(cormat[in_v[3],in_v[4]]>cormat[in_v[4],in_v[3]])
        {S[6,]<-c(in_v[3],in_v[4],1)
        }else{S[6,]<-c(in_v[4],in_v[3],1)}
    }
    #build initial gain table
    gain<-matrix(-Inf,nrow=n,ncol=(2*(n-2)))
    gain[ou_v,1]<-rowSums(cormat[ou_v,(tri[1,])])
    gain[ou_v,2]<-rowSums(cormat[ou_v,(tri[2,])])
    gain[ou_v,3]<-rowSums(cormat[ou_v,(tri[3,])])
    gain[ou_v,4]<-rowSums(cormat[ou_v,(tri[4,])])
    
    ntri<-4 #number of triangles
    gij<-matrix(nrow=1,ncol=ncol(gain))
    v<-matrix(nrow=1,ncol=ncol(gain))
    ve<-array()
    tr<-0
    for(e in 5:n)
    {
        if(length(ou_v)==1){
            ve<-ou_v
            v<-1
            w<-1
            tr<-which.max(gain[ou_v,])
        }else{
            for(q in 1:ncol(gain))
            {
                gij[,q]<-max(gain[ou_v,q])
                v[,q]<-which.max(gain[ou_v,q])
                tr<-which.max(gij)
            }
            ve<-ou_v[v[tr]]
            w<-v[tr]
        }
        #update vertex lists
        ou_v<-ou_v[-w]
        in_v[e]<-ve
        #update adjacency matrix
        for(u in 1:length(tri[tr,]))
        {
            cou<-6+((3*(e-5))+u)
            if(depend){
                if(cormat[ve,tri[tr,u]]>cormat[tri[tr,u],ve]){
                    S[cou,]<-cbind(ve,tri[tr,u],1)   
                }else{S[cou,]<-cbind(tri[tr,u],ve,1)}}else
                    S[cou,]<-cbind(ve,tri[tr,u],1)
        }
        #update 3-clique list
        separators[e-4,]<-tri[tr,]
        #update triangle list replacing 1 and adding 2 triangles
        tri[ntri+1,]<-cbind(rbind(tri[tr,c(1,3)]),ve)
        tri[ntri+2,]<-cbind(rbind(tri[tr,c(2,3)]),ve)
        tri[tr,]<-cbind(rbind(tri[tr,c(1,2)]),ve)
        #update gain table
        gain[ve,]<-0
        gain[ou_v,tr]<-rowSums(cormat[ou_v,tri[tr,],drop=FALSE])
        gain[ou_v,ntri+1]<-rowSums(cormat[ou_v,tri[ntri+1,],drop=FALSE])
        gain[ou_v,ntri+2]<-rowSums(cormat[ou_v,tri[ntri+2,],drop=FALSE])
        
        #update triangles
        ntri<-ntri+2
    }
    cliques<-rbind(in_v[1:4],(cbind(separators,in_v[5:ncol(cormat)])))
    
    L<-S
    if(depend)
    {W<-matrix(1:nrow(cormat),nrow=nrow(cormat),ncol=1)
    X<-matrix(1:nrow(cormat),nrow=nrow(cormat),ncol=1)
    Y<-matrix(0,nrow=nrow(cormat),ncol=1)
    Z<-cbind(W,X,Y)
    K<-rbind(L,Z)
    }else{
        L[,1]<-S[,2]
        L[,2]<-S[,1]
        K<-rbind(S,L)
    }
    
    x<-as.matrix(Matrix::sparseMatrix(i=K[,1],j=K[,2],x=K[,3]))
    diag(x)<-1
    
        for(r in 1:nrow(x))
            for(z in 1:ncol(x))
            {if(x[r,z]==1){x[r,z]<-tcormat[r,z]}}
    
    x<-as.matrix(x)
    colnames(x)<-colnames(cormat)
    rownames(x)<-colnames(cormat)
    
    return(list(A=x, separators=separators, cliques=cliques))
}
#----