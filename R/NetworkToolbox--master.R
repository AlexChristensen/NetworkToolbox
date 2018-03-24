#  NetworkToolbox Functions
#
#' Triangulated Maximally Filtered Graph
#' @description Applies the Triangulated Maximally Filtered Graph (TMFG) filtering method
#' (Please see and cite Massara et al., 2016)
#' @param data Can be a dataset or a correlation matrix
#' @param normal Should data be transformed to a normal distribution?
#' Input must be a dataset.
#' Defaults to FALSE.
#' Data is not transformed to be normal.
#' Set to TRUE if data should be transformed to be normal
#' (computes correlations using the \emph{cor_auto} function from the \emph{qgraph} package)
#' @param weighted Should network be weighted?
#' Defaults to TRUE.
#' Set to FALSE to produce an unweighted (binary) network
#' @param depend Is network a dependency (or directed) network?
#' Defaults to FALSE.
#' Set to TRUE to generate a TMFG-filtered dependency network
#' (output obtained from the \emph{depend} function)
#' @param na.data How should missing data be handled?
#' For "listwise" deletion the \emph{na.omit} function is applied.
#' Set to "fiml" for Full Information Maxmimum Likelihood (\emph{psych} package).
#' Full Information Maxmimum Likelihood is \strong{recommended} but time consuming
#' @return Returns a list of the adjacency matrix (A), separators (separators), and cliques (cliques)
#' @examples
#' weighted_TMFGnetwork<-TMFG(neoOpen)
#' 
#' 
#' unweighted_TMFGnetwork<-TMFG(neoOpen,weighted=FALSE)
#' 
#' @references
#' Christensen, A. P., Kenett, Y. N., Aste, T., Silvia, P. J., & Kwapil, T. R. (2018).
#' Network structure of the Wisconsin Schizotypy Scales-Short Forms: Examining psychometric network filtering approaches.
#' \emph{Behavior Research Methods}, 1-20, doi:10.3758/s13428-018-1032-9
#' 
#' Massara, G. P., Di Matteo, T., & Aste, T. (2016).
#' Network filtering for big data: Triangulated maximally filtered graph.
#' \emph{Journal of Complex Networks}, \emph{5}(2), 161-178.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @importFrom stats cor sd runif qt na.action qchisq
#' @importFrom utils capture.output
#' @export
#TMFG Filtering Method----
TMFG <-function (data, normal = FALSE, weighted = TRUE, depend = FALSE,
                 na.data = c("pairwise","listwise","fiml","none"))
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
        }else{data<-psych::corFiml(data)}
    }
    
    #corrlation matrix
    if(nrow(data)==ncol(data)){cormat<-data
    }else if(normal){cormat<-qgraph::cor_auto(data)
    }else{cormat<-cor(data)}
    
    n<-ncol(cormat)
    
        tcormat<-cormat
        cormat<-abs(cormat)
        
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
    
    if(weighted)
    {
        for(r in 1:nrow(x))
            for(z in 1:ncol(x))
            {if(x[r,z]==1){x[r,z]<-tcormat[r,z]}
            }
    }
    
    x<-as.matrix(x)
    colnames(x)<-colnames(cormat)
    rownames(x)<-colnames(cormat)
    
    return(list(A=x, separators=separators, cliques=cliques))
}
#----
#' Local/Global Sparse Inverse Covariance Matrix
#' @description Applies the Local/Global method to estimate the sparse inverse covariance matrix using a TMFG-filtered network
#' (see and cite Barfuss et al., 2016)
#' @param data Must be a dataset
#' @param partial Should the output network's connections be the partial correlation between two nodes given all other nodes?
#' Defaults to FALSE, which returns a sparse inverse covariance matrix.
#' Set to TRUE for a partial correlation matrix
#' @param normal Should data be transformed to a normal distribution?
#' Defaults to FALSE.
#' Data is not transformed to be normal.
#' Set to TRUE if data should be transformed to be normal
#' (computes correlations using the \emph{cor_auto} function from the \emph{qgraph} package)
#' @param na.data How should missing data be handled?
#' For "listwise" deletion the \emph{na.omit} function is applied.
#' Set to "fiml" for Full Information Maxmimum Likelihood (\emph{psych} package).
#' Full Information Maxmimum Likelihood is \strong{recommended} but time consuming
#' @param graphical Should network be checked for graphical modeling?
#' Defaults to FALSE.
#' Set to TRUE to determine if model is graphical
#' @return Returns a list containing the TMFG-filtered matrix (tmfg),
#' and the sparse TMFG-filtered inverse covariance matrix (logo)
#' @examples
#' LoGonet<-LoGo(neoOpen, partial = TRUE)$logo
#' 
#' TMFGnet<-LoGo(neoOpen)$tmfg
#' @references 
#' Barfuss, W., Massara, G. P., Di Matteo, T., & Aste, T. (2016).
#' Parsimonious modeling with information filtering networks.
#' \emph{Physical Review E}, \emph{94}(6), 062306.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @importFrom stats cov
#' @importFrom psych corFiml
#' @export
#LoGo Sparse Inverse Covariance Matrix----
LoGo <- function (data, partial = FALSE, normal = FALSE, 
                  na.data = c("pairwise","listwise","fiml","none"), graphical = FALSE)
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
        }else{data<-psych::corFiml(data)}
    }
    
    #corrlation matrix
    if(nrow(data)==ncol(data)){cormat<-data
    }else if(normal){cormat<-qgraph::cor_auto(data)
    }else{cormat<-cor(data)}
    
    #corrlation matrix
    if(nrow(data)==ncol(data)){
        cormat<-data
        S<-cor2cov(cormat,data)
    }else if(normal){
        #shrinkage covariance
        S<-cor2cov(cormat,data)
    }else{
        cormat<-cor(data)
        S<-cov(data)
    }
    
    tmfg<-TMFG(cormat)
    separators<-tmfg$separators
    cliques<-tmfg$cliques
    
    n<-ncol(S)
    Jlogo<-matrix(0,nrow=n,ncol=n)
    
    for(i in 1:nrow(cliques))
    {v<-cliques[i,]
    Jlogo[v,v]<-Jlogo[v,v]+solve(S[v,v])}
    
        for(i in 1:nrow(separators))
        {v<-separators[i,]
        Jlogo[v,v]<-Jlogo[v,v]-solve(S[v,v])}
    
    if(graphical)
    {
        par<-(-cov2cor(Jlogo))
        diag(par)<-0
        is.graphical(par)
    }
    
        if(partial)
        {
            Jlogo<-(-cov2cor(Jlogo))
            diag(Jlogo)<-0
            
            if(graphical)
            {is.graphical(Jlogo)}
        }
        
    colnames(Jlogo)<-colnames(data)
    row.names(Jlogo)<-colnames(data)
    
    if(!isSymmetric(Jlogo))
    {Jlogo<-as.matrix(Matrix::forceSymmetric(Jlogo))}
    
    return(list(logo=Jlogo,tmfg=tmfg$A))
}
#----
#' Planar Maximally Filtered Graph
#' @description Applies the Planar Maximally Filtered Graph (PMFG) filtering method
#' (see and cite Tumminello et al., 2005).
#' Currently very slow! (efficiency is being worked on)
#' @param data Can be a dataset or a correlation matrix
#' @param sparseList Should a sparse matrix be returned?
#' Defaults to FALSE.
#' Set to TRUE to return a sparse matrix list
#' @param normal Should data be transformed to a normal distribution?
#' Defaults to FALSE.
#' Data is not transformed to be normal.
#' Set to TRUE if data should be transformed to be normal
#' (computes correlations using the \emph{cor_auto} function from the \emph{qgraph} package)
#' @param weighted Should network be weighted?
#' Defaults to TRUE.
#' Set to FALSE to produce an unweighted (binary) network
#' @param na.data How should missing data be handled?
#' For "listwise" deletion the \emph{na.omit} function is applied.
#' Set to "fiml" for Full Information Maxmimum Likelihood (\emph{psych} package).
#' Full Information Maxmimum Likelihood is \strong{recommended} but time consuming
#' @param progBar Should progress bar be displayed?
#' Defaults to TRUE.
#' Set to FALSE for no progress bar
#' @return Returns a list of PMFG-filtered associaton matrix (pmfg).
#' If \emph{sparseList} = TRUE, then the list will contain a sparse list (List)
#' @examples
#' \dontrun{
#' weighted_PMFGnetwork<-PMFG(neoOpen)
#' }
#' @references
#' Tumminello, M., Aste, T., Di Matteo, T., & Mantegna, R. N. (2005).
#' A tool for filtering information in complex systems.
#' \emph{Proceedings of the National Academy of Sciences}, \emph{102}(30), 10421-10426.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#PMFG Filtering Method----
PMFG <- function (data, sparseList = FALSE, normal = FALSE, weighted = TRUE,
                  na.data = c("pairwise","listwise","fiml","none"), progBar = TRUE)
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
        }else{data<-psych::corFiml(data)}
    }
    
    #corrlation matrix
    if(nrow(data)==ncol(data)){cormat<-data
    }else if(normal){cormat<-qgraph::cor_auto(data)
    }else{cormat<-cor(data)}
    
    #create sparse data
    i<-as.vector(rep(1:ncol(data),ncol(data)))
    j<-sort(as.vector(rep(1:ncol(data),ncol(data))))
    w<-as.vector(cormat)
    
    kk<-which(i<j)
    
    ijw<-cbind(i[kk],j[kk],w[kk])
    
    ijw<-ijw[order(ijw[,3],decreasing=TRUE),]
    
    
    P<-Matrix::Matrix(0,nrow=ncol(data),ncol=ncol(data))
    
    as_graphnel <- function(graph) {
        
        if (!igraph::is_igraph(graph)) {
            stop("Not an igraph graph")
        }
        
        if ("name" %in% suppressWarnings(igraph::vertex_attr_names(graph)) &&
            is.character(suppressWarnings(igraph::V(graph)$name))) {
            name <- suppressWarnings(igraph::V(graph)$name)
        } else {
            name <- as.character(seq(suppressWarnings(igraph::vcount(graph))))    
        }
        
        edgemode <- "undirected"  
        
        if ("weight" %in% suppressWarnings(igraph::edge_attr_names(graph)) &&
            is.numeric(suppressWarnings(igraph::E(graph)$weight))) {
            al <- suppressWarnings(lapply(igraph::as_adj_edge_list(graph, "out"), as.vector))
            for (i in seq(along=al)) {
                edges <- suppressWarnings(igraph::ends(graph, al[[i]], names = FALSE))
                edges <- ifelse( edges[,2]==i, edges[,1], edges[,2])
                weights <- suppressWarnings(igraph::E(graph)$weight[al[[i]]])
                al[[i]] <- list(edges=edges, weights=weights)
            }
        } else {
            al <- suppressWarnings(igraph::as_adj_list(graph, "out"))
            al <- lapply(al, function(x) list(edges=as.vector(x)))
        }  
        
        names(al) <- name
        res <- graph::graphNEL(nodes=name, edgeL=al, edgemode=edgemode)
        
        ## Add graph attributes (other than 'directed')
        ## Are this "officially" supported at all?
        
        g.n <- suppressWarnings(igraph::graph_attr_names(graph))
        if ("directed" %in% g.n) {
            warning("Cannot add graph attribute `directed'")
            g.n <- g.n[ g.n != "directed" ]
        }
        for (n in g.n) {
            res@graphData[[n]] <- suppressWarnings(igraph::graph_attr(graph, n))
        }
        
        ## Add vertex attributes (other than 'name', that is already
        ## added as vertex names)
        
        v.n <- suppressWarnings(igraph::vertex_attr_names(graph))
        v.n <- v.n[ v.n != "name" ]
        for (n in v.n) {
            graph::nodeDataDefaults(res, attr=n) <- NA
            graph::nodeData(res, attr=n) <- suppressWarnings(igraph::vertex_attr(graph, n))
        }
        
        ## Add edge attributes (other than 'weight')
        
        e.n <- suppressWarnings(igraph::edge_attr_names(graph))
        e.n <- e.n[ e.n != "weight" ]
        if (length(e.n) > 0) {
            el <- suppressWarnings(igraph::as_edgelist(graph))
            el <- paste(sep="|", el[,1], el[,2])
            for (n in e.n) {
                graph::edgeDataDefaults(res, attr=n) <- NA
                res@edgeData@data[el] <- mapply(function(x,y) {
                    xx <- c(x,y); names(xx)[length(xx)] <- n; xx },
                    res@edgeData@data[el],
                    suppressWarnings(igraph::edge_attr(graph, n)),
                    SIMPLIFY=FALSE)
            }
        }
        
        res
    }
    
    for(ii in 1:pmin(6,nrow(ijw)))
    {
        P[ijw[ii,1],ijw[ii,2]]<-ijw[ii,3]
        P[ijw[ii,2],ijw[ii,1]]<-ijw[ii,3]
    }
    
    E<-6
    P1<-P
    
    if(progBar==TRUE)
    {pb <- txtProgressBar(max=(3*(ncol(data)-2)), style = 3)}
    
    while(E < 3*(ncol(data)-2))
    {
        ii<-ii+1
        
        P1[ijw[ii,1],ijw[ii,2]]<-ijw[ii,3]
        P1[ijw[ii,2],ijw[ii,1]]<-ijw[ii,3]
        
        graph<-suppressWarnings(convert2igraph(P1))
        
        g<-as_graphnel(graph)
        
            if(RBGL::boyerMyrvoldPlanarityTest(g)==TRUE)
            {
                P<-P1
                E<-E+1
                
                if(progBar==TRUE)
                {setTxtProgressBar(pb, E)}
            }else{P1<-P}
        
        if(ii>(ncol(data)*(ncol(data)-1)/2))
        {print("PMFG solution not found")}
        
    }
        
    if(progBar==TRUE)
    {close(pb)}
    
    pmfg<-as.matrix(P)
    
    if(!weighted)
    {pmfg<-ifelse(pmfg!=0,1,0)}
    
    if(sparseList)
    {
        pmfg<-ifelse(pmfg!=0,1,0)
    
        i<-as.vector(rep(1:ncol(data),ncol(data)))
        j<-sort(as.vector(rep(1:ncol(data),ncol(data))))
        w<-as.vector(pmfg)
    
        List<-cbind(i,j,w)
        for(i in 1:nrow(List))
        if(List[i,3]==0)
        {List[i,3]<-NA}
        List<-na.omit(List)
        
        return(list(pmfg=pmfg, List=List))
    }else{return(list(pmfg=pmfg))}
}
#----
#' Maximum Spanning Tree
#' @description Applies the Maximum Spanning Tree (MaST) filtering method
#' @param data Can be a dataset or a correlation matrix
#' @param normal Should data be transformed to a normal distribution?
#' Defaults to FALSE.
#' Data is not transformed to be normal.
#' Set to TRUE if data should be transformed to be normal
#' (computes correlations using the \emph{cor_auto} function from the \emph{qgraph} package)
#' @param weighted Should network be weighted?
#' Defaults to TRUE.
#' Set to FALSE to produce an unweighted (binary) network
#' @param depend Is network a dependency (or directed) network?
#' Defaults to FALSE.
#' Set TRUE to generate a MaST-filtered dependency network
#' (output obtained from the \emph{depend} function)
#' @param na.data How should missing data be handled?
#' For "listwise" deletion the \emph{na.omit} function is applied.
#' Set to "fiml" for Full Information Maxmimum Likelihood (\emph{psych} package).
#' Full Information Maxmimum Likelihood is \strong{recommended} but time consuming
#' @return A sparse association matrix
#' @examples
#' weighted_MaSTnetwork<-MaST(neoOpen)
#' 
#' 
#' unweighted_MaSTnetwork<-MaST(neoOpen,weighted=FALSE)
#' 
#' @references 
#' Adapted from: \url{https://www.mathworks.com/matlabcentral/fileexchange/23276-maximum-weight-spanning-tree--undirected}
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Maximum Spanning Tree----
MaST <- function (data, normal = FALSE, weighted = TRUE,
                  depend = FALSE, na.data = c("pairwise","listwise","fiml","none"))
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
        }else{data<-psych::corFiml(data)}
    }
    
    #corrlation matrix
    if(nrow(data)==ncol(data)){cormat<-data
    }else if(normal){cormat<-qgraph::cor_auto(data)
    }else{cormat<-cor(data)}
    
FIND_PathCompression <- function (temproot)
{
  ParentPointer[temproot]
  if(ParentPointer[temproot]!=temproot)
  {ParentPointer[temproot]<-FIND_PathCompression(ParentPointer[temproot])}
  parent<-ParentPointer[temproot]
}

corma<-abs(cormat)
nodeT<-0
nodeF<-0
weights<-0
wc<-0
n<-ncol(corma)
for (i in 1:n)
  for (j in 1:n)
    if (corma[i,j]!=0) #Figure out how to remove warning
    {
      wc<- wc+1
      nodeT[wc] <- i
      nodeF[wc] <- j
      weights[wc] <- corma[i,j]
    }
edgelist<-cbind(weights,nodeT,nodeF)
edgelist<-edgelist[order(edgelist[,1]),]
#Number of edges
e <- nrow(edgelist)
#Assign ParentPointer to each vertex
assign("ParentPointer",1:n)
#Assign a tree rank to each vertex
TreeRank<-matrix(0,nrow=1,ncol=n)
#MSTreeEdges and counter
MSTreeEdges<-matrix(0,nrow=n-1,ncol=3)
MSTcounter<-0
i<-e

while((MSTcounter<(n-1))&&(e>=1))
{
  #Find roots of the tree that the selected edge's two
  #vertices belong to. Also perform path compression.
  root1<-0
  root2<-0
  temproot<-0
  temproot<-as.numeric(edgelist[i,2])
  root1<-FIND_PathCompression(temproot)
  
  temproot<-as.numeric(edgelist[i,3])
  root2<-FIND_PathCompression(temproot)
  
  if(root1!=root2)
  {
    MSTcounter<-MSTcounter+1
    MSTreeEdges[MSTcounter,1:3]<-edgelist[i,]
    if(TreeRank[root1]>TreeRank[root2])
    {ParentPointer[root2]<-root1}else 
      if(TreeRank[root1]==TreeRank[root2])
      {TreeRank[root2]<-TreeRank[root2]+1
      ParentPointer[root1]<-root2}else if(TreeRank[root1]<TreeRank[root2])
      {ParentPointer[root1]<-root2}
  }
  i<-i-1
}
S<-MSTreeEdges
if(depend)
{
    for(i in 1:nrow(S))
        if(cormat[S[i,2],S[i,3]]>=cormat[S[i,3],S[i,2]])
        {S[i,1]<-cormat[S[i,2],S[i,3]]
        S[i,2]<-S[i,2]
        S[i,3]<-S[i,3]
        }else if(cormat[S[i,2],S[i,3]]<cormat[S[i,3],S[i,2]])
        {
            S[i,1]<-cormat[S[i,3],S[i,2]]
            S[i,2]<-S[i,3]
            S[i,3]<-S[i,2] 
        }
}
L<-S
if(depend)
{W<-matrix(1:nrow(cormat),nrow=nrow(cormat),ncol=1)
X<-matrix(1:nrow(cormat),nrow=nrow(cormat),ncol=1)
Y<-matrix(0,nrow=nrow(cormat),ncol=1)
Z<-cbind(Y,W,X)
K<-rbind(L,Z)
}else{L[,2]<-S[,3]
L[,3]<-S[,2]
K<-rbind(S,L)}
x<-as.matrix(Matrix::sparseMatrix(i=K[,2],j=K[,3],x=K[,1]))
diag(x)<-1
x<-as.matrix(x)
ifelse(x!=0,cormat,0)
if(!weighted)
{x<-ifelse(x!=0,1,0)}
x<-as.data.frame(x)
colnames(x)<-colnames(cormat)
x<-as.matrix(x)
return(x)
}
#----
#' ECO Neural Network Filter
#' @description Applies the ECO neural network filtering method
#' @param data Can be a dataset or a correlation matrix
#' @param weighted Should network be weighted?
#' Defaults to TRUE.
#' Set to FALSE to produce an unweighted (binary) network
#' @param directed Is the network directed?
#' Defaults to FALSE.
#' Set TRUE if the network is directed
#' @return A sparse association matrix
#' @examples
#' \dontrun{
#' weighted_undirected_ECOnetwork<-ECO(neoOpen)
#' 
#' unweighted_undirected_ECOnetwork<-ECO(neoOpen,weighted=FALSE)
#' 
#' weighted_directed_ECOnetwork<-ECO(neoOpen,directed=TRUE)
#' 
#' unweighted_directed_ECOnetwork<-ECO(neoOpen,weighted=FALSE,directed=TRUE)
#' }
#' @references 
#' Fallani, F. D. V., Latora, V., & Chavez, M. (2017).
#' A topological criterion for filtering information in complex brain networks.
#' \emph{PLoS Computational Biology}, \emph{13}(1), e1005305.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#ECO Neural Network Filter----
ECO <- function (data, weighted = TRUE, directed = FALSE)
{
    #corrlation matrix
    if(nrow(data)==ncol(data)){cormat<-data
    }else{cormat<-cor(data)}
    
    C<-cormat
    
  n<-ncol(C)
  S<-C
  if(directed)
  {
    numcon<-3*n
    ind<-which(C!=0)
  }else{C<-upper.tri(C,diag=TRUE)
  numcon<-1.5*n
  ind<-which(upper.tri(C,diag=TRUE)!=0)}
  
  S<-ifelse(C==1,S,0)
  
  if(numcon>length(ind))
  {
    stop("Input matrix is too sparse")
  }
  
  sorind<-matrix(0,nrow=length(ind),ncol=2)
  
  G<-S
  S<-abs(S)
  
  x<-S[ind]
  y<-ind
  h<-cbind(ind,S[ind])
  sorind<-h[order(-h[,2]),]
  C[sorind[(numcon+1):nrow(sorind),1]]<-0
  
  if(directed)
  {W<-C}else{W<-C+t(C)
  diag(W)<-1}
  J<-G+t(G)
  diag(J)<-1
  if(weighted)
  {
    W<-ifelse(W!=0,J,0)
  }
  W<-as.data.frame(W)
  colnames(W)<-colnames(data)
  W<-as.matrix(W)
  return(W)
}
#----
#' ECO+MaST Network Filter
#' @description Applies the ECO neural network filtering method combined with the MaST filtering method
#' @param data Can be a dataset or a correlation matrix
#' @param weighted Should network be weighted?
#' Defaults to TRUE.
#' Set to FALSE to produce an unweighted (binary) network
#' @return A sparse association matrix
#' @examples
#' weighted_ECOplusMaSTnetwork<-ECOplusMaST(neoOpen)
#' @references 
#' Fallani, F. D. V., Latora, V., & Chavez, M. (2017).
#' A topological criterion for filtering information in complex brain networks.
#' \emph{PLoS Computational Biology}, \emph{13}(1), e1005305.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#ECO Filter + MaST----
ECOplusMaST <- function (data, weighted = TRUE)
{
    #corrlation matrix
    if(nrow(data)==ncol(data)){cormat<-data
    }else{cormat<-cor(data)}
    
  if(weighted)
{
  a<-MaST(data,weighted=TRUE)
  b<-ECO(data,weighted=TRUE)
  k<-matrix(NA,nrow=nrow(a),ncol=ncol(a))
  for(i in 1:nrow(a))
    for(j in 1:ncol(a))
      if(a[i,j]==b[i,j]){k[i,j]<-a[i,j]}else k[i,j]<-a[i,j]+b[i,j]
}else
{
  a<-MaST(data,weighted=FALSE)
  b<-ECO(data,weighted=FALSE)
  k<-matrix(NA,nrow=nrow(a),ncol=ncol(a))
  for(i in 1:nrow(a))
    for(j in 1:ncol(a))
      if(a[i,j]==b[i,j])
      {k[i,j]<-a[i,j]}else k[i,j]<-a[i,j]+b[i,j]
}
  k<-as.data.frame(k)
  colnames(k)<-colnames(data)
  k<-as.matrix(k)
  return(k)
}
#----
#' Threshold Filter
#' @description Filters the network based on an r-value, alpha, adaptive alpha (see Perez & Pericchi, 2014),
#' bonferroni, false-discovery rate (FDR, \emph{fdrtool} package), or proportional density (fixed number of edges) value
#' @param data Can be a dataset or a correlation matrix
#' @param n Number of participants in sample.
#' Defaults to the number of rows in the data.
#' If input is a correlation matrix, then n \strong{must} be set
#' @param normal Should data be transformed to a normal distribution?
#' Defaults to FALSE.
#' Data is not transformed to be normal.
#' Set to TRUE if data should be transformed to be normal
#' (computes correlations using the \emph{cor_auto} function from the \emph{qgraph} package)
#' @param a When thresh = "alpha", "adaptive", and "bonferroni" an alpha threshold is applied (defaults to \strong{.05}).
#' For "adaptive", beta (Type II error) is set to \strong{a*5} for a medium effect size (\emph{r} = \strong{.3}).
#' When thresh = "FDR", a q-value threshold is applied (defaults to \strong{.10}).
#' When thresh = "proportional", a density threshold is applied (defaults to \strong{.15})
#' @param thresh Sets threshold. Defaults to "alpha".
#' Set to any value 0> \emph{r} >1 to retain values greater than set value,
#' "adaptive" for an adapative alpha based on sample size (Perez & Pericchi, 2014),
#' "bonferroni" for the bonferroni correction on alpha,
#' "FDR" for local false discovery rate,
#' and "proportional" for a fixed density of edges (keeps strongest correlations within density)
#' @param na.data How should missing data be handled?
#' For "listwise" deletion the \emph{na.omit} function is applied.
#' Set to "fiml" for Full Information Maxmimum Likelihood (\emph{psych} package).
#' Full Information Maxmimum Likelihood is \strong{recommended} but time consuming
#' @return Returns a list containing a filtered adjacency matrix (A) and the critical r value (r.cv)
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
#' Perez, M. E., & Pericchi, L. R. (2014).
#' Changing statistical significance with the amount of information: The adaptive \emph{a} significance level.
#' \emph{Statistics & Probability Letters}, \emph{85}, 20-24.
#' 
#' Strimmer, K. (2008).
#' fdrtool: A versatile R package for estimating local and tail area-based false discovery rates.
#' \emph{Bioinformatics}, \emph{24}(12), 1461-1462.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Threshold filtering----
threshold <- function (data, a,
                       thresh = c("alpha","adaptive","bonferroni","FDR","proportional"),
                       n = nrow(data), normal = FALSE, na.data = c("pairwise","listwise","fiml","none"))
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
        }else{data<-psych::corFiml(data)}
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
        
        adj.a <- function (a, n)
        {
            alpha<-a*(sqrt(pwr::pwr.r.test(r=.3,power=1-(a*5))$n*(log(pwr::pwr.r.test(r=.3,power=1-(a*5))$n)+qchisq((1-a),1))))/(sqrt(n*(log(n)+qchisq((1-a),1))))
            
            return(alpha)
        }
        
        a<-adj.a(a,nrow(data))
        
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
        
        fdrvec<-fdrtool::fdrtool(fdrvec,plot=FALSE,verbose=FALSE,statistic = "pvalue")$qval
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
#' Betwenness Centrality
#' @description Computes betweenness centrlaity of each node in a network
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted?
#' Defaults to TRUE.
#' Set to FALSE for unweighted measure of betwenness centrality
#' @return A vector of betweenness centrality values for each node in the network
#' @examples
#' A<-TMFG(neoOpen)$A
#' 
#' weighted_BC<-betweenness(A)
#' 
#' unweighted_BC<-betweenness(A,weighted=FALSE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Betweenness Centrality----
betweenness <- function (A, weighted = TRUE)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  if(!weighted)
  {B<-ifelse(A!=0,1,0)
  n<-ncol(B)
  I<-diag(60)
  d<-1
  NPd<-B
  NSPd<-NPd
  NSP<-NSPd
  diag(NSP)<-1
  L<-NSPd
  diag(L)<-1
  while (!is.na(which(NSPd!=0)[1]))
  {
    d<-d+1
    NPd<-as.matrix(NPd)%*%as.matrix(B)
    NSPd<-NPd*(L==0)
    NSP<-NSP+NSPd
    L<-L+d*(NSPd!=0)
  }
  L[!L]<-Inf
  diag(L)<-0
  NSP[!NSP]<-1
  Bt<-t(B)
  DP<-matrix(0,nrow=nrow(B),ncol=ncol(B))
  diam<-d-1
  
  for(d in diam:2)
  {
    DPd1<- (as.matrix(((L==d)*(1+DP)/NSP))%*%as.matrix(Bt))*((L==(d-1))*NSP)
    DP<-DP+DPd1
  }
  BC<-round(as.matrix(colSums(DP),ncol=ncol(A)),0)
  }else{
      G<-ifelse(1/A==Inf,0,1/A)
      
      if(any(G==-Inf))
      {G<-ifelse(G==-Inf,0,G)}
      
      if(any(!G==t(G)))
      {
      if(max(abs(G-t(G)))<1e-10){G<-(G+G)/2}
      }
      
      n<-ncol(G)
      
      BC<-matrix(0,nrow=n,ncol=1)
      
      for(u in 1:n)
      {
          D<-matrix(Inf,nrow=n,ncol=1)
          D[u]<-0
          NP<-matrix(0,nrow=n,ncol=1)
          NP[u]<-1
          S<-matrix(TRUE,nrow=n,ncol=1)
          P<-matrix(FALSE,nrow=n,ncol=n)
          Q<-matrix(0,nrow=n,ncol=1)
          q<-n
          
          G1<-G
          V<-u
          
          while(TRUE)
          {
              S[V]<-0
              G1[,V]<-0
              for(v in V)
              {
                  Q[q]<-v
                  q<-q-1
                  W<-which(G1[v,]!=0)
                  
                  for(w in W)
                  {
                      Duw<-D[v]+G1[v,w]
                      if(Duw<D[w])
                      {
                          D[w]<-Duw
                          NP[w]<-NP[v]
                          P[w,]<-0
                          P[w,v]<-1
                      }else if(Duw==D[w])
                      {
                          NP[w]<-NP[w]+NP[v]
                          P[w,v]<-1
                      }
                  }
              }
              
              
              minD<-suppressWarnings(min(D[S==TRUE]))
              if(length(minD)==0){break}else if(is.infinite(minD))
              {Q[1:q]<-ifelse(length(which(is.infinite(D)))==0,break,which(is.infinite(D)))
              break}
              V<-which(D==minD)
          }
          
          DP<-matrix(0,nrow=n,ncol=1)
          for(w in Q[1:n-1])
          {BC[w]<-BC[w]+DP[w]
          for(v in which(P[w,]!=0))
              DP[v]<-(DP[v]+(1+DP[w]))*NP[v]/NP[w]}
          
      }
  BC<-round(as.matrix(BC,ncol=ncol(A)),0)}
    BC<-as.data.frame(BC)
    row.names(BC)<-colnames(A)
    colnames(BC)<-"BC"
    
  return(BC)
}
#----
#' Randomized Shortest Paths Betweenness Centrality
#' @description Computes betweenness centrlaity based on randomized shortest paths of each node in a network
#' (Please see and cite Kivimaki et al., 2016)
#' @param A An adjacency matrix of network data
#' @param beta Sets the beta parameter.
#' Defaults to 0.01 (recommended).
#' Beta > 0.01 measure gets closer to weighted betweenness centrality (10) and beta < 0.01 measure gets closer to degree (.0001)
#' @return A vector of randomized shortest paths betweenness centrality values for each node in the network
#' @examples
#' A<-TMFG(neoOpen)$A
#' 
#' rspbc<-rspbc(A, beta=0.01)
#' @references 
#' Kivimaki, I., Lebichot, B., Saramaki, J., & Saerens, M. (2016).
#' Two betweenness centrality measures based on Randomized Shortest Paths.
#' \emph{Scientific Reports}, \emph{6}(19668), 1-15.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Randomized Shortest Paths Betweennesss Centrality----
rspbc <- function (A, beta = 0.01)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  n<-ncol(A)
  e<-matrix(1,nrow=n,ncol=1)
  I<-diag(1,nrow=nrow(A),ncol=ncol(A))
  degs<-as.matrix(A)%*%as.matrix(e)
  
  if(any(degs==0))
  {stop("Graph contains unconnected nodes")}
  
  D1<-matrix(0,nrow=nrow(I),ncol=ncol(I))
  for(i in 1:nrow(I))
    for(j in 1:ncol(I))
      if(I[i,j]==1)
      {D1[i,j]<-I[i,j]/degs[i]}
  
  Pref<-as.matrix(D1)%*%as.matrix(A)
  
  bets<-matrix(0,nrow=n,ncol=1)
  C<-1/A
  C<-as.matrix(C)
  C[is.infinite(C)]<-0
  W<-Pref*exp(-(beta)*C)
  rsums<-rowSums(W)
  
  
  Y<-I-W
  Z<-solve(Y,I)
  Zdiv<-1/Z
  Zdiv[Zdiv==Inf]<-0
  DZdiv<-matrix(0,nrow=nrow(Zdiv),ncol=ncol(Zdiv))
  diag(DZdiv)<-diag(Zdiv)
  
  bet<-diag(as.matrix(Z)%*%as.matrix(t(Zdiv-n*DZdiv))%*%as.matrix(Z))
  bet<-round(as.data.frame(bet),0)
  rownames(bet)<-colnames(A)
  bet<-as.matrix(bet)
  return(bet)
}
#----
#' Closeness Centrality
#' @description Computes closeness centrlaity of each node in a network
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted?
#' Defaults to TRUE.
#' Set to FALSE for unweighted measure of closeness centrality
#' @return A vector of closeness centrality values for each node in the network
#' @examples
#' A<-TMFG(neoOpen)$A
#'
#' weighted_LC<-closeness(A)
#' 
#' unweighted_LC<-closeness(A,weighted=FALSE)
#' @references
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export 
#Closeness Centrality----
closeness <- function (A, weighted = TRUE)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  if (!weighted){D<-distance(A,weighted=FALSE)}else if(weighted){D<-distance(A,weighted=TRUE)}
  C<-matrix(0,ncol=ncol(D))
  for(i in 1:ncol(D))
  {
    C[i]<-1/sum(D[,i])
  }
  LC<-t(as.data.frame(C)*100)
  row.names(LC)<-colnames(A)
  colnames(LC)<-"LC"
  LC<-round(as.data.frame(LC),3)
  LC<-as.matrix(LC)
  return(LC)
}
#----
#' Degree
#' @description Computes degree of each node in a network
#' @param A An adjacency matrix of network data
#' @return A vector of degree values for each node in the network.
#' If directed network, returns a list of in-degree (inDegree), out-degree (outDegree), and relative influence (relInf)
#' @examples
#' A<-TMFG(neoOpen)$A
#' 
#' deg<-degree(A)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Degree----
degree <- function (A)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  if(isSymmetric(A)==TRUE)
  {A<-ifelse(A!=0,1,0)
  Deg<-as.data.frame(colSums(A))
  colnames(Deg)<-c("Degree")
  Deg<-as.matrix(Deg)
  return(Deg)
  }else
  {A<-ifelse(A!=0,1,0)
  row.names(A)<-colnames(A)
  inDeg<-colSums(A)
  outDeg<-rowSums(A)
  relinf<-(outDeg-inDeg)/(outDeg+inDeg)
  
    if(all(inDeg==outDeg))
    {Deg<-as.matrix(inDeg)
    colnames(Deg)<-c("Degree")
    return(Deg)
    }else{return(list(inDegree=inDeg,outDegree=outDeg,relInf=relinf))}
  }
}
#----
#' Node Strength
#' @description Computes strength of each node in a network
#' @param A An adjacency matrix of network data
#' @return A vector of strength values for each node in the network.
#' If directed network, returns a list of in-strength (inStrength), out-strength (outStrength), and relative influence (relInf)
#' @examples
#' A<-TMFG(neoOpen)$A
#' 
#' str<-strength(A)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Node Strength----
strength <- function (A)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  if(isSymmetric(A)==TRUE)
  {
  strength<-round(as.data.frame(colSums(A)),2)
  colnames(strength)<-c("Strength")
  strength<-as.matrix(strength)
  return(strength)}else{
      row.names(A)<-colnames(A)
      inStr<-colSums(A)
      outStr<-rowSums(A)
      relinf<-(outStr-inStr)/(outStr+inStr)
      
      if(all(inStr==outStr))
      {Str<-as.matrix(inStr)
      colnames(Str)<-c("Strength")
      return(Str)
      }else{return(list(inStrength=inStr,outStrength=outStr,relInf=relinf))}
  }
}
#----
#' Eigenvector Centrality
#' @description Computes eigenvector centrality of each node in a network
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted?
#' Defaults to TRUE.
#' Set to FALSE for unweighted measure of eigenvector centrality
#' @return A vector of eigenvector centrality values for each node in the network
#' @examples
#' A<-TMFG(neoOpen)$A
#' 
#' weighted_EC<-eigenvector(A)
#' 
#' unweighted_EC<-eigenvector(A,weighted=FALSE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Eigenvector----
eigenvector <- function (A, weighted = TRUE)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  if (!weighted){A<-ifelse(A[]!=0,1,0)
  eigenvector<-eigen(A)
  eigenvector<-round(as.data.frame(abs(eigenvector$vectors[,1])),3)
  colnames(eigenvector)<-c("ECu")}else{eigenvector<-eigen(A)
  eigenvector<-round(as.data.frame(abs(eigenvector$vectors[,1])),3)
  colnames(eigenvector)<-c("ECw")}
  rownames(eigenvector)<-colnames(A)
  eigenvector<-as.matrix(eigenvector)
  return(eigenvector)
}
#----
#' Leverage Centrality
#' @description Computes leverage centrlaity of each node in a network (the degree of connected neighbors)
#' (Please see and cite Joyce et al., 2010)
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted?
#' Defaults to TRUE.
#' Set to FALSE for unweighted measure of leverage centrality
#' @return A vector of leverage centrality values for each node in the network
#' @examples
#' A<-TMFG(neoOpen)$A
#' 
#' weighted_lev<-leverage(A)
#'
#' unweighted_lev<-leverage(A, weighted=FALSE)
#' @references 
#' Joyce, K. E., Laurienti, P. J., Burdette, J. H., & Hayasaka, S. (2010).
#' A new measure of centrality for brain networks. 
#' \emph{PLoS One}, \emph{5}(8), e12200.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export 
#Leverage Centrality----
leverage <- function (A, weighted = TRUE)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  
  if(!weighted)
  {B<-ifelse(A!=0,1,0)}else{B<-A}
  con<-colSums(B)
  lev<-matrix(1,nrow=nrow(B),ncol=1)
  for(i in 1:ncol(B))
  {
    lev[i]<-(1/con[i])*sum((con[i]-con[which(B[,i]!=0)])/(con[i]+con[which(B[,i]!=0)]))
  }
  for(i in 1:nrow(lev))
  if(is.na(lev[i,]))
  {lev[i,]<-0}
  row.names(lev)<-colnames(A)
  return(lev)
}
#----
#' Node Impact
#' @description Computes impact measure or how much the average distance in the network changes with that node removed of each node in a network
#' (Please see and cite Kenett et al., 2011. See also \code{\link[networktools]{impact}} for a separate but related measure)
#' @param A An adjacency matrix of network data
#' @param progBar Defaults to FALSE.
#' Set to TRUE to see progress bar
#' @return A vector of node impact values for each node in the network
#' (impact > 0, greater ASPL when node is removed; impact < 0, lower ASPL when node is removed)
#' @examples
#' A<-TMFG(neoOpen)$A
#' 
#' nodeimpact<-impact(A)
#' @references 
#' Kenett, Y. N., Kenett, D. Y., Ben-Jacob, E., & Faust, M. (2011).
#' Global and local features of semantic networks: Evidence from the Hebrew mental lexicon.
#' \emph{PloS one}, \emph{6}(8), e23912.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Node Impact----
impact <- function (A, progBar = TRUE)
{
    allP<-pathlengths(A)$ASPL
    remove<-matrix(0,nrow=nrow(A),ncol=1)
    
    if(progBar)
    {pb <- txtProgressBar(max=ncol(A), style = 3)}
    
    for(j in 1:ncol(A))
    {
        remove[j,]<-(pathlengths(A[-j,-j])$ASPL)-allP
        
        if(progBar)
        {setTxtProgressBar(pb, j)}
    }
    
    if(progBar)
    {close(pb)}
    
    remove<-round(remove,3)
    colnames(remove)<-"Impact"
    row.names(remove)<-colnames(A)
    return(remove)
}
#----
#' Hybrid Centrality
#' @description Computes hybrid centrality of each node in a network
#' @param A An adjacency matrix of network data
#' @param BC How should the betweenness centrality be computed?
#' Defaults to "standard".
#' Set to "random" for rspbc.
#' Set to "average" for the average of "standard" and "random"
#' @param beta Beta parameter to be passed to the \emph{rspbc} function
#' @return A vector of hybrid centrality values for each node in the network
#' (higher values are more central, lower values are more peripheral)
#' @examples
#' A<-TMFG(neoOpen)$A
#' 
#' HC<-hybrid(A)
#' @references 
#' Pozzi, F., Di Matteo, T., & Aste, T. (2013).
#' Spread of risk across financial markets: Better to invest in the peripheries. 
#' \emph{Scientific Reports}, \emph{3}(1655), 1-7.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Hybrid Centality----
hybrid <- function (A, BC = c("standard","random","average"), beta)
{
  if(missing(BC))
  {BC<-"standard"
  }else{BC<-match.arg(BC)}

  if(missing(beta))
  {beta<-.01
  }else{beta<-beta}
    
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
    if(BC=="standard")
    {
        BCu<-betweenness(A,weighted=FALSE)
        BCw<-betweenness(A)
    }else if(BC=="random")
    {
        BCu<-rspbc(binarize(A),beta=beta)
        BCw<-rspbc(A,beta=beta)
    }else if(BC=="average")
    {
        BCu<-rowMeans(cbind(betweenness(A,weighted=FALSE),rspbc(binarize(A))))
        BCw<-rowMeans(cbind(betweenness(A),rspbc(A)))
    }
  CCu<-closeness(A,weighted=FALSE)
  CCw<-closeness(A)
  if(isSymmetric(A))
  {Deg<-degree(A)
  }else{Deg<-degree(A)$outDeg}
  if(isSymmetric(A))
  {Str<-strength(A)
  }else{Str<-strength(A)$outStr}
  ECu<-eigenvector(A,weighted=FALSE)
  ECw<-eigenvector(A)
  #levu<-leverage(A,weighted=FALSE)
  #levw<-leverage(A)
  #Eu<-PathLengths(A,weighted=FALSE)$ecc
  #Ew<-PathLengths(A)$ecc
  
  hybrid<-((rank(BCu,ties.method="max")+
            rank(BCw,ties.method="max")+
            rank(CCu,ties.method="max")+
            rank(CCw,ties.method="max")+
            rank(Deg,ties.method="max")+
            rank(Str,ties.method="max")+
            rank(ECu,ties.method="max")+
            rank(ECw,ties.method="max")+
            #rank(levu,ties.method="max")+
            #rank(levw,ties.method="max")-
            #rev(rank(Eu,ties.method="max"))+
            #rev(rank(Ew,ties.method="max"))-
            8)/(8*((ncol(A))-1)))
  hybrid<-round(as.data.frame(hybrid),3)
  colnames(hybrid)<-c("HC")
  rownames(hybrid)<-colnames(A)
  hybrid<-as.matrix(hybrid)
  return(hybrid)
}
#----
#' Network Adjusted Mean/Sum
#' @description The hybrid centrality is used to adjust the mean or sum score of participant's factor scores based on each node's centrality.
#' Each participant's response values are multipled by the corresponding hybrid centrality value (uses "average" for BC argument).
#' In this way, more central nodes contribute a greater score and less central nodes contribute a lesser score
#' @param data Must be a dataset
#' @param A Adjacency matrix that has already been filtered
#' @param adjusted Should adjusted values be the mean or sum score?
#' Defaults to "mean".
#' Set to "sum" for sum scores
#' @param factors Can be a vector of factor assignments or community detection algorithms
#' ("walktrap" or "louvain") can be used to determine the number of factors
#' (note that these factors must be at least 9 nodes or more for the TMFG method).
#' Defaults to 1 factor.
#' Set to "walktrap" for the walktrap algortihm.
#' Set to "louvain" for louvain community detection
#' @param standardize Should mean/sum scores be standardized?
#' Defaults to TRUE.
#' Set to FALSE for unstandardized mean/sum scores
#' @param method A network filtering method.
#' Defaults to "TMFG".
#' The "threshold" option is set to the default arguments (thresh = "alpha", a = .05).
#' To use additional arguments for "threshold", apply \emph{threshold} function with
#' desired arguments and set to "none".
#' Use "none" for an adjacency matrix that has already been filtered
#' @param normal Should data be transformed to a normal distribution?
#' Defaults to FALSE.
#' Data is not transformed to be normal.
#' Set to TRUE if data should be transformed to be normal
#' (computes correlations using the \emph{cor_auto} function from the \emph{qgraph} package)
#' @param na.data How should missing data be handled?
#' For "listwise" deletion the \emph{na.omit} function is applied.
#' Set to "fiml" for Full Information Maxmimum Likelihood (\emph{psych} package).
#' Full Information Maxmimum Likelihood is \strong{recommended} but time consuming
#' @param ... Additional arguments for community detection algorithms
#' @return Returns a list of the network adjusted score (NetAdjScore)
#' and the items associated with the specified or identified factors (FactorItems)
#' @examples
#' 
#' sumadj <- nams(neoOpen, adjusted = "sum")
#' 
#' knownfactors <- nams(neoOpen, method = "threshold",
#' factors = c(rep(1,6),rep(2,6),rep(3,6),rep(4,6),rep(5,6),rep(6,6)))
#' 
#' walkadj <- nams(neoOpen, method="threshold", factors = "walktrap")
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Network Adjusted Mean/Sum----
nams <- function (data, A,
                  adjusted = c("mean","sum"),
                  factors = c("walktrap","louvain"), standardize = TRUE,
                  method = c("TMFG","LoGo","PMFG","threshold","none"),
                  normal = FALSE, na.data = c("pairwise","listwise","fiml","none"), ...)
{
    if(missing(data))
    {stop("Dataset is required for analysis")}
    
    if(missing(A))
    {A<-NULL}
    
    if(missing(adjusted))
    {adjusted<-"mean"
    }else(adjusted<-match.arg(adjusted))
    
    if(missing(factors))
    {factors<-rep(1,ncol(data))}
    
    if(missing(method))
    {method<-"TMFG"
    }else{method<-match.arg(method)}
    
    if(missing(na.data))
    {na.data<-"none"
    }else{na.data<-match.arg(na.data)}
    
    
    if(method=="TMFG")
    {net<-TMFG(data,na.data=na.data,normal=normal)$A
    }else if(method=="LoGo")
    {net<-LoGo(data,na.data=na.data,normal=normal,partial=TRUE)$A
    }else if(method=="PMFG")
    {net<-PMFG(data,na.data=na.data,normal=normal)$pmfg
    }else if(method=="threshold")
    {net<-threshold(data,na.data=na.data,normal=normal)$A
    }else if(method=="none")
    {net<-A}
    
    if(!is.numeric(factors))
    {
        if(factors=="walktrap")
        {facts<-igraph::walktrap.community(convert2igraph(net))$membership
        }else if(factors=="louvain")
        {facts<-louvain(net,...)$community}
        
    }else{facts<-factors}
    
    fact<-matrix(0,nrow=nrow(data),ncol=max(facts))
    
    for(i in 1:max(facts))
    {
        
        if(method=="TMFG")
        {Ah<-suppressWarnings(TMFG(net[which(facts==i),which(facts==i)])$A)
        }else if(method=="PMFG")
        {Ah<-PMFG(net[which(facts==i),which(facts==i)])$pmfg
        }else if(method=="LoGo")
        {Ah<-LoGo(net[which(facts==i),which(facts==i)],partial=TRUE)$logo
        }else if(method=="threshold")
        {Ah<-threshold(net[which(facts==i),which(facts==i)])$A
        }else if(method=="none")
        {Ah<-net[which(facts==i),which(facts==i)]}
        
        hyb<-hybrid(Ah,BC="average")
        
        tdata<-t(data[,which(facts==i)])
        
        mat<-matrix(0,nrow=nrow(tdata),ncol=ncol(tdata))
        
        for(j in 1:ncol(tdata))
        {mat[,j]<-hyb*tdata[,j]}
        
        if(adjusted=="mean")
        {wei<-colMeans(mat)
        }else if(adjusted=="sum")
        {wei<-colSums(mat)}
        
        adj<-wei/mean(hyb)
        
        fact[,i]<-adj
    }
    
    if(standardize)
    {fact<-scale(fact)}
    
    fact<-as.data.frame(fact)
    
    for(l in 1:nrow(data))
    {row.names(fact)[l]<-paste("part",l,sep="")}
    
    for(k in 1:max(facts))
        colnames(fact)[k]<-paste("factor",k,sep="")
    
    
    ufacts<-unique(facts)
    
    matf<-matrix(0,nrow=length(ufacts),ncol=1)
    
    for(i in 1:length(ufacts))
    {matf[i,]<-paste(colnames(data)[which(facts==ufacts[i])],collapse = ", ")}
    
    matf<-noquote(matf)
    row.names(matf)<-ufacts
    colnames(matf)<-"Factors"
    
    return(list(NetAdjScore=fact,FactorItems=matf))
}
#----
#' List of Centrality Measures
#' @description Computes centrality measures of the network
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted?
#' Defaults to TRUE.
#' Set to FALSE for unweighted list of centrality measures
#' @return Returns a list of betweenness, closeness, degree (weighted = strength), eigenvector, and leverage centralities
#' @examples
#' A<-TMFG(neoOpen)$A
#' 
#' weighted_centralitylist<-centlist(A)
#' 
#' unweighted_centralitylist<-centlist(A, weighted=FALSE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Centrality List----
centlist <- function (A, weighted = TRUE)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  if(!weighted){BC<-betweenness(A,weighted=FALSE)
  CC<-closeness(A,weighted=FALSE)
  Deg<-degree(A)
  EC<-eigenvector(A,weighted=FALSE)
  lev<-leverage(A,weighted=FALSE)
  hyb<-hybrid(A)
  list(Betweenness=BC,Closeness=CC,Degree=Deg,Eigenvector=EC)}else{
    BC<-betweenness(A)
    CC<-closeness(A)
    Str<-strength(A)
    EC<-eigenvector(A)
    lev<-leverage(A)
    hyb<-hybrid(A)
    return(list(betweenness=BC,closeness=CC,strength=Str,eigenvector=EC,leverage=lev,hybrid=hybrid))}
}
#----
#' Distance
#' @description Computes distance matrix of the network
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted?
#' Defaults to FALSE.
#' Set to TRUE for weighted measure of distance
#' @return A distance matrix of the network
#' @examples
#' A<-TMFG(neoOpen)$A
#' 
#' unweighted_D<-distance(A)
#' 
#' weighted_D<-distance(A, weighted=TRUE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage, \emph{52}(3)}, 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Distance----
distance<-function (A, weighted = FALSE)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  if(!weighted)
  {B<-ifelse(A!=0,1,0)
  l<-1
  Lpath<-B
  D<-B
  Idx<-matrix(TRUE,nrow=nrow(B),ncol=ncol(B))
  while(any(Idx))
  {
    l<-l+1
    Lpath<-(as.matrix(Lpath))%*%(as.matrix(B))
    for(e in 1:nrow(Lpath))
      for(w in 1:ncol(Lpath))
        Idx[e,w]<-(Lpath[e,w]!=0&&(D[e,w]==0))
    D[Idx]<-l
  }
  
  D[!D]<-Inf
  diag(D)<-0
  }else if(weighted){
          G<-ifelse(1/A==Inf,0,1/A)
          
          if(any(G==-Inf))
          {G<-ifelse(G==-Inf,0,G)}
          
          if(any(!G==t(G)))
          {if(max(abs(G-t(G)))<1e-10)
          {G<-(G+G)/2}}
          
          n<-ncol(G)
          D<-matrix(Inf,nrow=n,ncol=n)
          diag(D)<-0
          B<-matrix(0,nrow=n,ncol=n)
          
          for(u in 1:n)
          {
              S<-matrix(TRUE,nrow=n,ncol=1)
              L1<-G
              V<-u
              while(TRUE)
              {
                  S[V]<-0
                  L1[,V]<-0
                  for(v in V)
                  {
                      W<-which(L1[v,]!=0)
                      d<-apply(rbind(D[u,W],(D[u,v]+L1[v,W])),2,min)    
                      wi<-apply(rbind(D[u,W],(D[u,v]+L1[v,W])),2,which.min)
                      D[u,W]<-d
                      ind<-W[wi==2]
                      B[u,ind]<-B[u,v]+1
                  }
                  
                  minD<-suppressWarnings(min(D[u,S==TRUE]))
                  if(length(minD)==0||is.infinite(minD)){break}
                  
                  V<-which(D[u,]==minD)
              }
          }
  }
    colnames(D)<-colnames(A)
    row.names(D)<-colnames(A)
    return(D)
}
#----
#' Characteristic Path Lengths
#' @description Computes global average shortest path length (ASPL),
#' local average shortest path length (ASPLi),
#' eccentricity (ecc),
#' and diameter (D) of a network
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted?
#' Defaults to FALSE.
#' Set to TRUE for weighted measures of ASPL, ASPLi, ecc, and D
#' @return Returns a list of ASPL, ASPLi, ecc, and D of a network
#' @examples
#' A<-TMFG(neoOpen)$A
#' 
#' unweighted_PL<-pathlengths(A)
#' 
#' weighted_PL<-pathlengths(A, weighted=TRUE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Path Lengths----
pathlengths <- function (A, weighted = FALSE)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  if(!weighted)
  {D<-distance(A,weighted=FALSE)}else if(weighted){D<-distance(A,weighted=TRUE)}
  n<-nrow(D)
  for(i in 1:ncol(D))
      for(j in 1:nrow(D))
      if(is.infinite(D[j,i]))
      {D[j,i]<-0}
  if(any(colSums(D)==0))
  {D<-D[,-(which(colSums(D)==0))]}
  aspli<-colSums(D*(D!=0))/(ncol(D)-1)
  aspl<-mean(aspli)
  Emat<-(D*(D!=0))
  ecc<-matrix(nrow=nrow(Emat),ncol=1)
  for(i in 1:nrow(Emat))
  {
    ecc[i,]<-max(Emat[i,])
  }
  d<-max(ecc)
  
  ecc<-as.data.frame(ecc)
  colnames(ecc)<-c("ecc")
  row.names(ecc)<-colnames(A)
  
  return(list(ASPL=aspl,ASPLi=aspli,ecc=ecc,diameter=d))
}
#----
#' Clustering Coefficient
#' @description Computes global clustering coefficient (CC) and local clustering coefficient (CCi)
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted?
#' Defaults to FALSE.
#' Set to TRUE for weighted measures of CC and CCi
#' @return Returns a list of CC and CCi
#' @examples
#' A<-TMFG(neoOpen)$A
#' 
#' unweighted_CC<-clustcoeff(A)
#' 
#' weighted_CC<-clustcoeff(A, weighted=TRUE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Clustering Coefficient----
clustcoeff <- function (A, weighted = FALSE)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  if(!weighted)
  {n<-ncol(A)
  A<-ifelse(A!=0,1,0)
  C<-matrix(0,nrow=n,ncol=1)
  
  for(i in 1:n)
  {
    v<-which(A[i,]!=0)
    k<-length(v)
    if(k >= 2)
    {
      S<-A[v,v]
      C[i]<-sum(S)/(k^2-k)
    }}
  C<-round(t(as.data.frame(C)),3)
  colnames(C)<-colnames(A)
  rownames(C)<-"[,1]"
  CCi<-C
  CC<-mean(C)
  }else{K<-colSums(A!=0)
  m<-A^(1/3)
  cyc<-diag(m%*%m%*%m)
  K[cyc==0]<-Inf
  C<-round(cyc/(K*(K-1)),3)
  CCi<-C
  CC<-mean(C)}
  return(list(CC=CC,CCi=CCi))
}
#----
#' Transitivity
#' @description Computes transitivity of a network
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted?
#' Defaults to FALSE.
#' Set to TRUE for a weighted measure of transitivity
#' @return Returns a value of transitivity
#' @examples
#' A<-TMFG(neoOpen)$A
#' 
#' trans<-transitivity(A, weighted=TRUE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Transitivity----
transitivity <- function (A, weighted = FALSE)
{
    if(!weighted)
    {
        A<-ifelse(A!=0,1,0)
        trans<-sum(diag(A%*%A%*%A))/((sum(A%*%A))-sum(diag(A%*%A)))
    }else if(weighted){
        K<-colSums(ifelse(A!=0,1,0))
        W<-A^(1/3)
        cyc<-diag(W%*%W%*%W)
        trans<-sum(cyc)/sum(K*(K-1))
    }
    
    return(trans)
}
#----
#' Louvain Community Detection Algorithm
#' @description Computes a vector of communities (community) and a global modularity measure (Q)
#' @param A An adjacency matrix of network data
#' @param gamma Defaults to 1.
#' Set to gamma > 1 to detect smaller modules and gamma < 1 for larger modules
#' @param M0 Input can be an initial community vector.
#' Defaults to none
#' @return Returns a list of community and Q
#' @examples
#' A<-TMFG(neoOpen)$A
#' 
#' modularity<-louvain(A)
#' @references
#' Blondel, V. D., Guillaume, J. L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks. 
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}(10), P10008.
#'  
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Louvain Community Detection----
louvain <- function (A, gamma, M0)
{
    if(missing(gamma))
    {gamma<-1
    }else{gamma<-gamma}
    
    if(missing(M0))
    {M0<-1:ncol(A)
    }else(M0<-M0)
    
    n<-ncol(A)
    s<-sum(A)
    
    if(min(A)<0)
    {warning("Matrix contains negative weights: absolute weights were used")
        A<-abs(A)}
    
    Mb<-M0
    M<-M0
    
    mat<-matrix(0,nrow=n,ncol=n)
    for(i in 1:n)
        for(j in 1:n)
        {mat[i,j]<-(colSums(A)[i]*rowSums(A)[j])/s}
    
    B<-A-(gamma*(mat))
    
    B<-(B+t(B))/2
    
    Hnm<-matrix(0,nrow=nrow(A),ncol=(length(unique(Mb))))
    

    for(m in 1:max(Mb))
    {
        if(!is.null(nrow(B[,which(Mb==m)])))
        {Hnm[,m]<-rowSums(B[,which(Mb==m)])
        }else{Hnm[,m]<-B[,which(Mb==m)]}
    }
    
    H<-colSums(Hnm)
    Hm<-rowSums(Hnm)
    
    Q0<-(-Inf)
    bsxfun<-matrix(0,nrow=n,ncol=n)
    diag(bsxfun)<-1
    Q<-sum(diag(as.matrix(B)*bsxfun))/s
    
    
    first_iter<-TRUE
    while(Q-Q0>0)
    {
        flag<-TRUE
        while(flag)
        {
            set.seed(0)
            flag<-FALSE
            for(u in sample(n))
            {
                ma<-Mb[u]
                dQ<-Hnm[u,] - Hnm[u,ma] + B[u,u]
                dQ[ma]<-0
                
                max_dQ<-max(dQ)
                mb<-which.max(dQ)
                
                if(max_dQ>0)
                {
                    flag<-TRUE
                    Mb[u]<-mb
                    
                    Hnm[,mb]<-Hnm[,mb]+B[,u]
                    Hnm[,ma]<-Hnm[,ma]-B[,u]
                    Hm[mb]<-Hm[mb]+H[u]
                    Hm[ma]<-Hm[ma]-H[u]
                }
            }
        }
        Mb<-match(Mb,unique(Mb))
        
        M0<-M
        if(first_iter)
        {
            M<-Mb
            first_iter<-FALSE
        }else{
            for(u in 1:n)
            {
                M[M0==u]<-Mb[u]
            }
        }
        
        n<-max(Mb)
        B1<-matrix(0,nrow=n,ncol=n)
        for(u in 1:n)
            for(v in u:n)
            {
                bm<-sum(sum(B[Mb==u,Mb==v]))
                B1[u,v]<-bm
                B1[v,u]<-bm
            }
        B<-B1
        
        Mb<-1:n
        Hnm<-B
        H<-colSums(B)
        
        Q0<-Q
        
        Q<-sum(diag(B))/s
        
    }
    return(list(community=M,Q=Q))
}
#----
#' Small-worldness Measure
#' @description Computes the small-worldness measure of a network
#' @param A An adjacency matrix of network data
#' @param iter Number of random (or lattice) networks to generate,
#' which are used to calculate the mean random ASPL and CC (or lattice)
#' @param progBar Defaults to FALSE.
#' Set to TRUE to see progress bar
#' @param method Defaults to "HG" (Humphries & Gurney, 2008).
#' Set to "rand" for the CC to be calculated using a random network or
#' set to "TJHBL" for (Telesford et al., 2011) where CC is calculated from a lattice network
#' @return Returns a list with value of small-worldness (swm) and random ASPL (rASPL).
#' For "rand", values > 1 indicate a small-world network (lrCCt = random CC).
#' For "HG", values > 3 indicate a small-world network (lrCCt = random Transitivity).
#' For "TJHBL" values near 0 indicate a small-world network (lrCCt = lattice CC),
#' while < 0 indicates a more regular network and > 0 indicates a more random network
#' @examples
#' A<-TMFG(neoOpen)$A
#'
#' swmHG <- smallworldness(A, method="HG")
#' 
#' swmRand <- smallworldness(A, method="rand")
#' 
#' swmTJHBL <- smallworldness(A, method="TJHBL")
#' @references 
#' Humphries, M. D., & Gurney, K. (2008).
#' Network 'small-world-ness': A quantitative method for determining canonical network equivalence.
#' \emph{PloS one}, \emph{3}(4), e0002051.
#' 
#' Telesford, Q. K., Joyce, K. E., Hayasaka, S., Burdette, J. H., & Laurienti, P. J. (2011).
#' The ubiquity of small-world networks.
#' \emph{Brain Connectivity}, \emph{1}(5), 367-375.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Small-worldness Measure----
smallworldness <- function (A, iter = 100, progBar = FALSE, method = c("HG","rand","TJHBL"))
{
    if(missing(method))
    {method<-"HG"
    }else{method<-match.arg(method)}
    
    mat<-matrix(0,nrow=nrow(A),ncol=ncol(A)) #Initialize bootstrap matrix
    asamps<-matrix(0,nrow=iter) #Initialize sample matrix
    csamps<-matrix(0,nrow=iter) #Initialize sample matrix
    if(progBar)
    {pb <- txtProgressBar(max=iter, style = 3)}
    for(i in 1:iter) #Generate array of bootstrapped samples
    {
        f<-round(runif(i,min=1,max=1000000),0)
        set.seed(f[round(runif(i,min=1,max=length(f)),0)])
        rand<-randnet(ncol(A),sum(ifelse(A!=0,1,0))/2)
        if(method=="TJHBL")
        {latt<-lattnet(ncol(A),sum(ifelse(A!=0,1,0))/2)}
        asamps[i,]<-pathlengths(rand)$ASPL
        if(method=="rand")
        {csamps[i,]<-clustcoeff(rand)$CC
        }else if(method=="HG"){csamps[i,]<-transitivity(rand)
        }else if(method=="TJHBL"){csamps[i,]<-clustcoeff(latt)$CC}else{stop("Method not available")}
        if(progBar)
        {setTxtProgressBar(pb, i)}
    }
    if(progBar)
    {close(pb)}
    
    nodes<-ncol(A)
    ASPL<-pathlengths(A)$ASPL
    CC<-clustcoeff(A)$CC
    trans<-transitivity(A)
    rASPL<-mean(asamps)
    
    if(method=="rand")
    {rCC<-mean(csamps)
    swm<-(CC/rCC)/(ASPL/rASPL)
    lrCCt<-rCC
    }else if(method=="HG")
    {rtrans<-mean(csamps)
    swm<-(trans/rtrans)/(ASPL/rASPL)
    lrCCt<-rtrans
    }else if(method=="TJHBL")
    {lCC<-mean(csamps)
    swm<-(rASPL/ASPL)-(CC/lCC)
    lrCCt<-lCC
    }
    
    return(list(swm=swm, rASPL=rASPL, lrCCt=lrCCt))
}
#----
#' Semantic Network Measures
#' @description Computes the average shortest path length (ASPL),
#' clustering coefficient(CC),
#' modularity (Q),
#' and small-worldness (S; defaults to "rand") 
#' @param A An adjacency matrix of network A
#' @param iter Number of iterations for the small-worldness measure
#' @return Returns a values for ASPL, CC, Q, and S
#' @examples
#' \dontrun{
#' lowO <- subset(animals, Group==1)[-1]
#' 
#' A<-PMFG(lowO)
#' 
#' globmeas<-semnetmeas(A)
#' }
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Semantic Network Measures----
semnetmeas <- function (A, iter)
{
    if(missing(iter))
    {iter<-100
    }else{iter<-iter}
    
    aspl<-pathlengths(A, weighted = FALSE)$ASPL
    cc<-clustcoeff(A, weighted = FALSE)$CC
    q<-louvain(A)$Q
    s<-smallworldness(A,iter=iter,method="rand",progBar = FALSE)
    swm<-s$swm
    raspl<-s$rASPL
    rCC<-s$lrCCt
    
    semnetmeas<-cbind(aspl,cc,q,swm,raspl,rCC)
    
    semnetmeas<-as.data.frame(semnetmeas)
    
    colnames(semnetmeas)<-c("ASPL","CC","Q","S","randASPL","randCC")
    
    semnetmeas<-as.matrix(semnetmeas)
    
    return(semnetmeas)
}
#----
#' Partial Bootstrapped Semantic Network Analysis
#' @description Bootstraps (without replacement) the nodes in the network and computes global network characteristics
#' @param data A set of data
#' @param method A network filtering method.
#' Defaults to "PMFG"
#' @param nodes Number of nodes (i.e., variables) to use in the bootstrap.
#' Defaults to 50% of nodes in the network.
#' Otherwise accepts the number of the nodes to be included
#' @param iter Number of bootstrap iterations.
#' Defaults to 1000 iterations
#' @param normal Should data be transformed to a normal distribution?
#' Defaults to FALSE.
#' Data is not transformed to be normal.
#' Set to TRUE if data should be transformed to be normal
#' (computes correlations using the \emph{cor_auto} function from the \emph{qgraph} package)
#' @param na.data How should missing data be handled?
#' For "listwise" deletion the \emph{na.omit} function is applied.
#' Set to "fiml" for Full Information Maxmimum Likelihood (\emph{psych} package).
#' Full Information Maxmimum Likelihood is \strong{recommended} but time consuming
#' @param seeds Seeds to use for random number generation.
#' Defaults to NULL.
#' Input seeds from previous run (\strong{see examples})
#' @param ... Additional arguments for filtering methods
#' @return Returns a list that includes the original semantic network measures (origmeas; ASPL, CC, Q, S),
#' the bootstrapped semantic network measures (bootmeas),
#' and the seeds used in the random number generator (Seeds)
#' @examples 
#' \dontrun{
#' lowO <- subset(animals, Group==1)[-1]
#' 
#' semPMFG<-semnetboot(lowO)
#' 
#' #Replicates previous analysis
#' repsemPMFG<-semnetboot(lowO, seeds = semPMFG$Seeds)
#' 
#' semLoGo<-semnetboot(lowO,method="LoGo")
#' 
#' semMaST<-semnetboot(lowO,method="MaST")
#' 
#' semThreshold<-semnetboot(lowO,method="threshold")
#' }
#' @references
#' Kenett, Y. N., Wechsler-Kashi, D., Kenett, D. Y., Schwartz, R. G., Ben Jacob, E., & Faust, M. (2013).
#' Semantic organization in children with cochlear implants: Computational analysis of verbal fluency.
#' \emph{Frontiers in Psychology}, {4}(543), 1-11.
#' @export
#Partial Bootstrapped Semantic Network Analysis----
semnetboot <- function (data, method = c("PMFG","TMFG","LoGo","MaST","threshold"),
                        normal = FALSE, nodes,
                        iter = 1000, na.data = c("pairwise","listwise","fiml","none"),
                        seeds = NULL, ...)
{
    if(missing(method))
    {method<-"PMFG"
    }else{method<-match.arg(method)}
    
    if(missing(nodes))
    {nodes<-round((ncol(data)/2),0)
    }else(nodes<-round(nodes,0))
    
    if(!is.null(seeds))
    {
        if(length(seeds)!=iter)
        {iter<-length(seeds)}
    }
    
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
        }else{data<-psych::corFiml(data)}
    }
    
    #corrlation matrix
    if(nrow(data)==ncol(data)){cormat<-data
    }else if(normal){cormat<-qgraph::cor_auto(data)
    }else{cormat<-cor(data)}
    
    mat<-matrix(0,nrow=nrow(data),ncol=nodes) #Initialize bootstrap matrix
    samps<-matrix(0,nrow=iter,ncol=4) #Initialize sample matrix
    Seeds<-vector(mode="numeric",length=iter) #Initialize seeds vector
    
    pb <- txtProgressBar(max=iter, style = 3)
    for(i in 1:iter) #Generate array of bootstrapped samples
    {
        if(is.null(seeds))
        {f<-round(runif(i,min=1,max=1000000),0)
        Seeds[i]<-f[round(runif(i,min=1,max=length(f)),0)][i]
        }else{Seeds<-as.vector(seeds)}
        set.seed(Seeds[i])
        mat<-data[sample(1:ncol(data),replace=TRUE)[1:nodes],]
        
        #missing data handling
        if(any(is.na(mat)))
        {
            if(na.data=="pairwise")
            {
                if(normal)
                {cormat<-qgraph::cor_auto(mat,missing=na.data)
                }else{cormat<-psych::cor2(mat,use=na.data)}
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
        }
        
        if(nrow(mat)==ncol(mat)){cormat<-mat
        }else if(normal){cormat<-qgraph::cor_auto(mat)
        }else{cormat<-cor(mat)}
        
        if(method=="TMFG")
        {samps[i,]<-suppressWarnings(semnetmeas(binarize(TMFG(cormat)$A),iter=10)[1:4])
        }else if(method=="LoGo")
        {l<-LoGo(mat,partial=TRUE)$logo
        samps[i,]<-suppressWarnings(semnetmeas(binarize(l),iter=10)[1:4])
        }else if(method=="MaST")
        {samps[i,]<-suppressWarnings(semnetmeas(binarize(MaST(cormat,...)),iter=10)[1:4])
        }else if(method=="PMFG")
        {samps[i,]<-suppressWarnings(semnetmeas(binarize(PMFG(cormat,progBar=FALSE)$pmfg),iter=10)[1:4])
        }else if(method=="threshold")
        {samps[i,]<-suppressWarnings(semnetmeas(binarize(threshold(cormat,...)$A),iter=10)[1:4])
        }else {stop("Method not available")}
        
        setTxtProgressBar(pb, i)
    }
    close(pb)
    
    if(method=="TMFG")
    {tru<-suppressWarnings(semnetmeas(TMFG(data,...)$A))
    }else if(method=="LoGo")
    {tru<-suppressWarnings(semnetmeas((-cov2cor(LoGo(data)$logo))))
    }else if(method=="MaST")
    {tru<-suppressWarnings(semnetmeas(MaST(data,...)))
    }else if(method=="PMFG")
    {tru<-suppressWarnings(semnetmeas(PMFG(data)))
    }else if(method=="MaST")
    {tru<-suppressWarnings(semnetmeas(MaST(data,...)))
    }else if(method=="ECO")
    {tru<-suppressWarnings(semnetmeas(threshold(data,...)$A))
    }else stop("Method not available")
    
    colnames(samps)<-c("ASPL","CC","Q","S")
    
    return(list(origmeas=tru,bootmeas=samps,Seeds=Seeds))
}
#----
#' Edge Replication
#' @description Computes the number of edges that replicate between two cross-sectional networks
#' @param A An adjacency matrix of network A
#' @param B An adjacency matrix of network B
#' @param corr Correlation method for assessing the relatonship between the replicated edge weights.
#' Defaults to "pearson".
#' Set to "spearman" for non-linear or monotonic associations.
#' Set to "kendall" for rank-order correlations
#' @param plot Should there be a plot of the replicated weights?
#' Defaults to FALSE.
#' Set to TRUE for a plot to be produced
#' @return Returns a list of the edges that replicated and their weights (replicatedEdges),
#' number of edges that replicate (replicated),
#' total number of edges (totalEdgesA & totalEdgesB),
#' the percentage of edges that replicate (percentageA & percentageB),
#' the density of edges (densityA & densityB),
#' the mean difference between edges that replicate (meanDifference),
#' the sd of the difference between edges that replicate (sdDifference),
#' and the correlation between the edges that replicate for both networks (correlation)
#' @examples
#' A<-TMFG(neoOpen)$A
#' 
#' B<-MaST(neoOpen)
#' 
#' edges<-edgerep(A,B)
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Edge Replication----
edgerep <- function (A, B, corr = c("pearson", "spearman", "kendall"), plot = FALSE)
{
    if(missing(corr))
    {corr<-"pearson"
    }else{corr<-match.arg(corr)}
    
if(!isSymmetric(A))
{
    if(all(rowSums(A)==colSums(A)))
    {A<-as.matrix(Matrix::forceSymmetric(A))
    }else{A<-A+t(A)
    warning("Adjacency matrix A was made to be symmetric")}
}
 
if(!isSymmetric(B))
{
    if(all(rowSums(B)==colSums(B)))
    {B<-as.matrix(Matrix::forceSymmetric(B))
    }else{B<-B+t(B)
    warning("Adjacency matrix B was made to be symmetric")}
}
           
    n<-ncol(A)
    
    diag(A)<-0
    diag(B)<-0
  
  repedges<-which(A&B!=0)
  repmat<-matrix(0,nrow=n,ncol=n)
  repmat[repedges]<-1
  count<-sum(repmat)/2
  
  repmat[lower.tri(repmat)]<-A[lower.tri(A)]
  repmat[upper.tri(repmat)]<-B[upper.tri(B)]
  
  i<-sort(c(rep(1:n,n)))
  j<-c(rep(1:n,n))
  lt<-as.vector(A)
  ut<-as.vector(B)
  
  replist<-cbind(i,j,lt,ut)
  replist<-replist[which(replist[,3]&replist[,4]!=0),]
  ord<-replist[order(replist[,3]),]
  sinlist<-ord[seq(1,nrow(ord),2),]
  sinlist<-sinlist[order(sinlist[,1],sinlist[,2]),]
  
  colnames(sinlist)<-c("nodeTo","nodeFrom","weight in A","weight in B")
  for(i in 1:nrow(sinlist))
      for(j in 1:length(colnames(A)))
      {
          if(sinlist[i,1]==j)
          {sinlist[i,1]<-colnames(A)[j]}
          if(sinlist[i,2]==j)
          {sinlist[i,2]<-colnames(A)[j]}
      }
  sinlist<-noquote(sinlist)
  sinlist[,c(3,4)]<-round(as.numeric(sinlist[,c(3,4)]),3)
  
  
  colnames(repmat)<-colnames(A)
  row.names(repmat)<-colnames(A)
  
  if(plot)
  {
    if(ncol(A)<=20)
    {corrplot::corrplot(repmat,method="color",addgrid.col = "grey",
                      tl.col="black",addCoef.col = "black",
                      title="A by B Edge Replcation",mar=c(3,4,2,2))
    }else{corrplot::corrplot(repmat,method="color",addgrid.col = "grey",
                           tl.col="black",
                           title="A by B Edge Replcation",mar=c(3,4,2,2))}
  }
                
                
  diag(A)<-0
  diag(B)<-0
  possibleA<-sum(ifelse(A!=0,1,0)/2)
  possibleB<-sum(ifelse(B!=0,1,0)/2)
  percentA<-count/possibleA
  percentB<-count/possibleB
  densityA<-possibleA/((ncol(A)^2-ncol(A))/2)
  densityB<-possibleB/((ncol(B)^2-ncol(B))/2)
  
  mat<-matrix(0,nrow=nrow(A),ncol=ncol(A))
  wc<-0
  wveca<-0
  wvecb<-0
  for(i in 1:ncol(A))
      for(j in 1:nrow(A))
          if(A[i,j]&&B[i,j]!=0)
          {
              mat[i,j]<-abs(A[i,j]-B[i,j])
              wc<-wc+1
              wveca[wc]<-A[i,j]
              wvecb[wc]<-B[i,j]
          }
  
  corr<-cor(wveca,wvecb,method=corr)
  
  m<-0
  vec<-0
  for(i in 1:nrow(A))
      for(j in 1:ncol(A))
          if(mat[i,j]!=0)
          {m<-m+1
              vec[m]<-mat[i,j]
              mvec<-mean(vec)
              svec<-sd(vec)}else if(all(mat==0)){mvec<-0
                                                 svec<-0}
  
  return(list(replicatedEdges=sinlist,
              replicated=count,
              totalEdgesA=possibleA,totalEdgesB=possibleB,
              percentageA=percentA,percentageB=percentB,
              densityA=densityA,densityB=densityB,
              meanDifference=mvec,sdDifference=svec,correlation=corr))
}
#----
#' Network Connectivity
#' @description Computes the average and standard deviation of the weights in the network
#' @param A An adjacency matrix of network A
#' @return Returns a list of the edge weights (weights),
#' the mean (mean),
#' the standard deviation (sd),
#' and the sum of the edge weights (total) in the network
#' @examples
#' A<-TMFG(neoOpen)$A
#' 
#' connectivity<-conn(A)
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Network Connectivity----
conn <- function (A)
{
    weights<-0
    wc<-0
    B<-A[lower.tri(A)]
    for(i in 1:length(B))
        if (B[i]!=0)
        {
            wc <- wc+1
            weights[wc] <- B[i]
        }
    tot<-sum(weights)
    mea<-mean(weights)
    s<-sd(weights)
    
    return(list(weights=weights,mean=mea,sd=s,total=tot))
}
#----
#' Bootstrapped Network Generalization
#' @description Bootstraps the sample to identify the most stable correlations.
#' Also produces a network that is penalizes low reliability edges.
#' This function is useful for overcoming the structural constraint of the IFN approach
#' @param data A set of data
#' @param method A network filtering method.
#' Defaults to "TMFG"
#' @param alpha Alpha threshold bootstrapped network generalization network
#' @param n Number of people to use in the bootstrap.
#' Defaults to full sample size
#' @param iter Number of bootstrap iterations.
#' Defaults to 100 iterations
#' @param normal Should data be transformed to a normal distribution?
#' Defaults to FALSE.
#' Data is not transformed to be normal.
#' Set to TRUE if data should be transformed to be normal
#' (computes correlations using the \emph{cor_auto} function from the \emph{qgraph} package)
#' @param na.data How should missing data be handled?
#' For "listwise" deletion the \emph{na.omit} function is applied.
#' Set to "fiml" for Full Information Maxmimum Likelihood (\emph{psych} package).
#' Full Information Maxmimum Likelihood is \strong{recommended} but time consuming
#' @param seeds Seeds to use for random number generation.
#' Defaults to NULL.
#' Input seeds from previous run (\strong{see examples})
#' @param cores Number of computer processing cores to use for bootstrapping samples.
#' Defaults to \emph{n} - 1 total number of cores.
#' Set to any number between 1 and maxmimum amount of cores on your computer
#' @param progBar Should progress bar be displayed?
#' Defaults to TRUE.
#' Set to FALSE for no progress bar
#' @param ... Additional arguments for filtering methods
#' @return Returns a list that includes the original filtered network (orignet),
#' correlation matrix of the mean bootstrapped network (bootmat),
#' unfiltered reliabilities of all of the connections (netrel),
#' seeds used in the random number generator (Seeds)
#' @examples
#' \dontrun{
#' bootTMFG<-bootgen(neoOpen)
#' 
#' #Replicates previous analysis
#' repbootTMFG<-bootgen(neoOpen, seeds = bootTMFG$Seeds)
#' 
#' bootLoGo<-bootgen(neoOpen,method="LoGo")
#' 
#' bootMaST<-bootgen(neoOpen,method="MaST")
#' 
#' bootThreshold<-bootgen(neoOpen,method="threshold")
#' }
#' @references
#' Perez, M. E., & Pericchi, L. R. (2014).
#' Changing statistical significance with the amount of information: The adaptive \emph{a} significance level.
#' \emph{Statistics & Probability Letters}, \emph{85}, 20-24.
#' 
#' Tumminello, M., Coronnello, C., Lillo, F., Micciche, S., & Mantegna, R. N. (2007).
#' Spanning trees and bootstrap reliability estimation in correlation-based networks.
#' \emph{International Journal of Bifurcation and Chaos}, \emph{17}(7), 2319-2329.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @importFrom stats cov2cor pt
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#Bootstrap Network Generalization----
bootgen <- function (data, method = c("MaST", "PMFG", "TMFG", "LoGo", "threshold"),
                     alpha = .05, n = nrow(data), iter = 100, normal = FALSE,
                     na.data = c("pairwise", "listwise", "fiml","none"),
                     seeds = NULL, cores, progBar = TRUE, ...)
{
    
    #arguments
    ##########################################################
    if(missing(method))
    {method<-"TMFG"
    }else{method<-match.arg(method)}
    
    if(!is.null(seeds))
    {
        if(length(seeds)!=iter)
        {iter<-length(seeds)}
    }
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
            {
                z[i,j]<--z[i,j]
            }
        }
    return(z)}
    
    #mean fisher's z
    zw <- function (z,iter)
    {sum((iter-3)*(z))/((iter-3)*iter)}
    
    #critical r value
    critical.r <- function(nrow, a){
        df <- nrow - 2
        critical.t <- qt( a/2, df, lower.tail = F )
        cvr <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
        return(cvr)}
    
    #adaptive alpha
    adapt <- function (n, a){
        alpha<-a*(sqrt(pwr::pwr.r.test(r=.3,power=1-(a*5))$n*(log(pwr::pwr.r.test(r=.3,power=1-(a*5))$n)+qchisq((1-a),1))))/(sqrt(n*(log(n)+qchisq((1-a),1))))
        return(alpha)}
    
    #pcor sig
    psig <- function (pcor, n, gp=(ncol(data)-2))
    {
        statistic <- pcor*sqrt((n-2-gp)/(1-pcor^2))
        p.value <- 2*pt(-abs(statistic),(n-2-gp))
        return(p.value)
    }
    
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
        }else{data<-psych::corFiml(data)}
    }
    ##########################################################
    
    #corrlation matrix
    if(nrow(data)==ncol(data)){cormat<-data
    }else if(normal){cormat<-qgraph::cor_auto(data)
    }else{cormat<-cor(data)}
    
    realmat<-cormat
    
    #general bootstrapping    
    ##########################################################
    sampslist<-list() #initialize sample list
    
    #if(method=="LoGo") #partial correlation sample list
    #{tsampslist<-list()}
    
    #Seeds
    Seeds<-vector(mode="numeric",length=iter) #initialize seeds vector

    for(i in 1:iter) #Generate array of bootstrapped samples
    {
        if(is.null(seeds))
        {f<-round(runif(i,min=1,max=1000000),0)
        Seeds[i]<-f[round(runif(i,min=1,max=length(f)),0)][i]
        }else{Seeds<-as.vector(seeds)}  #input replication seeds
    }
    
    #Parallel processing
    if(missing(cores))
    {
        cores <- parallel::detectCores()
        cores <- cores - 1 #not to overload your computer 
    }else{cores <- cores}
    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    
    if(progBar)
    {
        pb <- txtProgressBar(max=iter, style = 3) #progress bar
        progress <- function(num) setTxtProgressBar(pb, num)
        opts <- list(progress = progress)
    }else{opts<-list()}
    
    sampslist<-foreach::foreach(i=1:length(Seeds),
                         .packages = "NetworkToolbox",
                         .options.snow = opts)%dopar%
    {
        set.seed(Seeds[i])
        mat<-data[sample(1:n,replace=TRUE),]
        
        #missing data handling
        if(any(is.na(mat)))
        {
            if(na.data=="pairwise")
            {
                if(normal)
                {cormat<-qgraph::cor_auto(mat,missing=na.data)
                }else{cormat<-psych::cor2(mat,use=na.data)}
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
        }
        
        #correlation matrix
        if(nrow(mat)==ncol(mat)){cormat<-mat
        }else if(normal){cormat<-qgraph::cor_auto(mat)
        }else{cormat<-cor(mat)}
        
        if(method=="TMFG")
        {samps<-fish(TMFG(cormat)$A)
        }else if(method=="LoGo")
        {
            samps<-fish(suppressWarnings(LoGo(mat,partial=TRUE)$logo))
            #tsamps<-fish(TMFG(cormat)$A)
        }else if(method=="MaST")
        {samps<-fish(MaST(cormat,...))
        }else if(method=="PMFG")
        {samps<-fish(PMFG(cormat)$pmfg)
        }else if(method=="threshold")
        {samps<-fish(threshold(cormat,...)$A)
        }else{stop("Method not available")}
        
        #if(method=="LoGo")
        #{return(list(samps=samps,tsamps=tsamps))
        #}else{
            return(list(samps=samps))
            #}
        
    }
    if(progBar)
    {close(pb)}
    parallel::stopCluster(cl)
    ##########################################################
    
    #Original Networks
        if(method=="TMFG")
        {tru<-TMFG(data)$A
        }else if(method=="LoGo")
        {tru<-LoGo(data,partial=TRUE)$logo
        diag(tru)<-0
        }else if(method=="MaST")
        {tru<-MaST(data,...)
        }else if(method=="PMFG")
        {tru<-PMFG(data)$pmfg
        }else if(method=="MaST")
        {tru<-MaST(data,...)
        }else if(method=="threshold")
        {tru<-threshold(data,...)$A
        }else stop("Method not available")
    
    diag(tru)<-0
    
    #convert samples list from foreach to array
    samps<-array(0,dim=c(nrow=nrow(realmat),ncol=ncol(realmat),iter))
    
    for(i in 1:iter) #populate array
    {samps[,,i]<-sampslist[[i]]$samps}
    
    #convert partial samples list from foreach
    #if(method=="LoGo")
    #{
    #    tsamps<-array(0,dim=c(nrow=nrow(realmat),ncol=ncol(realmat),iter))
    #    
    #    for(i in 1:iter) #populate array
    #    {tsamps[,,i]<-sampslist[[i]]$tsamps}
    #}
    
    #Mean fisher matrix
    meanmat<-matrix(0,nrow=nrow(realmat),ncol=ncol(realmat)) #Initialize Mean matrix
    for(j in 1:nrow(realmat))
        for(k in 1:ncol(realmat))
        {meanmat[j,k]<-zw(samps[j,k,],iter)}
    
    #partial fisher mean matrix
    #if(method=="LoGo")
    #{
    #tmeanmat<-matrix(0,nrow=nrow(realmat),ncol=ncol(realmat)) #Initialize Mean matrix
    #for(j in 1:nrow(realmat))
    #    for(k in 1:ncol(realmat))
    #    {tmeanmat[j,k]<-zw(tsamps[j,k,],iter)}
    #}
    
    #convert fisher z to r
    meanmat<-psych::fisherz2r(meanmat)
    #if(method=="LoGo")
    #{tmeanmat<-psych::fisherz2r(tmeanmat)}

    #threshold value
        cvr<-critical.r(n,adapt(n,alpha))
        
        if(method=="LoGo")
        {
           sigmat<-meanmat
            
            a<-adapt(n,alpha)
            
            for(i in 1:nrow(meanmat))
                for(j in 1:ncol(meanmat))
                {sigmat[i,j]<-psig(meanmat[i,j],n)}
            
            #tmeanmat<-ifelse(abs(tmeanmat)>=cvr,tmeanmat,0)
            #meanmat<-ifelse(tmeanmat!=0,meanmat,0) #use correlation critical value
            meanmat<-ifelse(sigmat<=a,meanmat,0) #use partial critical value
        }else{meanmat<-ifelse(abs(meanmat)>=cvr,meanmat,0)}
    
    #return meanmat to bootmat
    #bootmat<-ifelse(meanmat!=0,realmat,0)
    bootmat<-meanmat
    colnames(bootmat)<-colnames(realmat)
    
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
    
    diag(bootmat)<-0
    
    if(!isSymmetric(bootmat))
    {bootmat<-as.matrix(Matrix::forceSymmetric(bootmat))}
    
    if(method=="LoGo")
    {is.graphical(bootmat)}
    
    return(list(orignet=tru,bootmat=bootmat,netrel=rel,Seeds=Seeds))
}
#----
#' Bootstrapped Network Generalization Plots
#' @description Generates reliability plots from the \emph{bootgen} function
#' @param object An output from the \emph{bootgen} function
#' @param bootmat Should the bootstrap generalization matrix be included?
#' Defaults to FALSE.
#' Set to TRUE to plot the original and bootmat networks' reliablities
#' @return Returns the following plots:
#' a reliability matrix for the original network
#' ("Original Network Correlation Reliabilities"
#' upper triangle = network reliabilites from the filtered network,
#' lower triangle = full, unfiltered matrix reliablities),
#' a plot of retained correlations on their reliability ("Original Network Correlation Strength on Reliability"),
#' If bootmat is set to TRUE, then also 
#' a reliability matrix for the bootgen network
#' ("bootgen Network Correlation Reliabilities"
#' upper triangle = network reliabilites from the bootgen filtered network,
#' lower triangle = full, unfiltered matrix reliablities),
#' a plot of retained correlations on their reliability ("bootgen Network Correlation Strength on Reliability"),
#' a frequency of the reliability of edges in both the original and bootgen network
#' @examples
#' \dontrun{
#' bootTMFG <- bootgen(neoOpen)
#' 
#' bootPlot <- bootPlot(bootTMFG)
#' }
#' @references
#' Wei, T. & Simko, V.(2017).
#' R package "corrplot": Visualization of a correlation matrix (Version 0.84).
#' Available from \url{https://github.com/taiyun/corrplot}
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @importFrom graphics abline plot text hist legend
#' @importFrom grDevices rgb
#' @importFrom stats lm na.omit
#' @export
#Bootstrap Network Generalization Plots----
bootgenPlot <- function (object, bootmat = FALSE)
{
    tru<-object$orignet
    rel<-object$netrel
    if(bootmat==TRUE)
    {
        boot<-object$bootmat
        brel<-rel
    }
    
    #create tru reliability plot matrix
    row.names(rel)<-colnames(rel)
    diag(rel)<-1
    upp<-matrix(0,nrow=nrow(rel),ncol=ncol(rel))
    for(i in 1:nrow(rel))
        for(j in 1:ncol(rel))
            if(rel[i,j]!=0&&tru[i,j]!=0)
            {upp[i,j]<-rel[i,j]}
    colnames(upp)<-colnames(rel)
    rel[upper.tri(rel)]<-upp[upper.tri(upp)]
    
    #reliablity on correlation plot
    x<-matrix(nrow=length(upp))
    y<-matrix(nrow=length(tru))
    diag(upp)<-0
    diag(tru)<-0
    wc<-0
    for(i in 1:nrow(tru))
        for(j in 1:ncol(tru))
            if((upp[i,j]!=0&&tru[i,j])!=0)
            {wc<-wc+1
            x[wc]<-upp[i,j]
            y[wc]<-tru[i,j]}
    xo<-na.omit(abs(x))
    yo<-na.omit(abs(y))
    
    dev.new()
    mar=c(2,2,2,2)
    cpo<-{plot(xo,yo,pch=16,ylab="Edge Weight Strength",xlab="Reliability",
               main="Original Network Edge Weight Strength on Reliability",xlim=c(0,1),ylim=range(yo))
        abline(lm(yo~xo))
        text(x=.05,y=max(yo-.05),labels = paste("r = ",round(cor(yo,xo),3)))}
    
    
    #plot reliability matrix
    dev.new()
    if(ncol(tru)<=20)
    {plt<-corrplot::corrplot(rel,method="color",
                             title="Original Network Edge Weight Reliabilities",
                             mar=c(2,2,2,2),tl.col="black",tl.cex=.75,
                             cl.lim = c(0,1),addgrid.col = "grey",addCoef.col = "black")
    }else if(ncol(tru)>20){
        plt<-corrplot::corrplot(rel,method="color",
                                title="Original Network Edge Weight Reliabilities",
                                mar=c(2,2,2,2),tl.col="black",tl.cex=.75,
                                cl.lim = c(0,1),addgrid.col = "grey")}
    
    
    ############################################################################################
    #create boot reliability plot matrix
    if(bootmat==TRUE)
    {
    row.names(brel)<-colnames(brel)
    diag(brel)<-1
    upp<-matrix(0,nrow=nrow(brel),ncol=ncol(brel))
    for(i in 1:nrow(brel))
        for(j in 1:ncol(brel))
            if(brel[i,j]!=0&&boot[i,j]!=0)
            {upp[i,j]<-brel[i,j]}
    colnames(upp)<-colnames(brel)
    brel[upper.tri(brel)]<-upp[upper.tri(upp)]
    
    #reliablity on correlation plot
    x<-matrix(nrow=length(upp))
    y<-matrix(nrow=length(boot))
    diag(upp)<-0
    diag(boot)<-0
    wc<-0
    for(i in 1:nrow(boot))
        for(j in 1:ncol(boot))
            if((upp[i,j]!=0&&boot[i,j])!=0)
            {wc<-wc+1
            x[wc]<-upp[i,j]
            y[wc]<-boot[i,j]}
    xo<-na.omit(abs(x))
    yo<-na.omit(abs(y))
    
    dev.new()
    mar=c(2,2,2,2)
    cpo<-{plot(xo,yo,pch=16,ylab="Edge Weight Strength",xlab="Reliability",
               main="bootgen Network Edge Weight Strength on Reliability",xlim=c(0,1),ylim=range(yo))
        abline(lm(yo~xo))
        text(x=.05,y=max(yo-.05),labels = paste("r = ",round(cor(yo,xo),3)))}
    
    
    #plot reliability matrix
    dev.new()
    if(ncol(boot)<=20)
    {plt<-corrplot::corrplot(brel,method="color",
                             title="bootgen Network Edge Weight Reliabilities",
                             mar=c(2,2,2,2),tl.col="black",tl.cex=.75,
                             cl.lim = c(0,1),addgrid.col = "grey",addCoef.col = "black")
    }else if(ncol(boot)>20){
        plt<-corrplot::corrplot(brel,method="color",
                                title="bootgen Network Edge Weight Reliabilities",
                                mar=c(2,2,2,2),tl.col="black",tl.cex=.75,
                                cl.lim = c(0,1),addgrid.col = "grey")}
    
    dev.new()
    hist(na.omit(as.vector(ifelse(tru!=0,rel,NA))),
         breaks=20,xlab="Reliability",main="Reliability Frequency",xlim=c(0,1),
         col = rgb(0,0,1,0.5))
    hist(na.omit(as.vector(ifelse(boot!=0,brel,NA))),breaks=20,
         add=TRUE, col = rgb(1,0,0,0.5))
    legend("topleft",c("orignet","bootgen","overlap"),fill=c(rgb(0,0,1,0.5),rgb(1,0,0,0.5),"maroon3"))
    
    
    }
}
#----
#' Bootstrapped Communities Likelihood
#' @description Bootstraps the sample with replace to compute walktrap reliability
#' @param data A set of data
#' @param normal Should data be transformed to a normal distribution?
#' Defaults to FALSE.
#' Data is not transformed to be normal.
#' Set to TRUE if data should be transformed to be normal
#' @param n Number of people to use in the bootstrap.
#' Defaults to full sample size
#' @param iter Number of bootstrap iterations.
#' Defaults to 100 iterations
#' @param filter Set filter method.
#' Defaults to "TMFG"
#' @param method Defaults to "walktrap".
#' Set to "louvain" for the louvain community detection algorithm
#' @param na.data How should missing data be handled?
#' For "listwise" deletion the \emph{na.omit} function is applied.
#' Set to "fiml" for Full Information Maxmimum Likelihood (\emph{psych} package).
#' Full Information Maxmimum Likelihood is \strong{recommended} but time consuming
#' @param steps Number of steps to use in the walktrap algorithm.
#' Defaults to 4.
#' Use a larger number of steps for smaller networks
#' @param seeds Seeds to use for random number generation.
#' Defaults to NULL.
#' Input seeds from previous run (\strong{see examples})
#' @param ... Additional arguments for network filtering methods
#' @return Returns a list with the factors and their proportion found across bootstrapped samples
#' (proportion; i.e., their likelihood) and the seeds used in the random number generator (Seeds)
#' @examples
#' \dontrun{
#' commTMFG<-commboot(neoOpen)
#' 
#' #Replicates previous analysis
#' repcommTMFG<-commboot(neoOpen, seeds = commTMFG$Seeds)
#' 
#' commThreshold<-commboot(neoOpen,filter="threshold")
#' }
#' @references
#' Blondel, V. D., Guillaume, J. L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks.
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}(10), P10008.
#'
#' Csardi, G., & Nepusz, T. (2006).
#' The igraph software package for complex network research.
#' \emph{InterJournal, Complex Systems}, \emph{1695}(5), 1-9.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Bootstrapped Community Reliability----
commboot <- function (data, normal = FALSE, n = nrow(data), iter = 100,
                      filter = c("TMFG","threshold"),
                      method = c("louvain","walktrap"),
                      na.data = c("pairwise","listwise","fiml","none"), steps = 4,
                      seeds = NULL, ...)
{
    if(missing(filter))
    {filter<-"TMFG"
    }else{filter<-match.arg(filter)}
    
    if(missing(method))
    {method<-"walktrap"
    }else{method<-match.arg(method)}
    
    if(!is.null(seeds))
    {
        if(length(seeds)!=iter)
        {iter<-length(seeds)}
    }
    
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
        }else{data<-psych::corFiml(data)}
    }
    
    #corrlation matrix
    if(nrow(data)==ncol(data)){cormat<-data
    }else if(normal){cormat<-qgraph::cor_auto(data)
    }else{cormat<-cor(data)}
    
    n<-nrow(data)
    col<-ncol(data)
    
    mat<-matrix(0,nrow=n,ncol=col) #Initialize bootstrap matrix
    comm<-matrix(0,nrow=iter,ncol=1) #Initialize community matrix
    Seeds<-vector(mode="numeric",length=iter) #Initiliaze seeds vector
    
    pb <- txtProgressBar(max=iter, style = 3)
    for(i in 1:iter) #Generate array of bootstrapped samples
    {
        if(is.null(seeds))
        {f<-round(runif(i,min=1,max=1000000),0)
        Seeds[i]<-f[round(runif(i,min=1,max=length(f)),0)][i]
        }else{Seeds<-as.vector(seeds)}
        set.seed(Seeds[i])
        mat<-data[sample(1:n,replace=TRUE),]
        if(any(colSums(mat)<=1)){stop("Increase sample size: not enough observations")}
        
        #missing data handling
        if(any(is.na(mat)))
        {
            if(na.data=="pairwise")
            {
                if(normal)
                {cormat<-qgraph::cor_auto(mat,missing=na.data)
                }else{cormat<-psych::cor2(mat,use=na.data)}
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
        }
        
        if(nrow(mat)==ncol(mat)){cormat<-mat
        }else if(normal){cormat<-qgraph::cor_auto(mat)
        }else{cormat<-cor(mat)}
        
        cormat<-cor(mat)
        
        if(method=="walktrap")
        {if(filter=="TMFG")
        {comm[i,]<-max(igraph::walktrap.community(convert2igraph(TMFG(cormat,...)$A),steps=steps)$membership)
        }else if(filter=="threshold")
        {comm[i,]<-max(igraph::walktrap.community(convert2igraph(threshold(cormat,...)$A),steps=steps)$membership)
        }}else if(method=="louvain")
        {if(filter=="TMFG")
        {comm[i,]<-max(suppressWarnings(louvain(TMFG(cormat,...)$A)$community))
        }else if(filter=="threshold")
        {comm[i,]<-max(suppressWarnings(louvain(threshold(cormat,...)$A)$community))}
        }
        setTxtProgressBar(pb, i)
    }
    close(pb)
    
    count<-0
    prop<-matrix(0,nrow=length(seq(from=min(comm),to=max(comm))),ncol=1)
    for(i in min(comm):max(comm))
    {
        count<-count+1
        prop[count,]<-length(which(comm==i))
    }
    
    prop<-round(prop/iter,3)
    prop<-cbind(seq(from=min(comm),to=max(comm)),prop)
    colnames(prop)<-c("Factors","Likelihood")
    
    return(list(proportions=prop,Seeds=Seeds))
}
#----
#' Dependency Network Approach
#' @description Generates a dependency matrix of the data (index argument is still in testing phase)
#' @param data A set of data
#' @param normal Should data be transformed to a normal distribution?
#' Defaults to FALSE. Data is not transformed to be normal.
#' Set to TRUE if data should be transformed to be normal
#' (computes correlations using the \emph{cor_auto} function from the \emph{qgraph} package)
#' @param na.data How should missing data be handled?
#' For "listwise" deletion the \emph{na.omit} function is applied.
#' Set to "fiml" for Full Information Maxmimum Likelihood (\emph{psych} package).
#' Full Information Maxmimum Likelihood is \strong{recommended} but time consuming
#' @param index Should correlation with the latent variable
#' (i.e., weighted average of all variables) be removed?
#' Defaults to FALSE.
#' Set to TRUE to remove common latent factor
#' @param fisher Should Fisher's Z-test be used to keep significantly higher influences (index only)?
#' Defaults to FALSE.
#' Set to TRUE to remove non-significant influences
#' @param progBar Should progress bar be displayed?
#' Defaults to TRUE.
#' Set to FALSE for no progress bar
#' @return Returns an adjacency matrix of dependencies
#' @examples
#' D<-depend(neoOpen)
#' 
#' Dindex<-depend(neoOpen,index=TRUE)
#' @references
#' Kenett, D. Y., Tumminello, M., Madi, A., Gur-Gershgoren, G., Mantegna, R. N., & Ben-Jacob, E. (2010).
#' Dominating clasp of the financial sector revealed by partial correlation analysis of the stock market.
#' \emph{PloS one}, \emph{5}(12), e15032.
#' 
#' Kenett, D. Y., Huang, X., Vodenska, I., Havlin, S., & Stanley, H. E. (2015).
#' Partial correlation analysis: Applications for financial markets.
#' \emph{Quantitative Finance}, \emph{15}(4), 569-578.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
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
        }else{data<-psych::corFiml(data)}
    }
    
    #corrlation matrix
    if(nrow(data)==ncol(data)){cormat<-data
    }else if(normal){cormat<-qgraph::cor_auto(data)
    }else{cormat<-cor(data)}
    
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
            }else{cordat<-psych::cor2(dat,use=na.data)}
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
            }else{dat<-psych::corFiml(dat)}
        }
        
        #corrlation matrix
        if(nrow(dat)==ncol(dat)){cordat<-dat
        }else if(normal){cordat<-qgraph::cor_auto(dat)
        }else{cordat<-cor(dat)}
        
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
#' Split sample
#' @description Randomly splits sample into equally divded sub-samples
#' @param data Must be a dataset
#' @param samples Number of samples to produce (Defaults to 10)
#' @param splits Number to divide the sample by (Defaults to 2; i.e., split-half)
#' @return A list containing the training (trainSample) and testing (testSample) sizes
#' (training sample is made to have slightly larger sizes in the case of uneven splits),
#' and their respective sample sizes (trainSize and testSize)
#' @examples
#' splithalf<-splitsamp(neoOpen)
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Split samples----
splitsamp <- function (data, samples = 10, splits = 2)
{
    data<-as.matrix(data)
    n<-nrow(data)
    spli<-n/splits
    
    if(spli%%1!=0)
    {
        spl<-array()
        for(i in 1:splits)
        {spl[i]<-round(spli,0)}
        dif<-sum(spl)-n
        if(dif<0)
        {
            for(i in seq(from=2,to=splits,by=2))
            {
                spl[i]<-spl[i]+1
                if(sum(spl)==n){break}
            }
        }else if(dif>0)
        {
            for(i in seq(from=2,to=splits,by=2))
            {
                spl[i]<-spl[i]-1
                if(sum(spl)==n){break}
            }
        }
    }else{spl<-rep(spli,splits)}

    samps<-list()
    for(i in 1:samples)
    {
        samps[[i]]<-as.data.frame(cbind(sample(1:n),data))
        samps[[i]]<-samps[[i]][order(samps[[i]]$V1),]
    }
    
    spl<-sort(spl,decreasing = TRUE)
    s<-splits/2
    
    tr<-seq(from=1,to=s,by=1)
    m<-c(1,cumsum(spl[tr]))
    
    if(s%%1!=0)
    {m<-c(m,ceiling(n/2))}

    count<-0

    train<-list()
    for(i in 1:samples)
        for(k in 1:(length(m)-1))
        {
            count<-count+1
            if(k==1)
            {train[[count]]<-samps[[i]][m[k]:m[k+1],-1]
            }else if(k>1){
                train[[count]]<-samps[[i]][(m[k]+1):m[k+1],-1]
            }
        }
    
    d<-diff(m)
    d[1]<-d[1]+1
    
    trainSize<-rep(d,samples)
    
    te<-seq(from=max(tr)+1,to=splits,by=1)
    if(s%%1!=0)
    {
        if(splits>3)
        {te<-seq(from=max(tr)+2,to=splits,by=1)
        }else{te<-3}
    }else if(s%%1==0){te<-seq(from=max(tr)+1,to=splits,by=1)}
    
    m<-m[length(m)]
    m[1]<-m[1]+1
    m<-c(m,cumsum(spl[te])+(m-1))
    
    if(s%%1!=0)
    {m<-c(m,n)}
    
    count<-0
    
    test<-list()
    for(i in 1:samples)
        for(k in 1:(length(m)-1))
        {
            count<-count+1
            if(k==1)
            {test[[count]]<-samps[[i]][m[k]:m[k+1],-1]
            }else if(k>1){
                test[[count]]<-samps[[i]][(m[k]+1):m[k+1],-1]
            }
        }
    
    d<-diff(m)
    d[1]<-d[1]+1
    
    testSize<-rep(d,samples)
    
    return(list(trainSamples=train,trainSize=trainSize,testSamples=test,testSize=testSize))
}
#----
#' Network Construction for splitsamp
#' @description Applies a network construction method to samples generated from the splitsamp function
#' @param object Object output from \emph{splitsamp} function
#' @param method A network filtering method.
#' Defaults to "TMFG"
#' @param bootmat Should \emph{bootgen} function be used?
#' Defaults to FALSE.
#' Set to TRUE to constuction bootgen networks
#' @param seeds Seeds to use for bootmat.
#' Defaults to NULL.
#' Input seeds from previous bootgen run
#' @param pB Should progress bar be displayed?
#' Defaults to TRUE.
#' Set to FALSE for no progress bar
#' @param ... Additional arguments to be passed to the network construction methods
#' @return Returns a list of networks for training (train) and testing (test) samples
#' @examples
#' samples <- splitsamp(neoOpen)
#' 
#' nets <- splitsampNet(samples, method="TMFG")
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Network construction for split sample----
splitsampNet <- function (object, method = c("MaST", "PMFG", "TMFG", "LoGo", "threshold"),
                          bootmat = FALSE, seeds = NULL, pB = TRUE, ...)
{
    if(missing(method))
    {method<-"TMFG"
    }else{method<-match.arg(method)}
    
    samps<-length(object$trainSamples)
    
    train<-list()
    test<-list()
    
    if(pB)
    {pb <- txtProgressBar(max=samps, style = 3)}
    
    for(i in 1:samps)
    {
        if(!bootmat)
        {
            if(method=="MaST")
            {
                train[[i]]<-MaST(object$trainSamples[[i]],...)
                test[[i]]<-MaST(object$testSamples[[i]],...)
            }else if(method=="PMFG")
            {
                train[[i]]<-PMFG(object$trainSamples[[i]],progBar = FALSE,...)$pmfg
                test[[i]]<-PMFG(object$testSamples[[i]],progBar = FALSE,...)$pmfg
            }else if(method=="TMFG")
            {
                train[[i]]<-TMFG(object$trainSamples[[i]],...)$A
                test[[i]]<-TMFG(object$testSamples[[i]],...)$A
            }else if(method=="LoGo")
            {
                train[[i]]<-LoGo(object$trainSamples[[i]],partial=TRUE,...)$logo
                test[[i]]<-LoGo(object$testSamples[[i]],partial=TRUE,...)$logo
            }else if(method=="threshold")
            {
                train[[i]]<-threshold(object$trainSamples[[i]],...)$A
                test[[i]]<-threshold(object$testSamples[[i]],...)$A
            }else{stop("Method not available")}
            
        }else{
            train[[i]]<-bootgen(object$trainSamples[[i]],method=method,
                                seeds=seeds,progBar = FALSE)$bootmat
            test[[i]]<-bootgen(object$testSamples[[i]],method=method,
                               seeds=seeds,progBar = FALSE)$bootmat
        }  
        if(pB)
        {setTxtProgressBar(pb, i)}
    }
    if(pB)
    {close(pb)}
    
    return(list(train=train,test=test))
}
#----
#' Statistics for splitsamp Networks
#' @description Applies a network statistics to the networks generated from the splitsampNet function
#' @param object Object output from \emph{splitsampNet} function
#' @param stats Which stats should be computed?
#' "edges" will compute statistics for the \emph{edgerep} function.
#' "centralities" will compute betweenness, closeness, and strength.
#' Both options are allowed
#' @param edges If stats = "edges, should the lower or upper replication percentage be used?
#' Defaults to "lower".
#' @return Returns a list of statistics for edges (Edges), centralities (Centralities), or both
#' @examples
#' samples <- splitsamp(neoOpen)
#' 
#' nets <- splitsampNet(samples, method="TMFG")
#' 
#' stats <- splitsampStats(nets, stats = "edges", edges = "lower")
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Stats for split sample----
splitsampStats <- function (object, stats = c("edges","centralities"),
                            edges = c("lower","upper"))
{
    if(missing(stats))
    {stats<-c("edges","centralities")
    }else{stats<-match.arg(stats,several.ok = TRUE)}
    
    if(missing(edges))
    {edges<-"lower"
    }else{edges<-match.arg(edges)}
    
    samps<-length(object$train)
    
    if(any(stats=="edges"))
    {
        edge <- list()
        repli <- vector(mode="numeric",length=samps)
        percent <- vector(mode="numeric",length=samps)
        numA <- vector(mode="numeric",length=samps)
        numB <- vector(mode="numeric",length=samps)
    }
    
    if(any(stats=="centralities"))
    {
        cent <- list()
        betw <- vector(mode="numeric",length=samps)
        clos <- vector(mode="numeric",length=samps)
        stre <- vector(mode="numeric",length=samps)
    }
    
    for(i in 1:samps)
    {
        if(any(stats=="edges"))
        {
            edge[[i]] <- edgerep(object$train[[i]],object$test[[i]])
            repli[i] <- edge[[i]]$replicated
            
            if(edges=="lower")
            {percent[i] <- ifelse(edge[[i]]$percentageA >= edge[[i]]$percentageB,edge[[i]]$percentageB,edge[[i]]$percentageA)
            }else if(edges=="upper")
            {percent[i] <- ifelse(edge[[i]]$percentageA <= edge[[i]]$percentageB,edge[[i]]$percentageB,edge[[i]]$percentageA)}
            
            numA[i] <- edge[[i]]$totalEdgesA
            numB[i] <- edge[[i]]$totalEdgesB
        }
        
        if(any(stats=="centralities"))
        {
            betw[i] <- cor(betweenness(object$train[[i]]),betweenness(object$test[[i]]))
            clos[i] <- cor(closeness(object$train[[i]]),closeness(object$test[[i]]))
            stre[i] <- cor(strength(object$train[[i]]),strength(object$test[[i]]))
        }
    }
    
    
    
    if(all(stats==c("edges","centralities")))
    {
        Edges<-list(EdgeRep=edge,Replicated=repli,Percent=percent,nEdgesTrain=numA,nEdgesTest=numB)
        Centralities<-list(BC=betw,LC=clos,Str=stre)
        
        return(list(Edges=Edges,Centralities=Centralities))
        
    }else if(any(stats=="edges"))
    {
        Edges<-list(EdgeRep=edge,Replicated=repli,Percent=percent,nEdgesTrain=numA,nEdgesTest=numB)
        
        return(Edges)
        
    }else if(any(stats=="centralities"))
    {
        Centralities<-list(BC=betw,LC=clos,Str=stre)
        
        return(Centralities)
    }
}
#----
#' Generates a Random Network
#' @description Generates a random network
#' @param nodes Number of nodes in random network
#' @param edges Number of edges in random network
#' @return Returns an adjacency matrix of a random network
#' @examples
#' rand <- randnet(10,27)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Random Network----
randnet <- function (nodes, edges)
{
    mat<-matrix(1,nrow=nodes,ncol=nodes)
    diag(mat)<-0
    ind<-ifelse(upper.tri(mat)==TRUE,1,0)
    i<-which(ind==1)
    rp<-sample(length(i))
    irp<-i[rp]
    
    rand<-matrix(0,nrow=nodes,ncol=nodes)
    rand[irp[1:edges]]<-1
    rand<-rand+t(rand)
    
    return(rand)
}
#----
#' Generates a Lattice Network
#' @description Generates a lattice network
#' @param nodes Number of nodes in lattice network
#' @param edges Number of edges in lattice network
#' @return Returns an adjacency matrix of a lattice network
#' @examples
#' latt <- lattnet(10,27)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Lattice Network----
lattnet <- function (nodes, edges)
{
    dlat<-matrix(0,nrow=nodes,ncol=nodes)
    lat<-matrix(0,nrow=nodes,ncol=nodes)
    
    for(i in 1:nodes)
    {
        if(i!=nodes)
        {dlat[i,i+1]<-1}
        if(i<nodes-1)
        {dlat[i,i+2]<-1}
    }
    lat<-dlat+t(dlat)
    
    over<-sum(lat)-edges
    if(over!=0)
    {rp<-sample(which(dlat==1))
    for(i in 1:over)
    {lat[rp[i]]<-0}}
    
    return(lat)   
}
#----
#' Binarize Network
#' @description Converts weighted adjacency matrix to a binarized adjacency matrix
#' @param A An adjacency matrix of network data (or an array of matrices)
#' @return Returns an adjancency matrix of 1's and 0's
#' @examples
#' net<-TMFG(neoOpen)$A
#' 
#' psyb<-binarize(net)
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Binarize function----
binarize <- function (A)
{
    bin<-ifelse(A!=0,1,0)
    
    return(bin)
}
#----
#' Regression Matrix
#' @description Computes regression such that one variable is regressed over all other variables
#' @param data A dataset
#' @param family Error distribution to be used in the regression model.
#' Defaults to "logistic".
#' Set to any family used in function \emph{glm} (see \code{\link{family}})
#' @param symmetric Should matrix be symmetric?
#' Defaults to TRUE, taking the mean of the two edge weights
#' (i.e., [\emph{i},\emph{j}] and [\emph{j},\emph{i}])
#' Set to FALSE for asymmwtric weights
#' (i.e., [\emph{i},\emph{j}] does not equal [\emph{j},\emph{i}])
#' @return A matrix of fully regressed coefficients where one variable is regressed over all others
#' @examples
#' #binarize responses
#' psyb <- ifelse(neoOpen>=4, 1, 0)
#' 
#' #perform logistic regression
#' mat <- reg(psyb)
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @importFrom stats glm
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
#' Convert Network(s) to igraph's Format
#' @description Converts single or multiple networks into \emph{igraph}'s format for network analysis
#' @param A Adjacency matrix (network matrix) or brain connectivity array (from \emph{convertConnBrainMat})
#' @param neural Is input a brain connectivity array (i.e., m x m x n)?
#' Defaults to FALSE.
#' Set to TRUE to convert each brain connectivity matrix
#' @return Returns a network matrix in \emph{igraph}'s format or
#' returns a list of brain connectivity matrices each of which have been convert to \emph{igraph}'s format
#' @examples
#' A <- TMFG(neoOpen)$A
#' 
#' igraphNetwork <- convert2igraph(A)
#' 
#' \dontrun{ 
#' neuralarray<-convertConnBrainMat()
#' 
#' igraphNeuralList <- convert2igraph(neuralarray, neural = TRUE)
#' }
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Convert matrices to igraph format----
convert2igraph <- function (A, neural = FALSE)
{
    if(!neural)
    {
        net <- igraph::as.igraph(qgraph::qgraph(A,DoNotPlot=TRUE))
        return(net)
    }else if(neural)
    {
        netlist <- list()
        
        n<-length(A)/nrow(A)/ncol(A)
        
        pb <- txtProgressBar(max=n, style = 3)
        
        for(i in 1:n)
        {
            netlist[[i]] <- igraph::as.igraph(qgraph::qgraph(A[,,i],DoNotPlot=TRUE))
            setTxtProgressBar(pb, i)
        }
        close(pb)
        
        return(netlist)
    }
}
#----
#' Import CONN Toolbox Brain Matrices to R format
#' @description Converts a Matlab brain z-score connectivity array (n x n x m) where \strong{n} is the n x n connectivity matrices and \strong{m} is the participant.
#' If you would like to simply import a connectivity array from Matlab, then see the examples
#' @param MatlabData Input for Matlab data file.
#' Defaults to interactive file choice
#' @param progBar Should progress bar be displayed?
#' Defaults to TRUE.
#' Set FALSE for no progress bar
#' @return Returns a list of correlation (rmat) and z-score (zmat) arrays of neural connectivity matrices (n x n x m)
#' @examples
#' \dontrun{
#' neuralarray<-convertConnBrainMat()
#' 
#' #Import correlation connectivity array from Matlab
#' library(R.matlab)
#' neuralarray<-readMat(file.choose())
#' }
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Convert CONN Toolbox Brain Matrices----
convertConnBrainMat <- function (MatlabData, progBar = TRUE)
{
    if(missing(MatlabData))
    {mat<-R.matlab::readMat(file.choose())
    }else{mat<-R.matlab::readMat(MatlabData)}
    
    if(!is.list(mat))
    {return(mat)
    }else
    
    #read in matlab data
    n1<-nrow(mat$Z) #determine number of rows
    n2<-ncol(mat$Z) #determine number of columns
    if(nrow(mat$Z)!=ncol(mat$Z))
    {warning("Row length does not match column length")}
    m<-length(mat$Z)/n1/n2 #determine number of participants
    
    #change row and column names
    coln1<-matrix(0,nrow=n1) #get row names
    for(i in 1:n1)
    {coln1[i,]<-mat$names[[i]][[1]][1,1]}

    coln2<-matrix(0,nrow=n2) #get column names
    for(i in 1:n2)
    {coln2[i,]<-mat$names2[[i]][[1]][1,1]}
    
    dat<-mat$Z
    if(progBar)
    {pb <- txtProgressBar(max=m, style = 3)}
    
    for(i in 1:m) #populate array
    {
        dat[,,i]<-psych::fisherz2r(mat$Z[,,i])
        for(j in 1:n1)
            for(k in 1:n2)
                if(is.na(dat[j,k,i]))
                {dat[j,k,i]<-0}
        if(progBar){setTxtProgressBar(pb, i)}
    }
    if(progBar){close(pb)}
    
    colnames(dat)<-coln2
    row.names(dat)<-coln1
    
    return(list(rmat=dat,zmat=mat$Z))
}
#----
#' Dependency Neural Networks
#' @description Applies the dependency network approach to neural network array
#' @param neuralarray Array from \emph{convertConnBrainMat} function
#' @param pB Should progress bar be displayed?
#' Defaults to TRUE.
#' Set FALSE for no progress bar
#' @param ... Additional arguments from \emph{depend} function
#' @return Returns an array of n x n x m dependency matrices
#' @examples
#' \dontrun{
#' neuralarray <- convertConnBrainMat()
#' 
#' dependencyneuralarray <- depna(neuralarray)
#' }
#' @references
#' Jacob, Y., Winetraub, Y., Raz, G., Ben-Simon, E., Okon-Singer, H., Rosenberg-Katz, K., ... & Ben-Jacob, E. (2016).
#' Dependency Network Analysis (DEPNA) Reveals Context Related Influence of Brain Network Nodes.
#' \emph{Scientific Reports}, \emph{6}, 27444.
#' 
#' Kenett, D. Y., Tumminello, M., Madi, A., Gur-Gershgoren, G., Mantegna, R. N., & Ben-Jacob, E. (2010).
#' Dominating clasp of the financial sector revealed by partial correlation analysis of the stock market.
#' \emph{PloS one}, \emph{5}(12), e15032.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Dependency Network Analysis----
depna <- function (neuralarray, pB = TRUE, ...)
{
    n<-length(neuralarray)/nrow(neuralarray)/ncol(neuralarray)
    
    for(i in 1:n)    
        if(nrow(neuralarray)!=ncol(neuralarray))
        {stop(paste("Participant ",i,"'s matrix is not symmetric",sep=""))}
    
    deparray<-neuralarray
    
    if(pB)
    {pb <- txtProgressBar(max=n, style = 3)}
    
    for(i in 1:n)
    {deparray[,,i]<-depend(neuralarray[,,i],progBar=FALSE,...)
    if(pB){setTxtProgressBar(pb, i)}}
    
    if(pB){close(pb)}
    
    return(deparray)
}
#----
#' Neural Network Filter
#' @description Applies a network filtering methodology to neural network array.
#' Removes edges from the neural network output from \emph{convertConnBrainMat} using a network filtering approach
#' @param neuralarray Array from \emph{convertConnBrainMat} function
#' @param method Filtering method to be applied
#' @param progBar Should progress bar be displayed?
#' Defaults to TRUE.
#' Set FALSE for no progress bar
#' @param ... Additional arguments from filtering methods
#' @return Returns an array of n x n x m filtered matrices
#' @examples
#' \dontrun{neuralarray <- convertConnBrainMat()
#' 
#' filteredneuralarray <- neuralnetfilter(neuralarray, method = "threshold", thresh = .50)
#' 
#' dependencyarray <- depna(neuralarray)
#' 
#' filtereddependencyarray <- neuralnetfilter(dependencyarray, method = "TMFG", depend = TRUE)
#' }
#' @references
#' Fallani, F. D. V., Latora, V., & Chavez, M. (2017).
#' A topological criterion for filtering information in complex brain networks.
#' \emph{PLoS Computational Biology}, \emph{13}(1), e1005305.
#' 
#' Massara, G. P., Di Matteo, T., & Aste, T. (2016).
#' Network filtering for big data: Triangulated maximally filtered graph.
#' \emph{Journal of Complex Networks}, \emph{5}(2), 161-178. 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Neural Network Filter----
neuralnetfilter <- function (neuralarray, method = c("TMFG","MaST","ECOplusMaST","ECO","threshold"),progBar = TRUE, ...)
{
    if(missing(method))
    {method<-"TMFG"
    }else{method<-match.arg(method)}
    
    n<-length(neuralarray)/nrow(neuralarray)/ncol(neuralarray)  
    
    for(i in 1:n)    
        if(nrow(neuralarray)!=ncol(neuralarray))
        {stop(paste("Participant ",i,"'s matrix is not symmetric",sep=""))}
    
    filarray<-neuralarray
    
    if(progBar)
    {pb <- txtProgressBar(max=n, style = 3)}
    
    for(i in 1:n)
    {
        if(method=="TMFG")
        {filarray[,,i]<-TMFG(neuralarray[,,i],...)$A
        }else if(method=="MaST")
        {filarray[,,i]<-MaST(neuralarray[,,i],...)
        }else if(method=="ECO")
        {filarray[,,i]<-ECO(neuralarray[,,i],...)
        }else if(method=="ECOplusMaST")
        {filarray[,,i]<-ECOplusMaST(neuralarray[,,i],...)
        }else if(method=="threshold")
        {filarray[,,i]<-threshold(neuralarray[,,i],...)$A
        }else{stop("Method not available")}
        if(progBar){setTxtProgressBar(pb, i)}
    }
    if(progBar){close(pb)}
    
    return(filarray)
}
#----
#' Local and Global Neural Network Characteristics 
#' @description Obtains a global or local network characteristic from neural network data
#' @param filarray Filtered array from \emph{neuralnetfilter} function
#' @param statistic A statistic to compute
#' @param progBar Should progress bar be displayed? Defaults to TRUE. Set FALSE for no progress bar
#' @param ... Additional arguments for statistics functions
#' @return Returns vector of global characteristics (rows = participants, columns = statistic) or a matrix of local characteristics (rows = ROIs, columns = participants)
#' @examples
#' \dontrun{
#' neuralarray <- convertConnBrainMat()
#' 
#' filteredneuralarray <- neuralnetfilter(neuralarray, method = "threshold", thres = .50)
#' 
#' ClusteringCoefficient <- neuralstat(filteredneuralarray, statistic = "CC")
#' 
#' AverageShortestPathLength <- neuralstat(filteredneuralarray, statistic = "ASPL")
#' 
#' Modularity <- neuralstat(filteredneuralarray, statistic = "Q")
#' 
#' Smallworldness <- neuralstat(filteredneuralarray, statistic = "S")
#' 
#' Trasitivity <- neuralstat(filteredneuralarray, statistic = "transitivity")
#' 
#' Connectivity <- neuralstat(filteredneuralarray, statistic = "conn")
#' 
#' BetweennessCentrality <- neuralstat(filteredneuralarray, statistic = "BC")
#' 
#' ClosenessCentrality <- neuralstat(filteredneuralarray, statistic = "LC")
#' 
#' Degree <- neuralstat(filteredneuralarray, statistic = "deg")
#' 
#' NodeStrength <- neuralstat(filteredneuralarray, statistic = "str")
#' 
#' Communities <- neuralstat(filteredneuralarray, statistic = "comm")
#' 
#' EigenvectorCentrality <- neuralstat(filteredneuralarray, statistic = "EC")
#' 
#' LeverageCentrality <- neuralstat(filteredneuralarray, statistic = "lev")
#' 
#' RandomShortestPathBC <- neuralstat(filteredneuralarray, statistic = "rspbc")
#' 
#' HybridCentrality <- neuralstat(filteredneuralarray, statistic = "hybrid")
#' 
#' NodeImpact <- neuralstat(filteredneuralarray, statistic = "impact")
#' }
#' @references
#' Blondel, V. D., Guillaume, J. L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks. 
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}(10), P10008.
#' 
#' Joyce, K. E., Laurienti, P. J., Burdette, J. H., & Hayasaka, S. (2010).
#' A new measure of centrality for brain networks. 
#' \emph{PLoS One}, \emph{5}(8), e12200.
#' 
#' Kenett, Y. N., Kenett, D. Y., Ben-Jacob, E., & Faust, M. (2011).
#' Global and local features of semantic networks: Evidence from the Hebrew mental lexicon.
#' \emph{PloS one}, \emph{6}(8), e23912.
#' 
#' Kivimaki, I., Lebichot, B., Saramaki, J., & Saerens, M. (2016).
#' Two betweenness centrality measures based on Randomized Shortest Paths.
#' \emph{Scientific Reports}, \emph{6}(19668), 1-15.
#' 
#' Pozzi, F., Di Matteo, T., & Aste, T. (2013).
#' Spread of risk across financial markets: Better to invest in the peripheries. 
#' \emph{Scientific Reports}, \emph{3}(1655), 1-7.
#' 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}(3), 1059-1069. 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @importFrom stats t.test
#' @export
#Neural Statistics----
neuralstat <- function (filarray, statistic = c("CC","ASPL","Q","S","transitivity","conn","BC","LC",
                                                "deg","inDeg","outDeg","degRI","str","inStr","outStr",
                                                "strRI","comm","EC","lev","rspbc","hybrid","impact"), progBar = TRUE, ...)
{
    statistic<-match.arg(statistic)
    
    n<-length(filarray)/ncol(filarray)/nrow(filarray)
    
    gstat<-matrix(0,nrow=n,ncol=1)
    lstat<-matrix(0,nrow=ncol(filarray),ncol=n)
    
    if(progBar)
    {pb <- txtProgressBar(max=n, style = 3)}
    
    for(i in 1:n)
    {
        if(statistic=="CC")
        {gstat[i,]<-clustcoeff(filarray[,,i],...)$CC
        gstat<-as.data.frame(gstat)
        colnames(gstat)<-"CC"
        row.names(gstat)[i]<-paste("sub",i,sep="")
        gstat<-as.matrix(gstat)
        }else if(statistic=="ASPL")
        {gstat[i,]<-pathlengths(filarray[,,i],...)$ASPL
        gstat<-as.data.frame(gstat)
        colnames(gstat)<-"ASPL"
        row.names(gstat)[i]<-paste("sub",i,sep="")
        gstat<-as.matrix(gstat)
        }else if(statistic=="Q")
        {gstat[i,]<-louvain(filarray[,,i],...)$Q
        gstat<-as.data.frame(gstat)
        colnames(gstat)<-"Q"
        row.names(gstat)[i]<-paste("sub",i,sep="")
        gstat<-as.matrix(gstat)
        }else if(statistic=="S")
        {gstat[i,]<-smallworldness(filarray[,,i],iter=10,method="HG",...)$swm
        gstat<-as.data.frame(gstat)
        colnames(gstat)<-"S"
        row.names(gstat)[i]<-paste("sub",i,sep="")
        gstat<-as.matrix(gstat)
        }else if(statistic=="transitivity")
        {gstat[i,]<-transitivity(filarray[,,i],...)
        gstat<-as.data.frame(gstat)
        colnames(gstat)<-"transitivity"
        row.names(gstat)[i]<-paste("sub",i,sep="")
        gstat<-as.matrix(gstat)
        }else if(statistic=="conn")
        {gstat[i,]<-conn(filarray[,,i],...)$total
        gstat<-as.data.frame(gstat)
        colnames(gstat)<-"connectivity"
        row.names(gstat)[i]<-paste("sub",i,sep="")
        gstat<-as.matrix(gstat)
        }else if(statistic=="BC")
        {lstat[,i]<-betweenness(filarray[,,i],...)[,1]
        lstat<-as.data.frame(lstat)
        colnames(lstat)[i]<-paste("sub",i,sep="")
        row.names(lstat)<-colnames(filarray)
        lstat<-as.matrix(lstat)
        }else if(statistic=="LC")
        {lstat[,i]<-closeness(filarray[,,i],...)[,1]
        lstat<-as.data.frame(lstat)
        colnames(lstat)[i]<-paste("sub",i,sep="")
        row.names(lstat)<-colnames(filarray)
        lstat<-as.matrix(lstat)
        }else if(statistic=="deg")
        {lstat[,i]<-degree(filarray[,,i])
        lstat<-as.data.frame(lstat)
        colnames(lstat)[i]<-paste("sub",i,sep="")
        row.names(lstat)<-colnames(filarray)
        lstat<-as.matrix(lstat)
        }else if(statistic=="indeg")
        {lstat[,i]<-degree(filarray[,,i])$inDegree
        lstat<-as.data.frame(lstat)
        colnames(lstat)[i]<-paste("sub",i,sep="")
        row.names(lstat)<-colnames(filarray)
        lstat<-as.matrix(lstat)
        }else if(statistic=="outdeg")
        {lstat[,i]<-degree(filarray[,,i])$outDegree
        lstat<-as.data.frame(lstat)
        colnames(lstat)[i]<-paste("sub",i,sep="")
        row.names(lstat)<-colnames(filarray)
        lstat<-as.matrix(lstat)
        }else if(statistic=="degRI")
        {lstat[,i]<-degree(filarray[,,i])$relInf
        lstat<-as.data.frame(lstat)
        colnames(lstat)[i]<-paste("sub",i,sep="")
        row.names(lstat)<-colnames(filarray)
        lstat<-as.matrix(lstat)
        }else if(statistic=="str")
        {lstat[,i]<-strength(filarray[,,i])
        lstat<-as.data.frame(lstat)
        colnames(lstat)[i]<-paste("sub",i,sep="")
        row.names(lstat)<-colnames(filarray)
        lstat<-as.matrix(lstat)
        }else if(statistic=="instr")
        {lstat[,i]<-strength(filarray[,,i])$inStrength
        lstat<-as.data.frame(lstat)
        colnames(lstat)[i]<-paste("sub",i,sep="")
        row.names(lstat)<-colnames(filarray)
        lstat<-as.matrix(lstat)
        }else if(statistic=="outstr")
        {lstat[,i]<-strength(filarray[,,i])$outStrength
        lstat<-as.data.frame(lstat)
        colnames(lstat)[i]<-paste("sub",i,sep="")
        row.names(lstat)<-colnames(filarray)
        lstat<-as.matrix(lstat)
        }else if(statistic=="strRI")
        {lstat[,i]<-strength(filarray[,,i])$relInf
        lstat<-as.data.frame(lstat)
        colnames(lstat)[i]<-paste("sub",i,sep="")
        row.names(lstat)<-colnames(filarray)
        lstat<-as.matrix(lstat)
        }else if(statistic=="comm")
        {lstat[,i]<-suppressWarnings(louvain(filarray[,,i],...)$community)
        lstat<-as.data.frame(lstat)
        colnames(lstat)[i]<-paste("sub",i,sep="")
        row.names(lstat)<-colnames(filarray)
        lstat<-as.matrix(lstat)
        }else if(statistic=="EC")
        {lstat[,i]<-eigenvector(filarray[,,i],...)[,1]
        lstat<-as.data.frame(lstat)
        colnames(lstat)[i]<-paste("sub",i,sep="")
        row.names(lstat)<-colnames(filarray)
        lstat<-as.matrix(lstat)
        }else if(statistic=="lev")
        {lstat[,i]<-leverage(filarray[,,i],...)[,1]
        lstat<-as.data.frame(lstat)
        colnames(lstat)[i]<-paste("sub",i,sep="")
        row.names(lstat)<-colnames(filarray)
        lstat<-as.matrix(lstat)
        }else if(statistic=="rspbc")
        {lstat[,i]<-rspbc(filarray[,,i],...)[,1]
        lstat<-as.data.frame(lstat)
        colnames(lstat)[i]<-paste("sub",i,sep="")
        row.names(lstat)<-colnames(filarray)
        lstat<-as.matrix(lstat)
        }else if(statistic=="hybrid")
        {lstat[,i]<-hybrid(filarray[,,i],...)[,1]
        lstat<-as.data.frame(lstat)
        colnames(lstat)[i]<-paste("sub",i,sep="")
        row.names(lstat)<-colnames(filarray)
        lstat<-as.matrix(lstat)
        }else if(statistic=="impact")
        {lstat[,i]<-impact(filarray[,,i],...)[,1]
        lstat<-as.data.frame(lstat)
        colnames(lstat)[i]<-paste("sub",i,sep="")
        row.names(lstat)<-colnames(filarray)
        lstat<-as.matrix(lstat)
        }
        
        if(progBar){setTxtProgressBar(pb, i)}
        
    }
    
    if(progBar){close(pb)}
    
    if(sum(gstat)!=0)
    {return(gstat)}
    
    if(sum(lstat)!=0)
    {return(lstat)}
}
#----
#' Neural Network Group Statistics Tests
#' @description Statistical test for group differences for global or local network characteristics of neural network data
#' (still in testing phase)
#' @param groups Participant vector of desired groups (\strong{see examples})
#' @param nstat A statistic vector (whole-network) or matrix (ROI) from the \emph{neuralstat} function
#' @param correction Multiple comparisons correction for ROI testing.
#' Defaults to local false discovery rate (i.e., "FDR")
#' @return Returns test statistics for given groups and statistics
#' @examples
#' \dontrun{
#' neuralarray <- convertConnBrainMat()
#' 
#' filteredneuralarray <- neuralnetfilter(neuralarray, method = "threshold", thres = .50)
#' 
#' AverageShortestPathLength <- neuralstat(filteredneuralarray, statistic = "ASPL")
#' 
#' Degree <- neuralstat(filteredneuralarray, statistic = "deg")
#' 
#' groups <- c(rep(1,30),rep(2,30))
#' 
#' WholeNetwork_t-test <- neuralstattest(groups, AverageShortestPathLength)
#' 
#' ROI_t-test <- neuralstattest(groups, Degree)
#' }
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @importFrom stats coef
#' @export
#Neural Groups----
neuralgrouptest <- function (groups, nstat, correction = c("bonferroni","FDR"))
{
    if(missing(correction))
    {correction<-"FDR"
    }else{correction<-match.arg(correction)}
    
    if(max(groups)==2)
    {
        if(ncol(nstat)==1){
            compare<-cbind(groups,nstat)
            group1<-compare[which(compare[,1]==1),2]
            group2<-compare[which(compare[,1]==2),2]
            return(t.test(group1,group2,var.equal=TRUE))
        }else if(ncol(nstat)>1)
        {
            compare<-cbind(groups,t(nstat))
            group1<-compare[which(compare[,1]==1),2:ncol(compare)]
            group2<-compare[which(compare[,1]==2),2:ncol(compare)]
            count1<-0
            count2<-0
            ROI<-array()
            unpval<-array()
            cpval<-array()
            tval<-array()
            oval<-array()
            g1<-array()
            g2<-array()
            for(i in 1:ncol(group1))
            {count1<-count1+1
            oval[count1]<-t.test(group1[,i],group2[,i],var.equal = TRUE)$p.value}
            if(t.test(group1[,i],group2[,i],var.equal = TRUE)$p.value<=.05)
            {count2<-count2+1
            ROI[count2]<-colnames(group1)[count1]
            tval[count2]<-t.test(group1[,i],group2[,i],var.equal = TRUE)$statistic
            unpval[count2]<-t.test(group1[,i],group2[,i],var.equal = TRUE)$p.value
            if(correction=="bonferroni")
            {cpval[count2]<-t.test(group1[,i],group2[,i],var.equal = TRUE)$p.value<=(.05/ncol(group1))}
            if(correction=="FDR")
            {cpval<-fdrtool::fdrtool(oval,plot=FALSE,verbose=FALSE,statistic = "pvalue")$qval
            cpval<-ifelse(cpval<=.10,cpval,0)}
            g1[count2]<-t.test(group1[,i],group2[,i],var.equal = TRUE)$estimate[1]
            g2[count2]<-t.test(group1[,i],group2[,i],var.equal = TRUE)$estimate[2]}
        }
        sig<-cbind(ROI,round(tval,3),round(unpval,5),round(cpval,5),round(g1,3),round(g2,3))
        sig<-as.data.frame(sig)
        colnames(sig)<-c("ROI","t-value","uncorrected p","corrected p","group 1 mean", "group 2 mean")
        
        if(all(is.na(sig)))
        {print("No significant differences")
        }else if(any(sig$`corrected p`==0))
        {sig$`corrected p`[which(sig$`corrected p`==0)]<-"ns"
        return(sig)}else{sig$`corrected p`[which(sig$`corrected p`==1)]<-"sig"
        return(sig)}
        }
    
    
    #else if(max(groups)>2)
    #{
    #    if(ncol(nstat)>1){
    #        compare<-cbind(groups,t(nstat))
    #        compare<-as.data.frame(compare)
    #        compare$groups<-as.factor(compare$groups)
    #        
    #        count1<-0
    #        count2<-0
    #        ROI<-array()
    #        unpval<-array()
    #        tpval<-array()
    #        df<-array()
    #        res<-array()
    #        f<-array()
    #        tgroup<-array()
    #        
    #        for(i in 2:ncol(compare))
    #        {count1<-count1+1
    #        if(anova(aov(compare[,i]~groups,compare))[1,5]<=.05)
    #        {count2<-count2+1
    #        ROI[count2]<-colnames(compare)[count1]
    #        test<-anova(aov(compare[,i]~groups,compare))
    #        df[count2]<-test[1,1]
    #        res[count2]<-test[2,1]
    #        f[count2]<-test[1,4]
    #        unpval[count2]<-test[1,5]}
    #        if(any(TukeyHSD(aov(compare[,i]~groups,compare))$groups[,4]<.05))
    #        {
    #        for(j in 1:length(which(TukeyHSD(aov(compare[,i]~groups,compare))$groups[,4]<.05)))
    #        tpval[count2]<-which(TukeyHSD(aov(compare[,i]~groups,compare))$groups[,4]<.05)
    #        tgroup[count2]<-row.names(TukeyHSD(aov(compare[,4]~groups,compare))$groups)[which(TukeyHSD(aov(compare[,4]~groups,compare))$groups[j,4]<.05)]
    #        }
    #        }
    #        
    #        
    #        
    #        return(t.test(group1,group2,var.equal=TRUE))
    #    }else if(ncol(nstat)!=1)
}
#----
#' Neural-Behavioral Correlation Test
#' @description Correlational test for global or local network characteristics of neural network data with behavioral data
#' (still in testing phase)
#' @param bstat Behavioral statistic for each participant with neural data
#' @param nstat A statistic vector (whole-network) or matrix (ROI) from the \emph{neuralstat} function
#' @return Returns correlation test statistics for given statistics
#' @examples
#' \dontrun{
#' neuralarray <- convertConnBrainMat()
#' 
#' filteredneuralarray <- neuralnetfilter(neuralarray, method = "threshold", thres = .50)
#' 
#' AverageShortestPathLength <- neuralstat(filteredneuralarray, statistic = "ASPL")
#' 
#' Degree <- neuralstat(filteredneuralarray, statistic = "deg")
#' 
#' bstat <- iq
#' 
#' WholeNetwork_corr <- neuralcorrtest(bstat, AverageShortestPathLength)
#' 
#' ROI_corr <- neuralcorrtest(bstat, Degree)
#' }
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @references
#' Ly, A., Verhagen, A. J., & Wagenmakers, E.-J. (2016).
#' Harold Jeffreys's default Bayes factor hypothesis tests: Explanation, extension, and application in psychology.
#' \emph{Journal of Mathematical Psychology}, \emph{72}, 19-32.
#' 
#' Wagenmakers, E. J., Verhagen, J., & Ly, A. (2016).
#' How to quantify the evidence for the absence of a correlation.
#' \emph{Behavior Research Methods}, \emph{48}(2), 413-426.
#' @importFrom stats cor.test
#' @export
#Neural Correlation----
neuralcorrtest <- function (bstat, nstat)
{
    bstat<-as.matrix(bstat)
    nstat<-as.matrix(nstat)
    a<-length(nstat)/nrow(nstat)/ncol(nstat)
    
    if(a>1)
    {
        corr<-matrix(0,nrow=nrow(nstat),ncol=ncol(nstat))
        for(i in 1:nrow(nstat))
            for(j in 1:ncol(nstat))
                if(i!=j)
        corr[i,j]<-cor(bstat,nstat[i,j,])
                
                critical.r <- function(iter, a)
                {
                    df <- iter - 2
                    critical.t <- qt( a/2, df, lower.tail = F )
                    cvr <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
                    return(cvr)
                }
                
        sigp<-ifelse(corr>=critical.r(ncol(nstat),.05),corr,0)
        sign<-ifelse(corr<=-(critical.r(ncol(nstat),.05)),corr,0)
        
        colnames(sigp)<-colnames(nstat)
        row.names(sigp)<-colnames(nstat)
        colnames(sign)<-colnames(nstat)
        row.names(sign)<-colnames(nstat)
        
        corr<-list(pos=sigp,neg=sigp)
    }
    
    
    if(a==1)
    {
    corBF <- function (n, r)
    {
        bf10JeffreysIntegrate <- function(n, r, alpha=1) {
            # Jeffreys' test for whether a correlation is zero or not
            # Jeffreys (1961), pp. 289-292
            # This is the exact result, see EJ
            ##
            if ( any(is.na(r)) ){
                return(NaN)
            }
            
            # TODO: use which
            if (n > 2 && abs(r)==1) {
                return(Inf)
            }
            
            hyperTerm <- Re(hypergeo::hypergeo((2*n-3)/4, (2*n-1)/4, (n+2*alpha)/2, r^2))
            logTerm <- lgamma((n+2*alpha-1)/2)-lgamma((n+2*alpha)/2)-lbeta(alpha, alpha)
            myResult <- sqrt(pi)*2^(1-2*alpha)*exp(logTerm)*hyperTerm
            return(myResult)
        }
        
        
        # 3.0 One-sided preparation
        
        mPlusMarginalBJeffreys <- function(n, r, alpha=1){
            # Ly et al 2014
            # This is the exact result with symmetric beta prior on rho
            # This is the contribution of one-sided test
            #
            #	
            if ( any(is.na(r)) ){
                return(NaN)
            }
            if (n > 2 && r>=1) {
                return(Inf)
            } else if (n > 2 && r<=-1){
                return(0)
            }
            
            hyperTerm <- Re(hypergeo::genhypergeo(U=c(1, (2*n-1)/4, (2*n+1)/4),
                                                  L=c(3/2, (n+1+2*alpha)/2), z=r^2))
            logTerm <- -lbeta(alpha, alpha)
            myResult <- 2^(1-2*alpha)*r*(2*n-3)/(n+2*alpha-1)*exp(logTerm)*hyperTerm
            return(myResult)
        }
        
        
        bfPlus0JeffreysIntegrate <- function(n, r, alpha=1){
            # Ly et al 2014
            # This is the exact result with symmetric beta prior on rho
            #	
            if ( any(is.na(r)) ){
                return(NaN)
            }
            if (n > 2 && r>=1) {
                return(Inf)
            } else if (n > 2 && r<=-1){
                return(0)
            }
            
            bf10 <- bf10JeffreysIntegrate(n, r, alpha)
            mPlus <- mPlusMarginalBJeffreys(n, r, alpha)
            
            if (is.na(bf10) || is.na(mPlus)){
                return(NA)
            }
            
            myResult <- bf10+mPlus	
            return(myResult)
        }
        
        bf<-bfPlus0JeffreysIntegrate(n,r)
        
        return(bf)
        
    }
    
    if(ncol(nstat)==1)
    {corr<-list(corr=cor.test(bstat,nstat),corrBF=corBF(bstat,nstat))
    }else if(ncol(nstat)>1)
    {
        corr<-matrix(0,nrow=nrow(nstat),ncol=4)
        
        for(i in 1:nrow(nstat))
        {corr[i,1]<-round(cor.test(bstat,nstat[i,])$estimate,5)
        corr[i,2]<-round(cor.test(bstat,nstat[i,])$p.value,5)}
        
        corr[,3]<-fdrtool::fdrtool(as.numeric(corr[,2]),statistic = "pvalue",plot=FALSE,verbose=FALSE)$qval
        
        for(i in 1:nrow(nstat))
        {corr[i,4]<-corBF(nrow(bstat),corr[i,1])}
        
        row.names(corr)<-row.names(nstat)
        colnames(corr)<-c("r","p-value","FDR","BF")
        corr<-corr[which(corr[,2]<=.05),]
    }
    }
    
    return(noquote(corr))
}
#----
#' Connectome-based Predictive Modeling--Internal Validation
#' @description Applies the Connectome-based Predictive Modeling approach to neural data.
#' This method predicts a behavioral statistic using neural connectivity from the sample.
#' \strong{Please cite Finn et al., 2015; Rosenberg et al., 2016; Shen et al., 2017}
#' @param neuralarray Array from \emph{convertConnBrainMat} function
#' @param bstat Behavioral statistic for each participant with neural data (a vector)
#' @param covar Covariates to be included in predicting relevant edges (\strong{time consuming}).
#' \strong{Must} be input as a list() (see examples)
#' @param thresh Sets an \strong{alpha} threshold for edge weights to be retained.
#' Defaults to .01
#' @param method Use "mean" or "sum" of edge strengths in the positive and negative connectomes.
#' Defaults to "mean"
#' @param model Regression model to use for fitting the data.
#' Defaults to "linear"
#' @param corr Correlation method for assessing the relatonship between the behavioral measure and edges between ROIs.
#' Defaults to "pearson".
#' Set to "spearman" for non-linear or monotonic associations
#' @param shen Are ROIs from Shen et al. 2013 atlas? Defaults to FALSE.
#' Set to TRUE for canonical networks plot
#' @param cores Number of computer processing cores to use when performing covariate analyses.
#' Defaults to \emph{n} - 1 total number of cores.
#' Set to any number between 1 and maxmimum amount of cores on your computer
#' @param progBar Should progress bar be displayed?
#' Defaults to TRUE.
#' Set to FALSE for no progress bar
#' @return Returns a list containing a matrix (correlation coefficient (r), p-value (p), Bayes Factor (BF), mean absolute error (mae), root mean square error (rmse)).
#' The list also contains the positive (posMask) and negative (negMask) masks, which can be visualized here: \url{http://bisweb.yale.edu/build/connviewer.html}
#' @references 
#' Finn, E. S., Shen, X., Scheinost, D., Rosenberg, M. D., Huang, J., Chun, M. M., Papademetris, X., Constable, R. T. (2015).
#' Functional connectome fingerprinting: Identifying individuals using patterns of brain connectivity.
#' \emph{Nature Neuroscience}, \emph{18}(11), 1664-1671.
#' 
#' Ly, A., Verhagen, A. J., & Wagenmakers, E.-J. (2016).
#' Harold Jeffreys's default Bayes factor hypothesis tests: Explanation, extension, and application in psychology.
#' \emph{Journal of Mathematical Psychology}, \emph{72}, 19-32.
#' 
#' Rosenberg, M. D., Finn, E. S., Scheinost, D., Papademetris, X., Shen, X., Constable, R. T., Chun, M. M. (2016).
#' A neuromarker of sustained attention from whole-brain functional connectivity.
#' \emph{Nature Neuroscience}, \emph{19}(1), 165-171.
#'
#' Shen, X. Finn, E. S., Scheinost, D., Rosenberg, M. D., Chun, M. M., Papademetris, X., Constable, R. T. (2017).
#' Using connectome-based predictive modeling to predict individual behavior from brain connectivity.
#' \emph{Nature Protocols}, \emph{12}(3), 506-518.
#' 
#' Wagenmakers, E. J., Verhagen, J., & Ly, A. (2016).
#' How to quantify the evidence for the absence of a correlation.
#' \emph{Behavior Research Methods}, \emph{48}(2), 413-426.
#' 
#' Wei, T. & Simko, V.(2017).
#' R package "corrplot": Visualization of a correlation matrix (Version 0.84).
#' Available from \url{https://github.com/taiyun/corrplot}
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @importFrom graphics par
#' @importFrom grDevices colorRampPalette dev.new
#' @importFrom foreach %dopar%
#' @export
#CPM Internal Validation----
cpmIV <- function (neuralarray, bstat, covar, thresh = .01, method = c("mean", "sum"),
                    model = c("linear","quadratic","cubic"),
                   corr = c("pearson","spearman"), shen = FALSE, cores, progBar = TRUE)
{
    if(missing(method))
    {method<-"mean"
    }else{method<-match.arg(method)}
    
    if(missing(model))
    {model<-"linear"
    }else{model<-match.arg(model)}
    
    if(missing(corr))
    {corr<-"pearson"
    }else{corr<-match.arg(corr)}
    
    if(is.list(neuralarray))
    {neuralarray<-neuralarray[[1]]}
    
    if(missing(covar))
    {covar<-NULL
    }else if(!is.list(covar))
    {stop("Covariates vectors must be input as a list: list()")}
    
    #functions list
    critical.r <- function(iter, a)
    {
        df <- iter - 2
        critical.t <- qt( a/2, df, lower.tail = F )
        cvr <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
        return(cvr)
    }
    
    bstat<-scale(bstat)
    bstat<-as.vector(bstat)
    
    #number of subjects
    no_sub<-length(neuralarray)/nrow(neuralarray)/ncol(neuralarray)
    #number of nodes
    no_node<-ncol(neuralarray)
    
    #initialize positive and negative behavior stats
    behav_pred_pos<-matrix(0,nrow=no_sub,ncol=1)
    behav_pred_neg<-matrix(0,nrow=no_sub,ncol=1)
    
    if(is.list(covar))
    {
        cvars<-do.call(cbind,covar,1)
        cvars<-scale(cvars)
    }
    
    pos_array <- array(0,dim=c(nrow=no_node,ncol=no_node,no_sub))
    neg_array <- array(0,dim=c(nrow=no_node,ncol=no_node,no_sub))
    
    
    #perform leave-out analysis
    if(progBar)
    {pb <- txtProgressBar(max=no_sub, style = 3)}
    
    for(leftout in 1:no_sub)
    {
        train_mats<-neuralarray
        train_mats<-train_mats[,,-leftout]
        ##initialize train vectors
        #vector length
        vctrow<-ncol(neuralarray)^2
        vctcol<-length(train_mats)/nrow(train_mats)/ncol(train_mats)
        train_vcts<-matrix(0,nrow=vctrow,ncol=vctcol)
        for(i in 1:vctcol)
        {train_vcts[,i]<-as.vector(train_mats[,,i])}
        
        #behavior stats
        train_behav<-bstat
        train_behav<-train_behav[-leftout]
        
        #correlate edges with behavior
        if(nrow(train_vcts)!=(no_sub-1))
        {train_vcts<-t(train_vcts)}
        
        rmat<-vector(mode="numeric",length=ncol(train_vcts))
        pmat<-vector(mode="numeric",length=ncol(train_vcts))
        
        if(is.list(covar))
        {
            if(missing(cores))
            {
                cores <- parallel::detectCores()
                cores <- cores - 1
            }else{cores <- cores}#not to overload your computer
            cl <- parallel::makeCluster(cores)
            doSNOW::registerDoSNOW(cl)
            
            pcorr<-suppressWarnings(
                foreach::foreach(i=1:ncol(train_vcts))%dopar%
            {
                temp<-cbind(train_vcts[,i],train_behav,cvars[-leftout,])
                ppcor::pcor.test(temp[,1],temp[,2],temp[,c(seq(from=3,to=2+ncol(cvars)))])
            }
            )
            parallel::stopCluster(cl)
            
            for(i in 1:length(pcorr))
            {
                rmat[i]<-pcorr[[i]]$estimate
                pmat[i]<-pcorr[[i]]$p.value
            }
            rmat<-ifelse(is.na(rmat),0,rmat)
            pmat<-ifelse(is.na(pmat),0,pmat)
        }else{rmat<-suppressWarnings(cor(train_vcts,train_behav,method=corr))}
        
        r_mat<-matrix(rmat,nrow=no_node,ncol=no_node)
        
        #set threshold and define masks
        pos_mask<-matrix(0,nrow=no_node,ncol=no_node)
        neg_mask<-matrix(0,nrow=no_node,ncol=no_node)
        
        if(!is.list(covar))
        {
            #critical r-value
            cvr<-critical.r((no_sub-1),thresh)
            pos_edges<-which(r_mat>=cvr)
            neg_edges<-which(r_mat<=(-cvr))
        }else
        {
            p_mat<-matrix(pmat,nrow=no_node,ncol=no_node)
            sig<-ifelse(p_mat<=thresh,r_mat,0)
            pos_edges<-which(r_mat>0&sig!=0)
            neg_edges<-which(r_mat<0&sig!=0)
        }
        
        
        pos_mask[pos_edges]<-1
        neg_mask[neg_edges]<-1
        
        pos_array[,,leftout] <- pos_mask
        neg_array[,,leftout] <- neg_mask
        
        #get sum of all edges in TRAIN subs (divide, if symmetric matrices)
        train_sumpos<-matrix(0,nrow=(no_sub-1),ncol=1)
        train_sumneg<-matrix(0,nrow=(no_sub-1),ncol=1)
        
        for(ss in 1:nrow(train_sumpos))
        {
            if(method=="sum")
            {
                train_sumpos[ss]<-sum(train_mats[,,ss]*pos_mask)/2
                train_sumneg[ss]<-sum(train_mats[,,ss]*neg_mask)/2
            }else if(method=="mean")
            {
                train_sumpos[ss]<-mean(train_mats[,,ss]*pos_mask)/2
                train_sumneg[ss]<-mean(train_mats[,,ss]*neg_mask)/2
            }
        }
        
        #generate regression formula with covariates
        #if(is.list(covar))
        #{cvar<-cvars[-leftout,]}
        
        #regressions----
        
        #build model on TRAIN subs
        if(model=="linear")
        {
            fit_pos<-coef(lm(train_behav~train_sumpos))
            fit_neg<-coef(lm(train_behav~train_sumneg))
        }else if(model=="quadratic")
        {
            quad_pos<-train_sumpos^2
            quad_neg<-train_sumneg^2
            
            fit_pos<-coef(lm(train_behav~train_sumpos+quad_pos))
            fit_neg<-coef(lm(train_behav~train_sumneg+quad_neg))
        }else if(model=="cubic")
        {
            cube_pos<-train_sumpos^3
            cube_neg<-train_sumneg^3
            
            quad_pos<-train_sumpos^2
            quad_neg<-train_sumneg^2
            
            fit_pos<-coef(lm(train_behav~train_sumpos+quad_pos+cube_pos))
            fit_neg<-coef(lm(train_behav~train_sumneg+quad_neg+cube_neg))
        }
        
        #run model on TEST sub
        test_mat<-neuralarray[,,leftout]
        if(method=="sum")
        {
            test_sumpos<-sum(test_mat*pos_mask)/2
            test_sumneg<-sum(test_mat*neg_mask)/2
        }else if(method=="mean")
        {
            test_sumpos<-mean(test_mat*pos_mask)/2
            test_sumneg<-mean(test_mat*neg_mask)/2
        }
        
        if(model=="linear")
        {
            behav_pred_pos[leftout]<-fit_pos[2]*test_sumpos+fit_pos[1]
            behav_pred_neg[leftout]<-fit_neg[2]*test_sumneg+fit_neg[1]
        }else if(model=="quadratic")
        {
            quad_post<-test_sumpos^2
            quad_negt<-test_sumneg^2
            
            behav_pred_pos[leftout]<-fit_pos[3]*quad_post+fit_pos[2]*test_sumpos+fit_pos[1]
            behav_pred_neg[leftout]<-fit_neg[3]*quad_negt+fit_neg[2]*test_sumneg+fit_neg[1]
        }else if(model=="cubic")
        {
            cube_post<-test_sumpos^3
            cube_negt<-test_sumneg^3
            
            quad_post<-test_sumpos^2
            quad_negt<-test_sumneg^2
            
            behav_pred_pos[leftout]<-fit_pos[4]*cube_post+fit_pos[3]*quad_post+fit_pos[2]*test_sumpos+fit_pos[1]
            behav_pred_neg[leftout]<-fit_neg[4]*cube_negt+fit_neg[3]*quad_negt+fit_neg[2]*test_sumneg+fit_neg[1]
        }
        
        if(progBar)
        {setTxtProgressBar(pb, leftout)}
    }
    if(progBar)
    {close(pb)}
    
    pos_mat <- matrix(0, nrow = no_node, ncol = no_node)
    neg_mat <- matrix(0, nrow = no_node, ncol = no_node)
    
    for(i in 1:no_node)
        for(j in 1:no_node)
    {
        pos_mat[i,j] <- sum(pos_array[i,j,])
        neg_mat[i,j] <- sum(neg_array[i,j,])
    }
    
    posmask <- ifelse(pos_mat==no_sub,1,0)
    negmask <- ifelse(neg_mat==no_sub,1,0)
    
    
    R_pos<-cor(behav_pred_pos,bstat)
    P_pos<-cor.test(behav_pred_pos,bstat)$p.value
    R_neg<-cor(behav_pred_neg,bstat)
    P_neg<-cor.test(behav_pred_neg,bstat)$p.value
    
    P_pos<-ifelse(round(P_pos,3)!=0,round(P_pos,3),noquote("< .001"))
    P_neg<-ifelse(round(P_neg,3)!=0,round(P_neg,3),noquote("< .001"))
    
    
    #plot positive
    dev.new()
    par(mar=c(5,5,4,2))
    plot(bstat,behav_pred_pos,xlab="Observed Score\n(Z-score)",ylab="Predicted Score\n(Z-score)",
         main="Positive Prediction",xlim=c(-3,3),ylim=c(-3,3),pch=16,col="darkorange2")
    abline(lm(behav_pred_pos~bstat))
    if(R_pos>=0)
    {text(x=-2,y=2,
          labels = paste("r = ",round(R_pos,3),"\np = ",P_pos))
    }else if(R_pos<0)
    {text(x=-2,y=-2,
          labels = paste("r = ",round(R_pos,3),"\np = ",P_pos))}
    #plot negative
    dev.new()
    par(mar=c(5,5,4,2))
    plot(bstat,behav_pred_neg,xlab="Observed Score\n(Z-score)",ylab="Predicted Score\n(Z-score)",
         main="Negative Prediction",xlim=c(-3,3),ylim=c(-3,3),pch=16,col="skyblue2")
    abline(lm(behav_pred_neg~bstat))
    if(R_neg>=0)
    {text(x=-2,y=2,
          labels = paste("r = ",round(R_neg,3),"\np = ",P_neg))
    }else if(R_neg<0)
    {text(x=-2,y=-2,
          labels = paste("r = ",round(R_neg,3),"\np = ",P_neg))}
    
    #shen plots----
    if(shen==TRUE)
    {
        shennets<-c(2,4,3,2,3,3,2,2,2,1,4,1,3,2,4,1,2,4,2,4,2,
                2,5,5,5,5,5,4,4,2,2,4,5,5,5,4,5,5,5,5,8,6,
                8,4,5,5,2,2,3,3,5,1,1,1,2,1,1,5,8,5,5,5,5,
                1,1,8,8,6,8,2,8,6,8,8,6,7,6,7,6,6,7,6,4,5,
                3,3,6,4,5,3,4,5,4,4,4,3,5,6,4,7,4,7,4,4,4,
                4,4,4,5,4,2,2,4,4,3,2,4,4,4,4,4,4,4,4,4,4,
                4,4,4,4,4,4,4,3,4,4,1,3,2,1,3,2,2,4,1,4,2,
                1,1,1,1,4,1,2,4,1,2,5,5,5,5,1,5,2,1,5,5,5,
                4,5,5,5,5,5,8,6,8,4,5,5,5,2,1,2,1,1,1,5,5,
                1,5,1,2,1,5,2,5,6,2,8,8,5,3,8,6,8,6,6,8,8,
                6,7,7,7,6,6,4,5,1,4,4,3,3,4,3,4,3,5,4,4,4,
                4,4,4,5,4,4,4,3,8,7,2,4,4,4,2,2,4,4,4,4,4,
                4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)
        
        #lobes
        rtlobenets<-c(rep("PFC",22),rep("Mot",11),rep("Ins",4),rep("Par",13),
                    rep("Tem",21),rep("Occ",11),rep("Lim",17),rep("Cer",20),
                    rep("Sub",9),rep("Bsm",5))
        
        ltlobenets<-c(rep("PFC",24),rep("Mot",10),rep("Ins",3),rep("Par",14),
                     rep("Tem",18),rep("Occ",14),rep("Lim",19),rep("Cer",21),
                     rep("Sub",8),rep("Bsm",4))
        
        lobenets<-c(rtlobenets,ltlobenets)
        
        pos_lobe<-posmask
        neg_lobe<-negmask
        colnames(pos_lobe)<-lobenets
        colnames(neg_lobe)<-lobenets
        
        ten<-c("PFC","Mot","Ins","Par","Tem","Occ","Lim","Cer","Sub","Bsm")
        
        colnames(pos_lobe)<-lobenets
        colnames(neg_lobe)<-lobenets
        
        poslobemat<-matrix(0,nrow=10,ncol=10)
        neglobemat<-matrix(0,nrow=10,ncol=10)
        
        for(i in 1:10)
            for(j in 1:10)
            {
                poslobemat[i,j]<-sum(pos_lobe[which(colnames(pos_lobe)==ten[i]),which(colnames(pos_lobe)==ten[j])])
                neglobemat[i,j]<-sum(neg_lobe[which(colnames(neg_lobe)==ten[i]),which(colnames(neg_lobe)==ten[j])])
            }
        
        ldiffmat<-(poslobemat-neglobemat)
        
        colnames(ldiffmat)<-c("PFC","Mot","Ins","Par","Tem","Occ","Lim","Cer","Sub","Bsm")
        row.names(ldiffmat)<-c("PFC","Mot","Ins","Par","Tem","Occ","Lim","Cer","Sub","Bsm")
        
        ldiffmat[upper.tri(ldiffmat)]<-0
        
        llim<-ifelse(abs(min(ldiffmat))>max(ldiffmat),abs(min(ldiffmat)),max(ldiffmat))
        
        colo<-colorRampPalette(c("skyblue2","white","darkorange2"))
        
        dev.new()
        corrplot::corrplot(ldiffmat,is.corr=FALSE,method="color",
                           tl.col="black",col = colo(100),na.label="square",
                           na.label.col = "white",addgrid.col="black",
                           title="Difference in the Number of Edges\nin Macroscale Regions",
                           mar=c(0,0,4,0),cl.length=3,cl.pos="b",cl.lim=c(-llim,llim))
    
        #canonical networks
        pos_nets<-posmask
        neg_nets<-negmask
        
        colnames(pos_nets)<-shennets
        colnames(neg_nets)<-shennets
    
        posnetmat<-matrix(0,nrow=max(shennets),ncol=max(shennets))
        negnetmat<-matrix(0,nrow=max(shennets),ncol=max(shennets))
    
        for(i in 1:max(shennets))
        for(j in 1:max(shennets))
            {
                posnetmat[i,j]<-sum(pos_nets[which(colnames(pos_nets)==i),which(colnames(pos_nets)==j)])
                negnetmat[i,j]<-sum(neg_nets[which(colnames(neg_nets)==i),which(colnames(neg_nets)==j)])
            }
    
        diffmat<-(posnetmat-negnetmat)
        
        colnames(diffmat)<-c("MF","FP","DM","SubC","MT","VI","VII","VA")
        row.names(diffmat)<-c("MF","FP","DM","SubC","MT","VI","VII","VA")

        diffmat[upper.tri(diffmat)]<-0
        
        dlim<-ifelse(abs(min(diffmat))>max(diffmat),abs(min(diffmat)),max(diffmat))
    
        colo<-colorRampPalette(c("skyblue2","white","darkorange2"))
        
        dev.new()
        corrplot::corrplot(diffmat,is.corr=FALSE,method="color",
                       tl.col="black",col = colo(100),na.label="square",
                       na.label.col = "white",addgrid.col="black",
                       title="Difference in the Number of Edges\nin the Canonical Networks",
                       mar=c(0,0,4,0),cl.length=3,cl.pos="b",cl.lim=c(-dlim,dlim))
    
    }
    
    bstat<-as.vector(bstat)
    behav_pred_pos<-as.vector(behav_pred_pos)
    behav_pred_neg<-as.vector(behav_pred_neg)
    
    corBF <- function (n, r)
    {
        bf10JeffreysIntegrate <- function(n, r, alpha=1) {
            # Jeffreys' test for whether a correlation is zero or not
            # Jeffreys (1961), pp. 289-292
            # This is the exact result, see EJ
            ##
            if ( any(is.na(r)) ){
                return(NaN)
            }
            
            # TODO: use which
            if (n > 2 && abs(r)==1) {
                return(Inf)
            }
            
            hyperTerm <- Re(hypergeo::hypergeo((2*n-3)/4, (2*n-1)/4, (n+2*alpha)/2, r^2))
            logTerm <- lgamma((n+2*alpha-1)/2)-lgamma((n+2*alpha)/2)-lbeta(alpha, alpha)
            myResult <- sqrt(pi)*2^(1-2*alpha)*exp(logTerm)*hyperTerm
            return(myResult)
        }
        
        
        # 3.0 One-sided preparation
        
        mPlusMarginalBJeffreys <- function(n, r, alpha=1){
            # Ly et al 2014
            # This is the exact result with symmetric beta prior on rho
            # This is the contribution of one-sided test
            #
            #	
            if ( any(is.na(r)) ){
                return(NaN)
            }
            if (n > 2 && r>=1) {
                return(Inf)
            } else if (n > 2 && r<=-1){
                return(0)
            }
            
            hyperTerm <- Re(hypergeo::genhypergeo(U=c(1, (2*n-1)/4, (2*n+1)/4),
                                                  L=c(3/2, (n+1+2*alpha)/2), z=r^2))
            logTerm <- -lbeta(alpha, alpha)
            myResult <- 2^(1-2*alpha)*r*(2*n-3)/(n+2*alpha-1)*exp(logTerm)*hyperTerm
            return(myResult)
        }
        
        
        bfPlus0JeffreysIntegrate <- function(n, r, alpha=1){
            # Ly et al 2014
            # This is the exact result with symmetric beta prior on rho
            #	
            if ( any(is.na(r)) ){
                return(NaN)
            }
            if (n > 2 && r>=1) {
                return(Inf)
            } else if (n > 2 && r<=-1){
                return(0)
            }
            
            bf10 <- bf10JeffreysIntegrate(n, r, alpha)
            mPlus <- mPlusMarginalBJeffreys(n, r, alpha)
            
            if (is.na(bf10) || is.na(mPlus)){
                return(NA)
            }
            
            myResult <- bf10+mPlus	
            return(myResult)
        }
        
        bf<-bfPlus0JeffreysIntegrate(n,r)
        
        return(bf)
        
    }
    
    pos_BF<-corBF(length(bstat),abs(R_pos))
    neg_BF<-corBF(length(bstat),abs(R_neg))
    
    for(i in 1:length(bstat))
    {
        #mae
        mae_pos<-sum(abs(behav_pred_pos[i]-bstat[i]))/length(bstat)
        mae_neg<-sum(abs(behav_pred_neg[i]-bstat[i]))/length(bstat)
        
        #rmse
        pos_rmse<-sqrt(sum((behav_pred_pos[i]-bstat[i])^2)/length(bstat))
        neg_rmse<-sqrt(sum((behav_pred_neg[i]-bstat[i])^2)/length(bstat))
    }
    
    results<-matrix(0,nrow=2,ncol=5)
    
    results[1,1]<-round(R_pos,3)
    results[1,2]<-P_pos
    results[1,3]<-round(pos_BF,3)
    results[1,4]<-round(mae_pos,3)
    results[1,5]<-round(pos_rmse,3)
    results[2,1]<-round(R_neg,3)
    results[2,2]<-P_neg
    results[2,3]<-round(neg_BF,3)
    results[2,4]<-round(mae_neg,3)
    results[2,5]<-round(neg_rmse,3)
    
    colnames(results)<-c("r","p","BF","mae","rmse")
    row.names(results)<-c("positive","negative")
    
    return(list(results=results,posMask=posmask,negMask=negmask))
}
#----
#' Connectome-based Predictive Modeling--External Validation
#' @description Applies the Connectome-based Predictive Modeling approach to neural data.
#' This method predicts a behavioral statistic using neural connectivity from the sample.
#' Results \strong{may} differ from \emph{Matlab} results because of robust GLM methodology.
#' This function is still in its \strong{testing phase}.
#' \strong{Please cite Finn et al., 2015; Rosenberg et al., 2016; Shen et al., 2017}
#' @param train_na Training dataset
#' (an array from \emph{convertConnBrainMat} function)
#' @param train_b Behavioral statistic for each participant for the \strong{training} neural data (a vector)
#' @param valid_na Validation dataset
#' (an array from \emph{convertConnBrainMat} function)
#' @param valid_b Behavioral statistic for each participant for the \strong{validation} neural data (a vector)
#' @param thresh Sets an \strong{alpha} threshold for edge weights to be retained.
#' Defaults to .01
#' @param overlap Should leave-one-out cross-validation be used?
#' Defaults to FALSE (use full dataset, no leave-one-out).
#' Set to TRUE to select edges that appear in every leave-one-out cross-validation network (\emph{time consuming})
#' @param progBar Should progress bar be displayed?
#' Defaults to TRUE.
#' Set to FALSE for no progress bar
#' @return Returns a list containing a matrix (r coefficient (r), p-value (p-value), Bayes Factor (BF), mean absolute error (mae), root mean square error (rmse)).
#' The list also contains the positive (posMask) and negative (negMask) masks used
#' @references 
#' Finn, E. S., Shen, X., Scheinost, D., Rosenberg, M. D., Huang, J., Chun, M. M., Papademetris, X., Constable, R. T. (2015).
#' Functional connectome fingerprinting: Identifying individuals using patterns of brain connectivity.
#' \emph{Nature Neuroscience}, \emph{18}(11), 1664-1671.
#' 
#' Ly, A., Verhagen, A. J., & Wagenmakers, E.-J. (2016).
#' Harold Jeffreys's default Bayes factor hypothesis tests: Explanation, extension, and application in psychology.
#' \emph{Journal of Mathematical Psychology}, \emph{72}, 19-32.
#' 
#' Rosenberg, M. D., Finn, E. S., Scheinost, D., Papademetris, X., Shen, X., Constable, R. T., Chun, M. M. (2016).
#' A neuromarker of sustained attention from whole-brain functional connectivity.
#' \emph{Nature Neuroscience}, \emph{19}(1), 165-171.
#'
#' Shen, X. Finn, E. S., Scheinost, D., Rosenberg, M. D., Chun, M. M., Papademetris, X., Constable, R. T. (2017).
#' Using connectome-based predictive modeling to predict individual behavior from brain connectivity.
#' \emph{Nature Protocols}, \emph{12}(3), 506-518.
#' 
#' Wagenmakers, E. J., Verhagen, J., & Ly, A. (2016).
#' How to quantify the evidence for the absence of a correlation.
#' \emph{Behavior Research Methods}, \emph{48}(2), 413-426.
#' 
#' Wei, T. & Simko, V.(2017).
#' R package "corrplot": Visualization of a correlation matrix (Version 0.84).
#' Available from \url{https://github.com/taiyun/corrplot}
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @importFrom MASS psi.bisquare
#' @export
#CPM External Validation----
cpmEV <- function (train_na, train_b, valid_na, valid_b,
                   thresh = .01, overlap = FALSE, progBar = TRUE)
{
    #number of nodes
    n_node<-ncol(train_na)
    
    #training data
    n_sub<-length(train_na)/nrow(train_na)/ncol(train_na)
    n_train_sub<-n_sub-1
    
    #validation data
    n_validation_sub<-length(valid_na)/nrow(valid_na)/ncol(valid_na)
    
    aa<-matrix(1,nrow=n_node,ncol=n_node)
    aa[lower.tri(aa,diag = TRUE)]<-0
    upp_id<-which(aa==1)
    n_edge<-length(upp_id)
    
    train_b<-scale(train_b)
    valid_b<-scale(valid_b)
    
    critical.r <- function(iter, a)
    {
        df <- iter - 2
        critical.t <- qt( a/2, df, lower.tail = F )
        cvr <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
        return(cvr)
    }
    
    if(overlap==TRUE)
    {
        if(progBar)
        {pb <- txtProgressBar(max=n_sub, style = 3)}
        
        pos_mask_all<-array(0,dim=c(n_node,n_node,n_sub))
        neg_mask_all<-array(0,dim=c(n_node,n_node,n_sub))
        
        for(excl_sub in 1:n_sub)
        {
            #exclude data from left-out subject
            train_mats_tmp<-train_na
            train_mats_tmp<-train_mats_tmp[,,-excl_sub]
            train_behav<-train_b
            train_behav<-train_behav[-excl_sub]
            
            #create n_train_sub x n_edge matrix
            vctrow<-ncol(train_na)^2
            vctcol<-length(train_mats_tmp)/nrow(train_mats_tmp)/ncol(train_mats_tmp)
            train_vect<-matrix(0,nrow=vctrow,ncol=vctcol)
            for(i in 1:vctcol)
            {train_vect[,i]<-as.vector(train_mats_tmp[,,i])}
            train_vect<-t(train_vect)
            upp_vect<-train_vect[,upp_id]
            
            #relate behavior to edge strength
            cp<-matrix(0,nrow=n_edge,ncol=1)
            cr<-matrix(0,nrow=n_edge,ncol=1)
            
            for(ii in 1:n_edge)
            {
                j<-summary(MASS::rlm(train_behav~upp_vect[,ii],psi = psi.bisquare))
                cr[ii]<-sign(j[4]$coefficients[6])*sqrt(((j[4]$coefficients[6]^2)/(n_train_sub-2))/(1+(j[4]$coefficients[6]^2)/(n_train_sub-2)))
            }
            
            #select edges based on threshold
            pos_edge<-matrix(0,nrow=1,ncol=n_edge)
            neg_edge<-matrix(0,nrow=1,ncol=n_edge)
            
            cvr<-critical.r((n_sub-1),.01)
            
            cp_pos<-which(cr>=cvr)
            cp_neg<-which(cr<=(-cvr))
            pos_edge[cp_pos]<-1
            neg_edge[cp_neg]<-1
            
            pos_mask<-matrix(0,nrow=n_node,ncol=n_node)
            neg_mask<-matrix(0,nrow=n_node,ncol=n_node)
            
            pos_mask[upp_id]<-pos_edge
            pos_mask<-pos_mask+t(pos_mask)
            neg_mask[upp_id]<-neg_edge
            neg_mask<-neg_mask+t(neg_mask)
            
            pos_mask_all[,,excl_sub]<-pos_mask
            neg_mask_all[,,excl_sub]<-neg_mask
            
            if(progBar)
            {setTxtProgressBar(pb, excl_sub)}
        }
        
        if(progBar)
        {close(pb)}
        
        pos_overlap<-matrix(0,nrow=n_node,ncol=n_node)
        neg_overlap<-matrix(0,nrow=n_node,ncol=n_node)
        
        for(i in 1:n_node)
            for(j in 1:n_node)
            {
                pos_overlap[i,j]<-sum(pos_mask_all[i,j,])
                neg_overlap[i,j]<-sum(neg_mask_all[i,j,])
            }
        
        pos_overlap<-ifelse(pos_overlap==n_sub,1,0)
        neg_overlap<-ifelse(neg_overlap==n_sub,1,0)
        
        
    }else if(overlap==FALSE)
    {
        if(progBar)
        {pb <- txtProgressBar(max=n_edge, style = 3)}
        
        #create n_train_sub x n_edge matrix
        vctrow<-ncol(train_na)^2
        vctcol<-length(train_mats_tmp)/nrow(train_mats_tmp)/ncol(train_mats_tmp)
        train_vect<-matrix(0,nrow=vctrow,ncol=vctcol)
        for(i in 1:vctcol)
        {train_vect[,i]<-as.vector(train_mats_tmp[,,i])}
        train_vect<-t(train_vect)
        upp_vect<-train_vect[,upp_id]
        
        #relate behavior to edge strength
        cr<-matrix(0,nrow=n_edge,ncol=1)
        
        for(ii in 1:n_edge)
        {
            j<-summary(MASS::rlm(train_behav~upp_vect[,ii],psi = psi.bisquare))
            b<-j[4]$coefficients[1:2]
            cr[ii]<-sign(j[4]$coefficients[6])*sqrt(((j[4]$coefficients[6]^2)/(n_train_sub-2))/(1+(j[4]$coefficients[6]^2)/(n_train_sub-2)))
            
            if(progBar)
            {setTxtProgressBar(pb, ii)}
        }
        
        if(progBar)
        {close(pb)}
        
        #select edges based on threshold
        pos_edge<-matrix(0,nrow=1,ncol=n_edge)
        neg_edge<-matrix(0,nrow=1,ncol=n_edge)
        
        cvr<-critical.r((n_sub-1),.01)
        
        cp_pos<-which(cr>=cvr)
        cp_neg<-which(cr<=(-cvr))
        pos_edge[cp_pos]<-1
        neg_edge[cp_neg]<-1
        
        pos_mask<-matrix(0,nrow=n_node,ncol=n_node)
        neg_mask<-matrix(0,nrow=n_node,ncol=n_node)
        
        pos_mask[upp_id]<-pos_edge
        pos_mask<-pos_mask+t(pos_mask)
        neg_mask[upp_id]<-neg_edge
        neg_mask<-neg_mask+t(neg_mask)
        
        pos_overlap<-pos_mask
        neg_overlap<-neg_mask
    }
    
    #sum edges for all subjects in the training set
    train_pos_sum<-matrix(0,nrow=n_sub,ncol=1)
    train_neg_sum<-matrix(0,nrow=n_sub,ncol=1)
    
    for(k in 1:n_sub)
    {
        train_pos_sum[k]<-sum(pos_overlap*train_na[,,k])
        train_neg_sum[k]<-sum(neg_overlap*train_na[,,k])
    }
    
    #build model with training data
    b_pos<-MASS::rlm(train_b~train_pos_sum,psi = psi.bisquare)
    b_neg<-MASS::rlm(train_b~train_neg_sum,psi = psi.bisquare)
    robGLM_fit<-MASS::rlm(train_b~train_pos_sum+train_neg_sum,psi = psi.bisquare)
    
    b_posc<-b_pos$coefficients
    b_negc<-b_neg$coefficients
    robGLM_fitc<-robGLM_fit$coefficients
    
    #generate predictions for validation set
    pred_pos<-matrix(0,nrow=n_validation_sub,ncol=1)
    pred_neg<-matrix(0,nrow=n_validation_sub,ncol=1)
    pred_glm<-matrix(0,nrow=n_validation_sub,ncol=1)
    
    validation_pos_sum<-matrix(0,nrow=n_validation_sub,ncol=1)
    validation_neg_sum<-matrix(0,nrow=n_validation_sub,ncol=1)
    
    for(vs in 1:n_validation_sub)
    {
        validation_pos_sum[vs]<-sum(pos_overlap*valid_na[,,vs])
        validation_neg_sum[vs]<-sum(neg_overlap*valid_na[,,vs])
        
        pred_pos[vs]<-(b_posc[2]*validation_pos_sum[vs])+b_posc[1]
        pred_neg[vs]<-(b_negc[2]*validation_neg_sum[vs])+b_negc[1]
        pred_glm[vs]<-robGLM_fitc[1]+(robGLM_fitc[2]*validation_pos_sum[vs])+(robGLM_fitc[3]*validation_neg_sum[vs])
    }
    
    corBF <- function (n, r)
    {
        bf10JeffreysIntegrate <- function(n, r, alpha=1) {
            # Jeffreys' test for whether a correlation is zero or not
            # Jeffreys (1961), pp. 289-292
            # This is the exact result, see EJ
            ##
            if ( any(is.na(r)) ){
                return(NaN)
            }
            
            # TODO: use which
            if (n > 2 && abs(r)==1) {
                return(Inf)
            }
            
            hyperTerm <- Re(hypergeo::hypergeo((2*n-3)/4, (2*n-1)/4, (n+2*alpha)/2, r^2))
            logTerm <- lgamma((n+2*alpha-1)/2)-lgamma((n+2*alpha)/2)-lbeta(alpha, alpha)
            myResult <- sqrt(pi)*2^(1-2*alpha)*exp(logTerm)*hyperTerm
            return(myResult)
        }
        
        
        # 3.0 One-sided preparation
        
        mPlusMarginalBJeffreys <- function(n, r, alpha=1){
            # Ly et al 2014
            # This is the exact result with symmetric beta prior on rho
            # This is the contribution of one-sided test
            #
            #	
            if ( any(is.na(r)) ){
                return(NaN)
            }
            if (n > 2 && r>=1) {
                return(Inf)
            } else if (n > 2 && r<=-1){
                return(0)
            }
            
            hyperTerm <- Re(hypergeo::genhypergeo(U=c(1, (2*n-1)/4, (2*n+1)/4),
                                                  L=c(3/2, (n+1+2*alpha)/2), z=r^2))
            logTerm <- -lbeta(alpha, alpha)
            myResult <- 2^(1-2*alpha)*r*(2*n-3)/(n+2*alpha-1)*exp(logTerm)*hyperTerm
            return(myResult)
        }
        
        
        bfPlus0JeffreysIntegrate <- function(n, r, alpha=1){
            # Ly et al 2014
            # This is the exact result with symmetric beta prior on rho
            #	
            if ( any(is.na(r)) ){
                return(NaN)
            }
            if (n > 2 && r>=1) {
                return(Inf)
            } else if (n > 2 && r<=-1){
                return(0)
            }
            
            bf10 <- bf10JeffreysIntegrate(n, r, alpha)
            mPlus <- mPlusMarginalBJeffreys(n, r, alpha)
            
            if (is.na(bf10) || is.na(mPlus)){
                return(NA)
            }
            
            myResult <- bf10+mPlus	
            return(myResult)
        }
        
        bf<-bfPlus0JeffreysIntegrate(n,r)
        
        return(bf)
        
    }
    
    BF_pos<-corBF(length(valid_b),r_pos)
    BF_neg<-corBF(length(valid_b),r_neg)
    BF_glm<-corBF(length(valid_b),r_glm)
    
    for(i in 1:length(valid_b))
    {
        #mae
        mae_pos<-sum(abs(pred_pos[i]-valid_b[i]))/length(valid_b)
        mae_neg<-sum(abs(pred_neg[i]-valid_b[i]))/length(valid_b)
        mae_glm<-sum(abs(pred_glm[i]-valid_b[i]))/length(valid_b)
        
        #rmse
        rmse_pos<-sqrt(sum((pred_pos[i]-valid_b[i])^2)/length(valid_b))
        rmse_neg<-sqrt(sum((pred_neg[i]-valid_b[i])^2)/length(valid_b))
        rmse_glm<-sqrt(sum((pred_glm[i]-valid_b[i])^2)/length(valid_b))
    }
    
    r_pos<-cor(valid_b,pred_pos)
    p_pos<-cor.test(valid_b,pred_pos)$p.value
    r_neg<-cor(valid_b,pred_neg)
    p_neg<-cor.test(valid_b,pred_neg)$p.value
    r_glm<-cor(valid_b,pred_glm)
    p_glm<-cor.test(valid_b,pred_glm)$p.value
    
    p_pos<-ifelse(round(p_pos,3)!=0,round(p_pos,3),noquote("<.001"))
    p_neg<-ifelse(round(p_neg,3)!=0,round(p_neg,3),noquote("<.001"))
    p_glm<-ifelse(round(p_glm,3)!=0,round(p_glm,3),noquote("<.001"))
    
    results<-matrix(0,nrow=3,ncol=5)
    
    results[1,1]<-round(r_pos,3)
    results[1,2]<-p_pos
    results[1,3]<-round(BF_pos,3)
    results[1,4]<-round(mae_pos,3)
    results[1,5]<-round(rmse_pos,3)
    results[2,1]<-round(r_neg,3)
    results[2,2]<-p_neg
    results[2,3]<-round(BF_neg,3)
    results[2,4]<-round(mae_neg,3)
    results[2,5]<-round(rmse_neg,3)
    results[3,1]<-round(r_glm,3)
    results[3,2]<-p_glm
    results[3,3]<-round(BF_glm,3)
    results[3,4]<-round(mae_glm,3)
    results[3,5]<-round(rmse_pos,3)
    
    colnames(results)<-c("r","p-value","BF","mae","rmse")
    row.names(results)<-c("positive","negative","full")
    
    #plot positive
    dev.new()
    par(mar=c(5,5,4,2))
    plot(valid_b,pred_pos,xlab="Observed Score\n(Z-score)",ylab="Predicted Score\n(Z-score)",
         main="Positive Prediction",xlim=c(-3,3),ylim=c(-3,3),pch=16,col="darkorange2")
    abline(lm(pred_pos~valid_b))
    if(r_pos>=0)
    {text(x=-2,y=2,
          labels = paste("r = ",round(r_pos,3),"\np = ",p_pos))
    }else if(r_pos<0)
    {text(x=-2,y=-2,
          labels = paste("r = ",round(r_pos,3),"\np = ",p_pos))}
    #plot negative
    dev.new()
    par(mar=c(5,5,4,2))
    plot(valid_b,pred_neg,xlab="Observed Score\n(Z-score)",ylab="Predicted Score\n(Z-score)",
         main="Negative Prediction",xlim=c(-3,3),ylim=c(-3,3),pch=16,col="skyblue2")
    abline(lm(pred_neg~valid_b))
    if(r_neg>=0)
    {text(x=-2,y=2,
          labels = paste("r = ",round(r_neg,3),"\np = ",p_neg))
    }else if(r_neg<0)
    {text(x=-2,y=-2,
          labels = paste("r = ",round(r_neg,3),"\np = ",p_neg))}
    #plot full
    dev.new()
    par(mar=c(5,5,4,2))
    plot(valid_b,pred_glm,xlab="Observed Score\n(Z-score)",ylab="Predicted Score\n(Z-score)",
         main="Full Prediction",xlim=c(-3,3),ylim=c(-3,3),pch=16,col="darkolivegreen2")
    abline(lm(pred_glm~valid_b))
    if(r_glm>=0)
    {text(x=-2,y=2,
          labels = paste("r = ",round(r_glm,3),"\np = ",p_glm))
    }else if(r_glm<0)
    {text(x=-2,y=-2,
          labels = paste("r = ",round(r_glm,3),"\np = ",p_glm))}
    
    return(list(results=results,posMask=pos_mask,negMask=neg_mask))
}
#----
#' Connectome-based Predictive Modeling--Fingerprinting
#' @description Applies the Connectome-based Predictive Modeling approach to neural data.
#' This method identifies individuals based on their specific connectivity patterns.
#' \strong{Please cite Finn et al., 2015; Rosenberg et al., 2016; Shen et al., 2017}
#' @param session1 Array from \emph{convertConnBrainMat} function
#' (first session)
#' @param session2 Array from \emph{convertConnBrainMat} function
#' (second session)
#' @param progBar Should progress bar be displayed?
#' Defaults to TRUE.
#' Set to FALSE for no progress bar
#' @return Returns a matrix containing the percentage and number of correctly identified subjects for sessions 1 and 2
#' @references 
#' Finn, E. S., Shen, X., Scheinost, D., Rosenberg, M. D., Huang, J., Chun, M. M., Papademetris, X., Constable, R. T. (2015).
#' Functional connectome fingerprinting: Identifying individuals using patterns of brain connectivity.
#' \emph{Nature Neuroscience}, \emph{18}(11), 1664-1671.
#' 
#' Rosenberg, M. D., Finn, E. S., Scheinost, D., Papademetris, X., Shen, X., Constable, R. T., Chun, M. M. (2016).
#' A neuromarker of sustained attention from whole-brain functional connectivity.
#' \emph{Nature Neuroscience}, \emph{19}(1), 165-171.
#'
#' Shen, X. Finn, E. S., Scheinost, D., Rosenberg, M. D., Chun, M. M., Papademetris, X., Constable, R. T. (2017).
#' Using connectome-based predictive modeling to predict individual behavior from brain connectivity.
#' \emph{Nature Protocols}, \emph{12}(3), 506-518.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#CPM Fingerprinting----
cpmFP <- function (session1, session2, progBar = TRUE)
{
    count1<-0
    count2<-0
    
    if(is.list(session1))
    {session1<-session1[[1]]}
    
    if(is.list(session2))
    {session2<-session2[[1]]}
    
    m<-nrow(session1)*ncol(session1)
    n<-length(session1)/m
    
    no_sub<-n
    
    if(isSymmetric(session1[,,1]))
    {sesh1<-matrix(as.vector(session1),nrow=m,ncol=n)}
    if(isSymmetric(session2[,,1]))
    {sesh2<-matrix(as.vector(session2),nrow=m,ncol=n)}
    
    tt_cor<-matrix(0,nrow=nrow(sesh1),ncol=1)
    
    if(progBar)
    {pb <- txtProgressBar(max=no_sub, style = 3)}
    
    for(i in 1:no_sub)
    {
        #session 1    
        tt_cor<-sesh2[,i]
        
        tt_to_all<-cor(tt_cor,sesh1)
        
        va_id<-which.max(tt_to_all)
        
        if(i == va_id)
        {count1<-count1+1}
        
        #session 2
        tt_cor<-sesh1[,i]
        
        tt_to_all<-cor(tt_cor,sesh2)
        
        va_id<-which.max(tt_to_all)
        
        if(i == va_id)
        {count2<-count2+1}
        
        if(progBar)
        {setTxtProgressBar(pb, i)}
    }
    
    if(progBar)
    {close(pb)}
    
    ident<-matrix(0,nrow=2,ncol=2)
    
    ident[1,1]<-round((count1/no_sub),3)
    ident[1,2]<-round((count2/no_sub),3)
    ident[2,1]<-round(count1,0)
    ident[2,2]<-round(count2,0)
    
    row.names(ident)<-c("percentage identified","number identified")
    colnames(ident)<-c("session 1","session 2")
    
    return(ident)
}
#----
#' Connectome-based Predictive Modeling--Fingerprinting Permutation
#' @description Applies the Connectome-based Predictive Modeling approach to neural data.
#' This method identifies individuals based on their specific connectivity patterns.
#' \strong{Please cite Finn et al., 2015; Rosenberg et al., 2016; Shen et al., 2017}
#' @param session1 Array from \emph{convertConnBrainMat} function
#' (first session)
#' @param session2 Array from \emph{convertConnBrainMat} function
#' (second session)
#' @param iter Number of iterations to perform. Defaults to 1,000
#' @param progBar Should progress bar be displayed?
#' Defaults to TRUE.
#' Set to FALSE for no progress bar
#' @return Returns a matrix containing the percentage and number of correctly identified subjects for sessions 1 and 2
#' @references 
#' Finn, E. S., Shen, X., Scheinost, D., Rosenberg, M. D., Huang, J., Chun, M. M., Papademetris, X., Constable, R. T. (2015).
#' Functional connectome fingerprinting: Identifying individuals using patterns of brain connectivity.
#' \emph{Nature Neuroscience}, \emph{18}(11), 1664-1671.
#' 
#' Rosenberg, M. D., Finn, E. S., Scheinost, D., Papademetris, X., Shen, X., Constable, R. T., Chun, M. M. (2016).
#' A neuromarker of sustained attention from whole-brain functional connectivity.
#' \emph{Nature Neuroscience}, \emph{19}(1), 165-171.
#'
#' Shen, X. Finn, E. S., Scheinost, D., Rosenberg, M. D., Chun, M. M., Papademetris, X., Constable, R. T. (2017).
#' Using connectome-based predictive modeling to predict individual behavior from brain connectivity.
#' \emph{Nature Protocols}, \emph{12}(3), 506-518.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#CPM Permutation Testing----
cpmFPperm <- function (session1, session2, iter = 1000, progBar = TRUE)
{
    rate<-matrix(nrow=iter,ncol=4)
    
    if(is.list(session1))
    {session1<-session1[[1]]}
    
    if(is.list(session2))
    {session2<-session2[[1]]}
    
    no_sub<-length(session1)/nrow(session1)/ncol(session1)
    
    m<-nrow(session1)*ncol(session1)
    n<-length(session1)/m
    
    if(isSymmetric(session1[,,1]))
    {sesh1<-matrix(as.vector(session1),nrow=m,ncol=n)}
    if(isSymmetric(session2[,,1]))
    {sesh2<-matrix(as.vector(session2),nrow=m,ncol=n)}
    
    if(progBar)
    {pb <- txtProgressBar(max=iter*no_sub, style = 3)}
    
    count3<-0
    
    for(j in 1:iter)
    {
        sub_order<-sample(no_sub)
        
        all_se1<-sesh1
        all_se2<-sesh2[,sub_order]
        
        count1<-0
        count2<-0
        
        for(i in 1:no_sub)
        {
            #session 1
            tt_cor<-sesh2[,i]
            
            tt_to_all<-cor(tt_cor,all_se1)
            
            va<-max(tt_to_all)
            va_id<-which.max(tt_to_all)
            
            if(i == va_id)
            {count1<-count1+1}
            
            #session 2
            tt_cor<-sesh1[,i]
            
            tt_to_all<-cor(tt_cor,sesh2)
            
            va<-max(tt_to_all)
            va_id<-which.max(tt_to_all)
            
            if(i == va_id)
            {count2<-count2+1}
            
            count3<-count3+1
            
            if(progBar)
            {setTxtProgressBar(pb, count3)}
        }
        
        rate[j,1]<-round(count1/no_sub,3)
        rate[j,2]<-round(count1,0)
        rate[j,3]<-round(count2/no_sub,3)
        rate[j,4]<-round(count2,0)
    }
    
    if(progBar)
    {close(pb)}
    
    colnames(rate)<-c("session 1 percent","session 1 count","session 2 percent","session 2 count")
    
    return(rate)
}
#----
#' Convert Correlation Matrix to Covariance Matrix
#' @description Converts a correlation matrix to a covariance matrix
#' @param cormat A correlation matrix
#' @param data The dataset the correlation matrix is from
#' @return Returns a covariance matrix
#' @examples
#' cormat <- cor(neoOpen)
#' 
#' covmat <- cor2cov(cormat,neoOpen)
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Corrlation to covariance----
cor2cov <- function (cormat, data)
{
    sds<-apply(data,2,sd)
    
    b<-sds%*%t(sds)
    
    S<-cormat*b
    
    return(S)
}
#----
#' Determines if Network is Graphical
#' @description Tests for whether the network is graphical.
#' Input must be a partial correlation network.
#' Function assumes that partial correlations were computed from a multivariate normal distribution
#' @param A A partial correlation network (adjacency matrix)
#' @return Returns a TRUE/FALSE for whether network is graphical
#' @examples
#' A <- LoGo(neoOpen, normal = TRUE, partial = TRUE)$logo
#' 
#' is.graphical(A)
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Is network graphical?----
is.graphical <- function (A)
{
    #make diagonal 1
    diag(A)<-1
    
    #covert partial correlations to covariance
    I<-diag(1, dim(A)[1])
    
    #compute covariance matrix
    error<-try(solve(I-A)%*%t(solve(I-A)),silent=TRUE)
    if(is.character(error))
    {ret<-"FALSE"
    }else{covmat<-solve(I-A)%*%t(solve(I-A))}
    
    #covert to inverse covariance
    error<-try(solve(covmat),silent=TRUE)
    if(is.character(error))
    {ret<-"FALSE"
    }else{inv<-solve(covmat)
    
    #reduce small values to 0
    check<-zapsmall(inv)
    
    
        error<-try(any(eigen(check)$values<0),silent=TRUE)
        if(is.character(error))
        {ret<-"FALSE"
        }else if(any(eigen(check)$values<0))
        {ret<-"FALSE"
        }else{ret<-"TRUE"}
    }
    cat("Is network graphical:",ret)
}
#----
#NEO-PI-3 Openness to Experience Data----
#'NEO-PI-3 Openness to Experience Data
#' 
#' A response matrix (n = 802) of NEO-PI-3's Openness to Experience
#' from Christensen, Cotter, & Silvia (under review).
#' 
#' @docType data
#' 
#' @usage data(neoOpen)
#' 
#' @format A 802x48 response matrix
#' 
#' @keywords datasets
#' 
#' @references
#' Christensen, A.P., Cotter, K.N., & Silvia, P.J. (under review).
#' Reopening Openness to Experience:
#' A network analysis of four Openness to Experience inventories.
#' \url{http://doi.org/10.17605/OSF.IO/954A7}
#' @examples 
#' data(neoOpen)
"neoOpen"
#----
#Animal Verbal Fluency Data----
#'Animal Verbal Fluency Data
#'
#' Verbal fluency reponses for the animal category. The group variable specifies
#' people's level of NEO-FFI-3 Openness to Experience (low = 1, moderate = 2, high = 3).
#' 
#' A binary verbal fluency response matrix (n = 354) from
#' Christensen, Kenett, Cotter, Beaty, & Silvia (under review).
#' 
#' @docType data
#' 
#' @usage data(animals)
#' 
#' @format A 802x48 response matrix
#' 
#' @keywords datasets
#' 
#' @references
#' Christensen, A.P., Kenett, Y. N., Cotter, K.N., Beaty, R. E., & Silvia, P.J. (under review).
#' Remotely close associations:
#' Openness to Experience and semantic memory structure
#' @examples 
#' data(animals)
"animals"
#----
#NEO-PI-3 for Resting-state Data----
#'NEO-PI-3 for Resting-state Data
#'
#' NEO-PI-3 Openness to Experience associated with resting-state data (\emph{n} = 144).
#' 
#' Behavioral data of NEO-PI-3 associated with each connectivity matrix (open)
#' 
#' @docType data
#' 
#' @usage data(behavOpen)
#' 
#' @format behavOpen (vector, length = 144)
#' 
#' @keywords datasets
#' 
#' @references
#' Beaty, R. E., Chen, Q., Christensen, A. P., Qiu, J., Silvia, P. J., & Schacter, D. L. (2018).
#' Brain networks of the imaginative mind: Dynamic functional connectivity of default and cognitive control networks relates to Openness to Experience.
#' \emph{Human Brain Mapping}, \emph{39}(2), 811-821.
#' 
#' Beaty, R. E., Kenett, Y. N., Christensen, A. P., Rosenberg, M. D., Benedek, M., Chen, Q., ... & Silvia, P. J. (2018).
#' Robust prediction of individual creative ability from brain functional connectivity.
#' \emph{Proceedings of the National Academy of Sciences}, 201713532.
#' @examples 
#' data(behavOpen)
"behavOpen"
#----