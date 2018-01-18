#  NetworkToolbox Functions
#
#' Triangulated Maximally Filtered Graph
#' @description Applies the Triangulated Maximally Filtered Graph (TMFG) filtering method
#' @param data Can be a dataset or a correlation matrix
#' @param binary Is dataset dichotomous? Defaults to FALSE. Set TRUE if dataset is dichotomous (tetrachoric correlations are computed)
#' @param weighted Should network be weighted? Defaults to TRUE. Set FALSE to produce an unweighted (binary) network
#' @param depend Is network a dependency (or directed) network? Defaults to FALSE. Set TRUE to generate a TMFG-filtered dependency network
#' @return Returns a list of the adjacency matrix (A), separators (separators), and cliques (cliques)
#' @examples
#' weighted_TMFGnetwork<-TMFG(hex)
#' 
#' weighted_binary_TMFGnetwork<-TMFG(hexb,binary=TRUE)
#' 
#' unweighted_TMFGnetwork<-TMFG(hex,weighted=FALSE)
#' 
#' unweighted_binary_TMFGnetwork<-TMFG(hexb,binary=TRUE,weighted=FALSE)
#' @references 
#' Massara, G. P., Di Matteo, T., & Aste, T. (2016).
#' Network filtering for big data: Triangulated maximally filtered graph.
#' \emph{Journal of Complex Networks}, \emph{5}(2), 161-178.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @importFrom stats cor sd runif qt
#' @export
#TMFG Filtering Method----
TMFG <-function (data, binary = FALSE, weighted = TRUE, depend = FALSE)
{
    if(nrow(data)==ncol(data)){cormat<-data}else
        if(binary){cormat<-psych::tetrachoric(data)$rho}else{cormat<-cor(data)}
    n<-ncol(cormat)
    if(n<9){print("Matrix is too small")}
    nodeTO<-array()
    nodeFROM<-array()
    nodeWEIGHT<-array()
    count<-0
    for(i in 1:nrow(cormat))
        for(j in 1:ncol(cormat))
            if(cormat[i,j] != 0)
            {
                count<-count+1
                nodeTO[count]<-i
                nodeFROM[count]<-j
                nodeWEIGHT[count]<-cormat[i,j]
            }
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
            {if(x[r,z]==1){x[r,z]<-cormat[r,z]}
            }
    }
    
    x<-as.matrix(x)
    colnames(x)<-colnames(cormat)
    rownames(x)<-colnames(cormat)
    return(list(A=x, separators=separators, cliques=cliques))
}
#----
#' Local/Global Sparse Inverse Covariance Matrix
#' @description Applies the Local/Global method to estimate the sparse inverse covariance matrix
#' @param data Must be a dataset
#' @param separators Defaults to separators obtained from the TMFG function. Requires a list of separators
#' @param cliques Defaults to cliques obtained from the TMFG function. Requires a list of cliques
#' @return Returns a sparse TMFG-filtered inverse covariance matrix
#' @examples
#' LoGonet<-LoGo(hex)
#' @references 
#' Barfuss, W., Massara, G. P., Di Matteo, T., & Aste, T. (2016).
#' Parsimonious modeling with information filtering networks.
#' \emph{Physical Review E}, \emph{94}(6), 062306.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @importFrom stats cov
#' @export
#LoGo Sparse Inverse Covariance Matrix----
LoGo <- function (data, separators = TMFG(data)$separators, cliques = TMFG(data)$cliques)
{
    S<-cov(data)
    
    n<-ncol(S)
    Jlogo<-matrix(0,nrow=n,ncol=n)
    
    for(i in 1:nrow(cliques))
    {v<-cliques[i,]
    Jlogo[v,v]<-Jlogo[v,v]+solve(S[v,v])}
    
    for(i in 1:nrow(separators))
    {
        v<-separators[i,]
        Jlogo[v,v]<-Jlogo[v,v]-solve(S[v,v])
    }
    
    colnames(Jlogo)<-colnames(data)
    row.names(Jlogo)<-colnames(data)
    
    return(Jlogo)
}
#----
#' Maximum Spanning Tree
#' @description Applies the Maximum Spanning Tree (MaST) filtering method
#' @param data Can be a dataset or a correlation matrix
#' @param binary Is dataset dichotomous? Defaults to FALSE. Set TRUE if dataset is dichotomous (tetrachoric correlations are computed)
#' @param weighted Should network be weighted? Defaults to TRUE. Set FALSE to produce an unweighted (binary) network
#' @param depend Is network a dependency (or directed) network? Defaults to FALSE. Set TRUE to generate a MaST-filtered dependency network
#' @return A sparse association matrix
#' @examples
#' weighted_MaSTnetwork<-MaST(hex)
#' 
#' weighted_binary_MaSTnetwork<-MaST(hexb,binary=TRUE)
#' 
#' unweighted_MaSTnetwork<-MaST(hex,weighted=FALSE)
#' 
#' unweighted_binary_MaSTnetwork<-MaST(hexb,binary=TRUE,weighted=FALSE)
#' @references 
#' Adapted from: \url{https://www.mathworks.com/matlabcentral/fileexchange/23276-maximum-weight-spanning-tree--undirected}
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Maximum Spanning Tree----
MaST <- function (data, binary = FALSE, weighted = TRUE, depend = FALSE)
{
FIND_PathCompression <- function (temproot=temproot)
{
  ParentPointer[temproot]
  if(ParentPointer[temproot]!=temproot)
  {ParentPointer[temproot]<-FIND_PathCompression(ParentPointer[temproot])}
  parent<-ParentPointer[temproot]
}
if(nrow(data)==ncol(data)){cormat<-data}else
  if(binary){cormat<-psych::tetrachoric(data)$rho}else{cormat<-cor(data)}
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
#' @param weighted Should network be weighted? Defaults to TRUE. Set FALSE to produce an unweighted (binary) network
#' @param binary Is dataset dichotomous? Defaults to FALSE. Set TRUE if dataset is dichotomous (tetrachoric correlations are computed)
#' @param directed Is the network directed? Defaults to FALSE. Set TRUE if the network is directed
#' @return A sparse association matrix
#' @examples
#' weighted_undirected_ECOnetwork<-ECO(hex)
#' 
#' unweighted_undirected_ECOnetwork<-ECO(hex,weighted=FALSE)
#' 
#' weighted_directed_ECOnetwork<-ECO(hex,directed=TRUE)
#' 
#' unweighted_directed_ECOnetwork<-ECO(hex,weighted=FALSE,directed=TRUE)
#' 
#' weighted_undirected_binary_ECOnetwork<-ECO(hexb,binary=TRUE)
#' 
#' unweighted_undirected_binary_ECOnetwork<-ECO(hexb,weighted=FALSE,binary=TRUE)
#' 
#' weighted_directed_binary_ECOnetwork<-ECO(hexb,directed=TRUE,binary=TRUE)
#' 
#' unweighted_directed_binary_ECOnetwork<-ECO(hexb,weighted=FALSE,directed=TRUE,binary=TRUE)
#' @references 
#' Fallani, F. D. V., Latora, V., & Chavez, M. (2017).
#' A topological criterion for filtering information in complex brain networks.
#' \emph{PLoS Computational Biology}, \emph{13}(1), e1005305.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#ECO Neural Network Filter----
ECO <- function (data, weighted = TRUE, binary = FALSE, directed = FALSE)
{
  if(nrow(data)==ncol(data)){C<-data}else
    if(binary){C<-psych::tetrachoric(data)$rho}else{C<-cor(data)}
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
  
  x<-S[ind]
  y<-ind
  h<-cbind(ind,S[ind])
  sorind<-h[order(-h[,2]),]
  C[sorind[(numcon+1):nrow(sorind),1]]<-0
  
  if(directed)
  {W<-C}else{W<-C+t(C)
  diag(W)<-1}
  J<-S+t(S)
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
#' @param binary Is dataset dichotomous? Defaults to FALSE. Set TRUE if dataset is dichotomous (tetrachoric correlations are computed)
#' @param weighted Should network be weighted? Defaults to TRUE. Set FALSE to produce an unweighted (binary) network
#' @return A sparse association matrix
#' @examples
#' weighted_ECOplusMaSTnetwork<-ECOplusMaST(hex)
#' 
#' weighted_binary_ECOplusMaSTnetwork<-ECOplusMaST(hexb,binary=TRUE)
#' 
#' unweighted_binary_ECOplusMaSTnetwork<-ECOplusMaST(hex,weighted=FALSE)
#' 
#' unweighted_binary_ECOplusMaSTnetwork<-ECOplusMaST(hexb,binary=TRUE,weighted=FALSE)
#' @references 
#' Fallani, F. D. V., Latora, V., & Chavez, M. (2017).
#' A topological criterion for filtering information in complex brain networks.
#' \emph{PLoS Computational Biology}, \emph{13}(1), e1005305.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#ECO Filter + MaST----
ECOplusMaST <- function (data, weighted = TRUE, binary = FALSE)
{
  if(weighted&&!binary)
{
  a<-MaST(data,weighted=TRUE)
  b<-ECO(data,weighted=TRUE)
  k<-matrix(NA,nrow=nrow(a),ncol=ncol(a))
  for(i in 1:nrow(a))
    for(j in 1:ncol(a))
      if(a[i,j]==b[i,j]){k[i,j]<-a[i,j]}else k[i,j]<-a[i,j]+b[i,j]
}else if(weighted&&binary)
{
  a<-MaST(data,weighted=TRUE,binary=TRUE)
  b<-ECO(data,weighted=TRUE,binary=TRUE)
  k<-matrix(NA,nrow=nrow(a),ncol=ncol(a))
  for(i in 1:nrow(a))
    for(j in 1:ncol(a))
      if(a[i,j]==b[i,j]){k[i,j]<-a[i,j]}else k[i,j]<-a[i,j]+b[i,j]
}else if(!weighted&&binary)
{
  a<-MaST(data,weighted=FALSE,binary=TRUE)
  b<-ECO(data,weighted=FALSE,binary=TRUE)
  k<-matrix(NA,nrow=nrow(a),ncol=ncol(a))
  for(i in 1:nrow(a))
    for(j in 1:ncol(a))
      if(a[i,j]==b[i,j]){k[i,j]<-a[i,j]}else k[i,j]<-a[i,j]+b[i,j]         
}else if(!weighted&&!binary)
{
  a<-MaST(data,weighted=FALSE,binary=FALSE)
  b<-ECO(data,weighted=FALSE,binary=FALSE)
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
#' @description Filters the network based on an r-value, alpha, boneferroni, or false-discovery rate (FDR)
#' @param data Can be a dataset or a correlation matrix
#' @param binary Is dataset dichotomous? Defaults to FALSE. Set TRUE if dataset is dichotomous (tetrachoric correlations are computed)
#' @param thresh Sets threshold (defaults to \emph{r} = .10). Set to "alpha" to use an alpha value, "bonferroni" for the bonferroni correction, and "FDR" for local false discovery rate
#' @param a Defaults to .05. Applied when thresh = "alpha" and "bonferroni"
#' @return Returns a list containing a filtered adjacency matrix (A) and the critical r value (r.cv)
#' @examples
#' threshnet<-threshold(hex)
#' 
#' alphanet<-threshold(hex, thresh = "alpha")
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @importFrom grDevices dev.off
#' @importFrom utils capture.output
#' @export
#Threshold filtering----
threshold <- function (data, binary = FALSE, thresh = c(.10,"alpha","bonferroni","FDR"), a = .05)
{
    if(missing(thresh))
    {thresh<-.10
    }else{thresh<-match.arg(thresh)}
    
    if(nrow(data)==ncol(data)){cormat<-data}else
        if(binary){cormat<-psych::tetrachoric(data)$rho}else{cormat<-cor(data)}
    
    critical.r <- function(nrow, a){
        df <- nrow - 2
        critical.t <- qt( a/2, df, lower.tail = F )
        cvr <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
        return(cvr)}
    
    if(thresh=="alpha")
    {thr<-critical.r(nrow(data),a)
    }else if(thresh=="bonferroni")
    {thr<-critical.r(nrow(data),(a/((ncol(cormat)^2)-(ncol(cormat))/2)))
    }else if(thresh=="FDR")
    {
        fdrmat<-matrix(0,nrow=((ncol(cormat)^2)-(ncol(cormat))),ncol=3)
        w<-0
        for(i in 1:ncol(cormat))
            for(j in 1:ncol(cormat))
                if(i!=j)
        {
            w<-w+1
            fdrmat[w,1]<-i
            fdrmat[w,2]<-j
            fdrmat[w,3]<-cormat[i,j]
        }
        
        fdrmat[,3]<-fdrtool::fdrtool(fdrmat[,3],plot=FALSE,cutoff.method=
                                  "locfdr",statistic = "correlation")$qval
    }else thr<-thresh
    
    if(!thresh=="FDR")
    {cormat<-ifelse(cormat>=thr,cormat,0)
    }else if(thresh=="FDR")
    {
        fdrmat[,3]<-ifelse(fdrmat[,3]<=a,fdrmat[,3],0)
        fdrmat<-as.matrix(Matrix::sparseMatrix(i=fdrmat[,1],j=fdrmat[,2],x=fdrmat[,3]))
        cormat<-ifelse(fdrmat!=0,cormat,0)
        thr<-min(cormat[cormat!=0])
    }
    
    diag(cormat)<-1
    
    return(list(A=cormat, r.cv=thr))
}
#----
#' Betwenness Centrality
#' @description Computes betweenness centrlaity of each node in a network
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted? Defaults to TRUE. Set to FALSE for unweighted measure of betwenness centrality
#' @return A vector of betweenness centrality values for each node in the network
#' @examples
#' A<-TMFG(hex)$A
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
#' @param A An adjacency matrix of network data
#' @param beta Sets the beta parameter. Defaults to 0.01 (recommended). Beta > 0.01 measure gets closer to weighted betweenness centrality (10) and beta < 0.01 measure gets closer to degree (.0001)
#' @return A vector of randomized shortest paths betweenness centrality values for each node in the network
#' @examples
#' A<-TMFG(hex)$A
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
#' @param weighted Is the network weighted? Defaults to TRUE. Set to FALSE for unweighted measure of closeness centrality
#' @return A vector of closeness centrality values for each node in the network
#' @examples
#' A<-TMFG(hex)$A
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
#' @return A vector of degree values for each node in the network. If directed network, returns a list of in-degree (inDegree), out-degree (outDegree), and relative influence (relInf)
#' @examples
#' A<-TMFG(hex)$A
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
  return(list(inDegree=inDeg,outDegree=outDeg,relInf=relinf))
  }
}
#----
#' Node Strength
#' @description Computes strength of each node in a network
#' @param A An adjacency matrix of network data
#' @return A vector of strength values for each node in the network. If directed network, returns a list of in-strength (inStrength), out-strength (outStrength), and relative influence (relInf)
#' @examples
#' A<-TMFG(hex)$A
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
      return(list(inStrength=inStr,outStrength=outStr,relInf=relinf))
  }
}
#----
#' Eigenvector Centrality
#' @description Computes eigenvector centrality of each node in a network
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted? Defaults to TRUE. Set to FALSE for unweighted measure of eigenvector centrality
#' @return A vector of eigenvector centrality values for each node in the network
#' @examples
#' A<-TMFG(hex)$A
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
#' @description Computes leverage centrlaity of each node in a network
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted? Defaults to TRUE. Set to FALSE for unweighted measure of leverage centrality
#' @return A vector of leverage centrality values for each node in the network
#' @examples
#' A<-TMFG(hex)$A
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
#' @description Computes impact measure of each node in a network
#' @param A An adjacency matrix of network data
#' @return A vector of node impact values for each node in the network
#' @examples
#' A<-TMFG(hex)$A
#' 
#' nodeimpact<-impact(A)
#' @references 
#' Kenett, Y. N., Kenett, D. Y., Ben-Jacob, E., & Faust, M. (2011).
#' Global and local features of semantic networks: Evidence from the Hebrew mental lexicon.
#' \emph{PloS one}, \emph{6}(8), e23912.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Node Impact----
impact <- function (A)
{
    allP<-pathlengths(A)$ASPL
    remove<-matrix(0,nrow=nrow(A),ncol=1)
    for(j in 1:ncol(A))
    {remove[j,]<-(pathlengths(A[-j,-j])$ASPL)-allP}
    remove<-round(remove,3)
    colnames(remove)<-"Impact"
    row.names(remove)<-colnames(A)
    return(remove)
}
#----
#' Hybrid Centrality
#' @description Computes hybrid centrality of each node in a network
#' @param A An adjacency matrix of network data
#' @return A vector of hybrid centrality values for each node in the network (higher values are more central, lower values are more peripheral)
#' @examples
#' A<-TMFG(hex)$A
#' 
#' HC<-hybrid(A)
#' @references 
#' Pozzi, F., Di Matteo, T., & Aste, T. (2013).
#' Spread of risk across financial markets: Better to invest in the peripheries. 
#' \emph{Scientific Reports}, \emph{3}(1655), 1-7.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Hybrid Centality----
hybrid <- function (A)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  BCu<-betweenness(A,weighted=FALSE)
  BCw<-betweenness(A)
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
            8)/(8*((ncol(A))-8)))
  hybrid<-round(as.data.frame(hybrid),3)
  colnames(hybrid)<-c("HC")
  rownames(hybrid)<-colnames(A)
  hybrid<-as.matrix(hybrid)
  return(hybrid)
}
#----
#' List of Centrality Measures
#' @description Computes centrality measures of the network
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted? Defaults to TRUE. Set to FALSE for unweighted list of centrality measures
#' @return Returns a list of betweenness, closeness, degree (weighted = strength), eigenvector, and leverage centralities
#' @examples
#' A<-TMFG(hex)$A
#' 
#' weighted_centralitylist<-centlist(A)
#' 
#' unweighted_centralitylist<-centlist(A,weighted=FALSE)
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
  list(Betweenness=BC,Closeness=CC,Degree=Deg,Eigenvector=EC)}else{
    BC<-betweenness(A)
    CC<-closeness(A)
    Str<-strength(A)
    EC<-eigenvector(A)
    lev<-leverage(A)
    return(list(betweenness=BC,closeness=CC,strength=Str,eigenvector=EC,leverage=lev))}
}
#----
#' Distance
#' @description Computes distance matrix of the network
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted? Defaults to FALSE. Set to TRUE for weighted measure of distance
#' @return A distance matrix of the network
#' @examples
#' A<-TMFG(hex)$A
#' 
#' unweighted_D<-distance(A)
#' 
#' weighted_D<-distance(A,weighted=TRUE)
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
#' @description Computes global average shortest path length (ASPL), local average shortest path length (ASPLi), eccentricity (ecc), and diameter (D) of a network
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted? Defaults to FALSE. Set to TRUE for weighted measures of ASPL, ASPLi, ecc, and D
#' @return Returns a list of ASPL, ASPLi, ecc, and D of a network
#' @examples
#' A<-TMFG(hex)$A
#' 
#' unweighted_PL<-pathlengths(A)
#' 
#' weighted_PL<-pathlengths(A,weighted=TRUE)
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
#' @param weighted Is the network weighted? Defaults to FALSE. Set to TRUE for weighted measures of CC and CCi
#' @return Returns a list of CC and CCi
#' @examples
#' A<-TMFG(hex)$A
#' 
#' unweighted_CC<-clustcoeff(A)
#' 
#' weighted_CC<-clustcoeff(A,weighted=TRUE)
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
#' @param weighted Is the network weighted? Defaults to FALSE. Set to TRUE for a weighted measure of transitivity
#' @return Returns a value of transitivity
#' @examples
#' A<-TMFG(hex)$A
#' 
#' trans<-transitivity(A,weighted=TRUE)
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
#' @param gamma Defaults to 1. Set to gamma > 1 to detect smaller modules and gamma < 1 for larger modules
#' @param M0 Defaults to none. Input can be an initial community vector
#' @param method Defaults to "modularity". Set to "potts" for Potts model
#' @return Returns a list of community and Q
#' @examples
#' A<-TMFG(hex)$A
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
louvain <- function (A, gamma = 1, M0 = 1:ncol(A), method = c("modularity","potts"))
{
    if(missing(method))
    {method<-"modularity"
    }else{method<-match.arg(method)}
    
    n<-ncol(A)
    s<-sum(A)
    
    if(min(A)<0)
    {warning("Matrix contains negative weights: absolute weights were used")
        A<-abs(A)}
    
    Mb<-unique(M0)
    M<-Mb
    
    if(method=="modularity") 
    {mat<-matrix(0,nrow=n,ncol=n)
    for(i in 1:n)
        for(j in 1:n)
        {mat[i,j]<-(colSums(A)[i]*rowSums(A)[j])/s}
    
    B<-A-(gamma*(mat))
    }else if(method=="potts"){
        B<-A-gamma*!A
    }
    
    
    B<-(B+t(B))/2
    
    Hnm<-matrix(0,nrow=n,ncol=n)
    
    for(m in 1:max(Mb))
    {Hnm[,m]<-B[,Mb==m]}
    
    H<-colSums(Hnm)
    Hm<-rowSums(Hnm)
    
    Q0<-(-Inf)
    bsxfun<-matrix(0,nrow=n,ncol=n)
    diag(bsxfun)<-1
    Q<-sum(diag(B*bsxfun))/s
    
    
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
#' @param iter Number of random (or lattice) networks to generate, which are used to calculate the mean random ASPL and CC (or lattice)
#' @param progBar Defaults to FALSE. Set to TRUE to see progress bar
#' @param method Defaults to "HG" (Humphries & Gurney, 2008). Set to "rand" for the CC to be calculated using a random network or
#' set to "TJHBL" for (Telesford et al., 2011) where CC is calculated from a lattice network
#' @return Returns a value of small-worldness.
#' For "rand", values > 1 indicate a small-world network.
#' For "HG", values > 3 indicate a small-world network.
#' For "TJHBL" values near 0 indicate a small-world network
#' while < 0 indicates a more regular network and > 0 
#' indicates a more random network
#' @examples
#' A<-TMFG(hex)$A
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
    }else if(method=="HG")
    {rtrans<-mean(csamps)
    swm<-(trans/rtrans)/(ASPL/rASPL)
    }else if(method=="TJHBL")
    {lCC<-mean(csamps)
    swm<-(rASPL/ASPL)-(CC/lCC)}
    
    return(swm)
}
#----
#' Semantic Network Measures
#' @description Computes the average shortest path length (ASPL), clustering coefficient(CC),
#' modularity (Q), and small-worldness (S) 
#' @param A An adjacency matrix of network A
#' @return Returns a values for ASPL, CC, Q, and S
#' @examples
#' A<-TMFG(hex)$A
#' 
#' connectivity<-semnetmeas(A)
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Semantic Network Measures----
semnetmeas <- function (A)
{
    aspl<-pathlengths(A)$ASPL
    cc<-clustcoeff(A)$CC
    q<-louvain(A)$Q
    s<-smallworldness(A,iter=100,progBar = FALSE,method="rand")
    
    semnetmeas<-cbind(aspl,cc,q,s)
    
    semnetmeas<-as.data.frame(semnetmeas)
    
    colnames(semnetmeas)<-c("ASPL","CC","Q","S")
    
    semnetmeas<-as.matrix(semnetmeas)
    
    return(semnetmeas)
}
#----
#' Edge Replication
#' @description Computes the number of edges that replicate between two cross-sectional networks
#' @param A An adjacency matrix of network A
#' @param B An adjacency matrix of network B
#' @return Returns a list of the number of edges that replicate (replicated), total number of edges (possibleA & possibleB), the percentage of edges that replicate (percentageA & percentageB), the density of edges (densityA & densityB), the mean difference between edges that replicate (meanDifference), the sd of the difference between edges that replicate (sdDifference), and the correlation between the edges that replicate for both networks (correlation)
#' @examples
#' A<-TMFG(hex)$A
#' 
#' B<-MaST(hex)
#' 
#' edges<-edgerep(A,B)
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Edge Replication----
edgerep <- function (A, B)
{
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
    warning("Adjacency matrix A was made to be symmetric")}
}

count<-0
  for(i in 1:ncol(A))
    for(j in 1:nrow(A))
        if(i!=j)
            {if(A[i,j]&&B[i,j]!=0){count<-count+1}}
                count<-count/2
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
  corr<-cor(wveca,wvecb)
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
  
  return(list(replicated=count,
              totalEdgesA=possibleA,totalEdgesB=possibleB,
              percentageA=percentA,percentageB=percentB,
              densityA=densityA,densityB=densityB,
              meanDifference=mvec,sdDifference=svec,correlation=corr))
}
#----
#' Network Connectivity
#' @description Computes the average and standard deviation of the weights in the network
#' @param A An adjacency matrix of network A
#' @return Returns a list of the edge weights (weights), the mean (mean), the standard deviation (sd), and the sum of the edge weights (total) in the network
#' @examples
#' A<-TMFG(hex)$A
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
#' @description Bootstraps the sample to identify the most stable correlations
#' @param data A set of data
#' @param method A network filtering method. Defaults to "TMFG"
#' @param binary Is dataset dichotomous? Defaults to FALSE. Set TRUE if dataset is dichotomous (tetrachoric correlations are computed)
#' @param n Number of people to use in the bootstrap. Defaults to full sample size
#' @param iter Number of bootstrap iterations. Defaults to 1000 iterations
#' @param depend Is network a dependency (or directed) network? Defaults to FALSE. Set TRUE to generate a TMFG-filtered dependency network
#' @param ... Additional arguments for filtering methods
#' @return Returns a list that includes the original filtered network (orignet),
#' correlation matrix of the mean bootstrapped network (bootmat),
#' reliabilities of the connections in the original network (netrel),
#' reliabilities of the connections in the bootstrapped network (bootrel),
#' a plot of the bootrel reliability matrix (netrel; upper triangle = actual network reliabilites, bootrel; lower triangle = overall network reliablities),
#' a plot of included correlations on their reliability (ConR)
#' @examples
#' \dontrun{
#' prepTMFG<-bootgen(hex)
#' 
#' prepLoGo<-bootgen(hex,method="LoGo")
#' 
#' prepMaST<-bootgen(hex,method="MaST")
#' 
#' prepECO<-bootgen(hex,method="ECO")
#' 
#' prepECOplusMaST<-bootgen(hex,method="ECOplusMaST")
#' 
#' prepThreshold<-bootgen(hex,method="threshold")
#' }
#' @references
#' Tumminello, M., Coronnello, C., Lillo, F., Micciche, S., & Mantegna, R. N. (2007).
#' Spanning trees and bootstrap reliability estimation in correlation-based networks.
#' \emph{International Journal of Bifurcation and Chaos}, \emph{17}(7), 2319-2329.
#' 
#' Wei, T. & Simko, V.(2017).
#' R package "corrplot": Visualization of a correlation matrix (Version 0.84).
#' Available from \url{https://github.com/taiyun/corrplot}
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @importFrom graphics abline plot text
#' @importFrom stats lm na.omit cov2cor
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#Bootstrap Network Generalization----
bootgen <- function (data, method = c("TMFG","LoGo","MaST","ECOplusMaST","ECO","threshold"), binary = FALSE, n = nrow(data), iter = 1000, depend = FALSE, ...)
{
    if(missing(method))
    {method<-"TMFG"
    }else{method<-match.arg(method)}
    
    if(nrow(data)==ncol(data)){stop("Input must be a dataset")}else
        if(binary){realmat<-psych::tetrachoric(data)$rho}else{realmat<-cor(data)}
    mat<-matrix(0,nrow=n,ncol=ncol(data)) #Initialize bootstrap matrix
    samps<-array(0,c(nrow=nrow(realmat),ncol=ncol(realmat),iter)) #Initialize sample matrix
    if(depend)
    {dsamps<-array(0,c(nrow=nrow(realmat),ncol=ncol(realmat),iter))} #Initialize dsample matrix
    
    #fisher's z
    fish <- function (r)
    {z<-.5*log((1+abs(r))/(1-abs(r)))
    if(nrow(r)>1&&ncol(r)>1&&length(r)>1)
    {diag(z)<-0}
    return(z)}
    
    pb <- txtProgressBar(max=iter, style = 3)
    for(i in 1:iter) #Generate array of bootstrapped samples
    {
        f<-round(runif(i,min=1,max=1000000),0)
        set.seed(f[round(runif(i,min=1,max=length(f)),0)])
        mat<-data[round(runif(n,min=1,max=n),0),]
        cormat<-cor(mat)
        if(!depend){
        if(method=="TMFG")
        {samps[,,i]<-fish(TMFG(cormat)$A)
        }else if(method=="LoGo")
        {l<-fish((-cov2cor(LoGo(mat))))
        samps[,,i]<-l
        }else if(method=="MaST")
        {samps[,,i]<-fish(MaST(cormat,...))
        }else if(method=="ECOplusMaST")
        {samps[,,i]<-fish(ECOplusMaST(cormat,...))
        }else if(method=="ECO")
        {samps[,,i]<-fish(ECO(cormat,...))
        }else if(method=="threshold")
        {samps[,,i]<-fish(threshold(cormat,...)$A)
        }else {stop("Method not available")}
        }else{if(method=="TMFG")
        {samps[,,i]<-fish(TMFG(cormat)$A)
        dsamps[,,i]<-TMFG(depend(cormat,progBar=FALSE),depend=TRUE)$A
        }else stop("Method not available")
        }
        setTxtProgressBar(pb, i)
    }
    close(pb)
    
    if(!depend)
    {
        if(method=="TMFG")
        {tru<-TMFG(data,...)$A
        }else if(method=="LoGo")
        {tru<-(-cov2cor(LoGo(data)))
        diag(tru)<-0
        }else if(method=="MaST")
        {tru<-MaST(data,...)
        }else if(method=="ECOplusMaST")
        {tru<-ECOplusMaST(data,...)
        }else if(method=="MaST")
        {tru<-MaST(data,...)
        }else if(method=="ECO")
        {tru<-ECO(data,...)
        }else if(method=="threshold")
        {tru<-threshold(data,...)$A
        }else stop("Method not available")
    }
    
    
    zw <- function (z,iter)
    {
        sum((iter-3)*(z))/((iter-3)*iter)
    }
    
    #Mean matrix
    meanmat<-matrix(0,nrow=nrow(realmat),ncol=ncol(realmat)) #Initialize Mean matrix
    for(j in 1:nrow(realmat))
        for(k in 1:ncol(realmat))
        {meanmat[j,k]<-zw(samps[j,k,],iter)}
    
    if(!method=="LoGo")
    {meanmat<-psych::fisherz2r(meanmat)}
    
    #Set alpha
    if(n<=iter)
    {a<-1/iter
    }else if(n>iter)
    {a<-.05}
    
        critical.r <- function(iter, a){
            df <- iter - 2
            critical.t <- qt( a/2, df, lower.tail = F )
            cvr <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
            return(cvr)}
    if(!method=="LoGo")
    {for(x in 1:nrow(meanmat))
        for(y in 1:ncol(meanmat))
            if(meanmat[x,y]<=critical.r(iter,(a/((ncol(data)^2)-(ncol(data)/2)))))
            {meanmat[x,y]<-0}
    }
    
    #return meanmat to bootmat
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
            
            reprel<-rel
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
            wc<-0
            for(i in 1:nrow(cormat))
                for(j in 1:ncol(cormat))
                    if((upp[i,j]!=0&&tru[i,j])!=0)
                    {wc<-wc+1
                    x[wc]<-upp[i,j]
                    y[wc]<-tru[i,j]}
            xo<-na.omit(abs(x))
            yo<-na.omit(abs(y))
            
            mar=c(2,2,2,2)
            cpo<-{plot(xo,yo,pch=16,ylab="Correlation Strength",xlab="Reliability",
                       main="Bootstrapped Correlation Strength on Reliability",xlim=c(0,1),ylim=range(yo))
                abline(lm(yo~xo))
                text(x=.05,y=max(yo-.05),labels = paste("r = ",round(cor(yo,xo),3)))}
            
            
            #plot reliability matrix
            if(ncol(realmat)<=20)
            {plt<-corrplot::corrplot(rel,method="color",
                                     title="Bootstrapped Correlation Reliabilities",
                                     mar=c(2,2,2,2),tl.col="black",tl.cex=.75,
                                     cl.lim = c(0,1),addgrid.col = "grey",addCoef.col = "black")
            }else if(ncol(realmat)>20){
                plt<-corrplot::corrplot(rel,method="color",
                                        title="Bootstrapped Correlation Reliabilities",
                                        mar=c(2,2,2,2),tl.col="black",tl.cex=.75,
                                        cl.lim = c(0,1),addgrid.col = "grey")}
            
            if(depend)
            {   dtru<-TMFG(depend(data,progBar=FALSE),depend=TRUE)$A
                
                #mean matrix
                dmeanmat<-matrix(0,nrow=nrow(realmat),ncol=ncol(realmat)) #Initialize Mean matrix
                for(j in 1:nrow(realmat))
                    for(k in 1:ncol(realmat))
                    {dmeanmat[j,k]<-mean(dsamps[j,k,])}
                    
                    for(x in 1:nrow(dmeanmat))
                        for(y in 1:ncol(dmeanmat))
                            if(dmeanmat[x,y]<=critical.r(iter))
                            {dmeanmat[x,y]<-0}
                #bootmat
                dbootmat<-dmeanmat
                colnames(dbootmat)<-colnames(dbootmat)
                #reliability count
                dsamp<-array(0,c(nrow=nrow(realmat),ncol=ncol(realmat),iter))
                drel<-matrix(0,nrow=nrow(realmat),ncol=ncol(realmat))
                
                for(j in 1:nrow(realmat))
                for(k in 1:ncol(realmat))
                    for(l in 1:iter)
                        if(dsamps[j,k,l]!=0)
                        {dsamp[j,k,l]<-1}
                #reliability
                for(j in 1:nrow(realmat))
                    for(k in 1:ncol(realmat))
                        drel[j,k]<-sum(dsamp[j,k,])/iter
                    colnames(drel)<-colnames(data)
                #reliablity plot
                dreprel<-drel
                row.names(drel)<-colnames(drel)
                diag(drel)<-1
                #reliablity on correlation plot
                dx<-matrix(nrow=length(drel))
                dy<-matrix(nrow=length(dtru))
                dwc<-0
                for(i in 1:nrow(cormat))
                    for(j in 1:ncol(cormat))
                        if((drel[i,j]!=0&&dtru[i,j])!=0)
                        {dwc<-dwc+1
                        dx[dwc]<-drel[i,j]
                        dy[dwc]<-dtru[i,j]}
            dxo<-na.omit(dx)
            dyo<-na.omit(dy)
            
            mar=c(2,2,2,2)
            dcpo<-{plot(dxo,dyo,pch=16,ylab="Dependency Strength",xlab="Reliability",
                       main="Dependency Strength on Reliability",xlim=c(0,1),ylim=range(dyo))
                abline(lm(dyo~dxo))
                text(x=.05,y=max(dyo-.05),labels = paste("r = ",round(cor(dyo,dxo),3)))}
            
            #plot reliability matrix
            if(ncol(realmat)<=20)
            {dplt<-corrplot::corrplot(drel,method="color",
                                     title="Bootstrapped Dependency Reliabilities",
                                     mar=c(2,2,2,2),tl.col="black",tl.cex=.75,
                                     cl.lim = c(0,1),addgrid.col = "grey",addCoef.col = "black")
            }else if(ncol(realmat)>20){
                dplt<-corrplot::corrplot(drel,method="color",
                                        title="Bootstrapped Dependency Reliabilities",
                                        mar=c(2,2,2,2),tl.col="black",tl.cex=.75,
                                        cl.lim = c(0,1),addgrid.col = "grey")}
            }
    
    diag(bootmat)<-1        
    
    if(!depend)
    {return(list(orignet=tru,bootmat=bootmat,netrel=upp,bootrel=reprel,plotrel=plt,ConR=cpo))
    }else{orignet<-(list(orignet=tru,bootmat=bootmat,bootrel=reprel,plotrel=plt,ConR=cpo))
          depnet<-(list(orignet=dtru,bootmat=dbootmat,bootrel=dreprel,plotrel=dplt,DonR=dcpo))
          return(list(undirected=orignet,directed=depnet))}
}
#----
#' Bootstrapped Communities Likelihood
#' @description Bootstraps the sample with replace to compute walktrap reliability
#' @param data A set of data
#' @param binary Is dataset dichotomous? Defaults to FALSE. Set TRUE if dataset is dichotomous (tetrachoric correlations are computed)
#' @param n Number of people to use in the bootstrap. Defaults to full sample size
#' @param iter Number of bootstrap iterations. Defaults to 100 iterations
#' @param filter Set filter method. Defaults to "TMFG"
#' @param method Defaults to "louvain". Set to "walktrap" for the walktrap algorithm
#' @param steps Number of steps to use in the walktrap algorithm. Defaults to 4. Use a larger number of steps for smaller networks
#' @param ... Additional arguments for network filtering methods
#' @return The factors and their proportion found across bootstrapped samples (i.e., their likelihood)
#' @examples
#' commTMFG<-commboot(hex)
#' 
#' commThreshold<-commboot(hex,filter="threshold")
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
commboot <- function (data, binary = FALSE, n = nrow(data), iter = 100, filter = c("TMFG","threshold"), method = c("louvain","walktrap"), steps = 4, ...)
{
    if(missing(filter))
    {filter<-"TMFG"
    }else{filter<-match.arg(filter)}
    
    if(missing(method))
    {method<-"louvain"
    }else{method<-match.arg(method)}
    
    col<-ncol(data)
    if(nrow(data)==ncol(data)){stop("Input must be a dataset")}else
        if(binary){realmat<-psych::tetrachoric(data)$rho}else{realmat<-cor(data)}
    mat<-matrix(0,nrow=n,ncol=col) #Initialize bootstrap matrix
    comm<-matrix(0,nrow=iter,ncol=1) #Initialize community matrix
    pb <- txtProgressBar(max=iter, style = 3)
    for(i in 1:iter) #Generate array of bootstrapped samples
    {
        f<-round(runif(i,min=1,max=1000000),0)
        set.seed(f[round(runif(i,min=1,max=length(f)),0)])
        mat<-data[round(runif(n,min=1,max=n),0),]
        if(any(colSums(mat)<=1)){stop("Increase sample size: not enough observations")}
        cormat<-cor(mat)
        if(method=="walktrap")
        {if(filter=="TMFG")
        {comm[i,]<-max(igraph::walktrap.community(igraph::as.igraph(qgraph::qgraph(TMFG(cormat,...)$A,DoNotPlot=TRUE)),steps=steps)$membership)
        }else if(filter=="threshold")
        {comm[i,]<-max(igraph::walktrap.community(igraph::as.igraph(qgraph::qgraph(threshold(cormat,...)$A,DoNotPlot=TRUE)),steps=steps)$membership)
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
    
    return(prop)
}
#----
#' Dependency Matrix
#' @description Generates a dependency matrix of the data
#' @param data A set of data
#' @param binary Is dataset dichotomous? Defaults to FALSE. Set TRUE if dataset is dichotomous (tetrachoric correlations are computed)
#' @param index Should correlation with the latent variable (i.e., weighted average of all variables) be removed? Defaults to FALSE. Set to TRUE to remove common latent factor
#' @param fisher Should Fisher's Z-test be used to keep significantly higher influences (index only)? Defaults to FALSE. Set to TRUE to remove non-significant influences
#' @param progBar Should progress bar be displayed? Defaults to TRUE. Set FALSE for no progress bar.
#' @return Returns an adjacency matrix of dependencies
#' @examples
#' D<-depend(hex)
#' 
#' binaryD<-depend(hexb,binary=TRUE)
#' 
#' Dindex<-depend(hex,index=TRUE)
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
depend <- function (data, binary = FALSE, index = FALSE, fisher = FALSE, progBar = TRUE)
{
    if(nrow(data)==ncol(data)){cormat<-data}else
        if(binary){cormat<-psych::tetrachoric(data)$rho}else{cormat<-cor(data)}
    
    inter<-((ncol(cormat)*(ncol(cormat)-1)*(ncol(cormat)-2)))
    
    if(index)
    {
        m<-rowMeans(data)
        dat<-cbind(data,m)
        if(nrow(dat)==ncol(dat)){cordat<-dat}else
            if(binary){cordat<-psych::tetrachoric(dat)$rho}else{cordat<-cor(dat)}
        
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
        {zsd<-1/sqrt(n-3)
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
#' @return A list containing the training (trainSample) and testing (testSample) sizes (training sample is made to have slightly larger sizes in the case of uneven splits), and their respective sample sizes (trainSize and testSize)
#' @examples
#' splithalf<-splitsamp(hex)
#' @references 
#' Forbes, M. K., Wright, A. G. C., Markon, K. E., & Krueger, R. F. (2017a).
#' Evidence that psychopathology symptom networks do not replicate.
#' \emph{Journal of Abnormal Psychology}, \emph{126}(7), 969-988.
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
#' net<-TMFG(hex)$A
#' 
#' hexb<-binarize(net)
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Binarize function----
binarize <- function (A)
{
    bin<-ifelse(A!=0,1,0)
    
    return(bin)
}
#----
#' Kullback-Leibler Divergence
#' @description Estimates the Kullback-Leibler Divergence which measures how one probability distribution diverges from a second distribution (equivalent means are assumed). Matrices \strong{must} be positive definite for accurate measurement
#' @param base Full or base model (e.g., a correlation or covariance matrix of the data)
#' @param test Reduced or testing model (e.g., a sparse correlation or covariance matrix)
#' @param basedata Full or base dataset to be compared
#' @param testdata Testing dataset to be compared
#' @param corr Are models correlation matrices? Defaults to TRUE. Set to FALSE for covariance matrices
#' @return A value greater than 0. Values between 0 and 1 suggests the reduced or testing model is near the full mode
#' @examples
#' A1 <- cov(hex)
#' 
#' A2 <- solve(LoGo(hex))
#' 
#' kld_value <- kld(A1, A2, hex, corr = FALSE)
#' 
#' @references 
#' Kullback, S., & Leibler, R. A. (1951).
#' On information and sufficiency.
#' \emph{The Annals of Mathematical Statistics}, \emph{22}(1), 79-86.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Kullback-Leibler Divergence----
kld <- function (base, test, basedata, testdata, corr = TRUE)
{
    if(missing(testdata))
    {testdata<-basedata}
    if(corr) 
    {
        cor2cov <- function (A,data)
        {
            sds<-apply(data,2,sd)
            b<-sds%*%t(sds)
            S<-A*b
            return(S)
        }
        
        if(any(eigen(base)$values<0))
        {
            base<-matrix(Matrix::nearPD(base,corr=TRUE)$mat@x,nrow=nrow(base),ncol=ncol(base))
            warning("Base was forced to be positive definite")
        }
        
        if(any(eigen(test)$values<0))
        {
            test<-matrix(Matrix::nearPD(test,corr=TRUE)$mat@x,nrow=nrow(test),ncol=ncol(test))
            warning("Test was forced to be positive definite")
        }
        
        base<-cor2cov(base,basedata)
        test<-cor2cov(test,testdata)
    }else{
        if(any(eigen(base)$values<0))
        {
            base<-matrix(Matrix::nearPD(base)$mat@x,nrow=nrow(base),ncol=ncol(base))
            warning("Base was forced to be positive definite")
        }
        
        if(any(eigen(test)$values<0))
        {
            test<-matrix(Matrix::nearPD(test)$mat@x,nrow=nrow(test),ncol=ncol(test))
        }
    }
    
    
    kld<-.5*(log(det(test)/det(base)) + sum(diag(MASS::ginv(test)%*%base)) - nrow(base))
    
    return(kld)
}
#----
#' Root Mean Square Error
#' @description Computes the root mean square error of a sparse model to a full model
#' @param base Base (or full) model to be evaulated against
#' @param model Reduced (or testing) model (e.g., a sparse correlation or covariance matrix)
#' @param test Data that is not the base dataset (i.e., comparing a training model to test data)
#' @return RMSE value. Lower values suggest more similarity between the full and sparse model
#' @examples
#' A <- solve(LoGo(hex))
#' 
#' root <- rmse(hex, A)
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Root Mean Square Error----
rmse <- function (base, model, test)
{
    if(missing(test))
    {test<-base}
    
    if(nrow(base)==ncol(base))
    {stop("Base must be a dataset")}
    
    if(nrow(test)==ncol(test))
    {stop("Test must be a dataset")}
    
    
    cor2cov <- function (A, data)
    {
        sds<-apply(data,2,sd)
        
        b<-sds%*%t(sds)
        
        S<-A*b
        
        return(S)
    }
    
    if(all(diag(model)==1))
    {mod<-cor2cov(model,test)
    }else if(all(diag(model)==0))
    {
        diag(model)<-1
        mod<-cor2cov(model,test)
    }else{mod<-model}
    
    root <- sqrt(sum((cov(base)-mod)^2)/(ncol(base)^2))
    
    return(root)
}
#----
#' Import CONN Toolbox Brain Matrices to R format and convert to correlations
#' @description Converts a Matlab brain z-score connectivity file (n x n x m) where \strong{n} is the n x n connectivity matrices and \strong{m} is the participant
#' @param MatlabData Input for Matlab data file. Defaults to interactive file choice
#' @param progBar Should progress bar be displayed? Defaults to TRUE. Set FALSE for no progress bar
#' @return Returns an array of correlation connectivity (n x n x m)
#' @examples
#' \dontrun{neuralarray<-convertConnBrainMat()}
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Convert CONN Toolbox Brain Matrices----
convertConnBrainMat <- function (MatlabData = file.choose(), progBar = TRUE)
{
    
    mat<-R.matlab::readMat(file.choose()) #read in matlab data
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
    
    return(dat)
}
#----
#' Dependency Neural Networks
#' @description Applies the dependency network approach to neural network array
#' @param neuralarray Array from \emph{convertConnBrainMat} function
#' @param pB Should progress bar be displayed? Defaults to TRUE. Set FALSE for no progress bar
#' @param ... Additional arguments from \emph{depend} function
#' @return Returns an array of n x n x m dependency matrices
#' @examples
#' \dontrun{neuralarray <- convertConnBrainMat()
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
#' @description Applies a network filtering methodology to neural network array. Removes edges from the neural network output from \emph{convertConnBrainMat} using a network filtering approach
#' @param neuralarray Array from \emph{convertConnBrainMat} function
#' @param progBar Should progress bar be displayed? Defaults to TRUE. Set FALSE for no progress bar
#' @param method Filtering method to be applied (e.g., "MaST")
#' @param ... Additional arguments from filtering methods
#' @return Returns an array of n x n x m filtered matrices
#' @examples
#' \dontrun{neuralarray <- convertConnBrainMat()
#' 
#' filteredneuralarray <- neuralnetfilter(neuralarray, method = "threshold", thres = .50)
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
neuralnetfilter <- function (neuralarray, progBar = TRUE, method = c("TMFG","MaST","ECOplusMaST","ECO","threshold"), ...)
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
        {gstat[i,]<-smallworldness(filarray[,,i],iter=10,...)
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
        {lstat[,i]<-degree(filarray[,,i])$inDegree
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
        {lstat[,i]<-strength(filarray[,,i])$inStrength
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
        {lstat[,i]<-louvain(filarray[,,i],...)$community
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
#' Group-wise Neural Network Statistics Tests
#' @description Statistical test for group differences for global or local network characteristics of neural network data (\strong{only t-tests})
#' @param groups Participant list divided into the desired groups (\strong{see examples})
#' @param nstat A statistic vector (whole-network) or matrix (ROI) from the \emph{neuralstat} function
#' @param correction Multiple comparisons correction for ROI testing. Defaults to local false discovery rate (i.e., "FDR")
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
#' @export
#Neural Groups----
neuralstattest <- function (groups, nstat, correction = c("bonferroni","FDR"))
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
            {cpval<-fdrtool::fdrtool(oval,plot=FALSE,statistic = "pvalue")$qval}
            g1[count2]<-t.test(group1[,i],group2[,i],var.equal = TRUE)$estimate[1]
            g2[count2]<-t.test(group1[,i],group2[,i],var.equal = TRUE)$estimate[2]}
        }
        sig<-cbind(ROI,round(tval,3),round(unpval,5),round(cpval,5),round(g1,3),round(g2,3))
        if(any(sig==0))
        {sig[which(sig==0)]<-"ns"}else{sig[which(sig==1)]<-"sig"}
        colnames(sig)<-c("ROI","t-value","uncorrected p","corrected p","group 1 mean", "group 2 mean")
        if(all(is.na(sig)))
        {print("No significant differences")
        }else{return(sig)}
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
#HEXACO Openness data----
#' HEXACO Openness to Experience Response Matrix
#' 
#' A response matrix (n = 802) of HEXACO's Openness to Experience
#' from Christensen, Cotter, & Silvia (under review).
#' 
#' @docType data
#' 
#' @usage data(hex)
#' 
#' @format A 802x16 response matrix
#' 
#' @keywords datasets
#' 
#' @references
#' 
#' Christensen, A.P., Cotter, K.N., Silvia, P.J. (under review).
#' Nomological network of Openness to Experience:
#' A network analysis of four Openness to Experience inventories.
#' \url{http://doi.org/10.17605/OSF.IO/954A7}
#' 
#' @examples 
#' 
#' data(hex)
"hex"
#----
#HEXACO Openness data----
#' HEXACO Openness to Experience Response Matrix (Binarized)
#' 
#' A response matrix (n = 802) of HEXACO's Openness to Experience
#' from Christensen, Cotter, & Silvia (under review).
#' 
#' @docType data
#' 
#' @usage data(hexb)
#' 
#' @format A 802x16 response matrix
#' 
#' @keywords datasets
#' 
#' @references
#' 
#' Christensen, A.P., Cotter, K.N., Silvia, P.J. (under review).
#' Nomological network of Openness to Experience:
#' A network analysis of four Openness to Experience inventories.
#' \url{http://doi.org/10.17605/OSF.IO/954A7}
#' 
#' @examples 
#' 
#' data(hexb)
"hexb"
