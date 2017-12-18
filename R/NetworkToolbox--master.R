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
    gij<-matrix(nrow=1,ncol=272)
    v<-matrix(nrow=1,ncol=272)
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
    diag(x)<-0
    
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
#' @return Returns a sparse TMFG-filtered matrix of the inverse covariance
#' @examples
#' 
#' LoGonet<-LoGo(hex)
#' 
#' @references 
#' Barfuss, W., Massara, G. P., Di Matteo, T., & Aste, T. (2016).
#' Parsimonious modeling with information filtering networks.
#' \emph{Physical Review E}, \emph{94}(6), 062306.
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
diag(x)<-0
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
#' @description Filters the network based on an r-value or alpha
#' @param data Can be a dataset or a correlation matrix
#' @param binary Is dataset dichotomous? Defaults to FALSE. Set TRUE if dataset is dichotomous (tetrachoric correlations are computed)
#' @param thresh Sets threshold (defaults to \emph{r} = .10). Set to "alpha" to use an alpha value
#' @param a Defaults to .05. Only applied when thresh = "alpha"
#' @return Returns a filtered adjacency matrix
#' @examples
#' 
#' threshnet<-threshold(hex)
#' 
#' alphanet<-threshold(hex, thresh = "alpha")
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Threshold filtering----
threshold <- function (data, binary = FALSE, thresh = .10, a = .05)
{
    if(nrow(data)==ncol(data)){cormat<-data}else
        if(binary){cormat<-psych::tetrachoric(data)$rho}else{cormat<-cor(data)}
    
    if(thresh=="alpha")
    {
        critical.r <- function(nrow, a){
            df <- nrow - 2
            critical.t <- qt( a/2, df, lower.tail = F )
            cvr <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
            return(cvr)}
        
        thresh<-critical.r(nrow(data),a)
    }
    
    cormat<-ifelse(cormat>=thresh,cormat,0)
    
    return(cormat)
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
#' Epskamp, S., Cramer, A. O., Waldorp, L. J., Schmittmann, V. D., & Borsboom, D. (2012).
#' qgraph: Network visualizations of relationships in psychometric data.
#' \emph{Journal of Statistical Software}, \emph{48}(4), 1-18.
#' 
#' Opsahl, T., Agneessens, F., & Skvoretz, J. (2010).
#' Node centrality in weighted networks: Generalizing degree and shortest paths.
#' \emph{Social Networks}, \emph{32}(3), 245-251.
#' 
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
  BC<-round(as.matrix(colSums(DP),ncol=60),0)}else{BC<-qgraph::centrality_auto(A)$node.centrality[,1]
  BC<-round(as.matrix(BC,ncol=60),0)}
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
#' \dontrun{
#' A<-TMFG(hex)$A
#'
#' weighted_LC<-closeness(A)
#' 
#' unweighted_LC<-closeness(A,weighted=FALSE)
#' }
#' @references
#' Epskamp, S., Cramer, A. O., Waldorp, L. J., Schmittmann, V. D., & Borsboom, D. (2012).
#' qgraph: Network visualizations of relationships in psychometric data.
#' Journal of Statistical Software, 48(4), 1-18.
#' 
#' Opsahl, T., Agneessens, F., & Skvoretz, J. (2010).
#' Node centrality in weighted networks: Generalizing degree and shortest paths.
#' Social networks, 32(3), 245-251.
#' 
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
  if (!weighted)
  {D<-distance(A,weighted=FALSE)
  C<-matrix(0,ncol=ncol(D))
  for(i in 1:ncol(D))
  {
    C[i]<-1/sum(D[,i])
  }
  LC<-t(as.data.frame(C)*100)
  rownames(LC)<-colnames(A)
  colnames(LC)<-c("LCu")}else{LC<-qgraph::centrality_auto(A)$node.centrality[,2]*10}
  LC<-round(as.data.frame(LC),3)
  LC<-as.matrix(LC)
  return(LC)
}
#----
#' Degree
#' @description Computes degree of each node in a network
#' @param A An adjacency matrix of network data
#' @return A vector of degree values for each node in the network
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
#' @return A vector of strength values for each node in the network
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
#' nodeimp<-impact(A)
#'
#' @references 
#' Kenett, Y. N., Kenett, D. Y., Ben-Jacob, E., & Faust, M. (2011).
#' Global and local features of semantic networks: Evidence from the Hebrew mental lexicon.
#' \emph{PloS one}, \emph{6}(8), e23912.
#' 
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
  Deg<-degree(A)
  Str<-strength(A)
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
#' @description Computes distance matrix of the network (Weighted not coded)
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
  return(D)}else{print("Weighted not coded.")}
}
#----
#' Characteristic Path Lengths
#' @description Computes global average shortest path length (ASPL), local average shortest path length (ASPLi), eccentricity (ecc), and diameter (D) of a network (Weighted not coded)
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
  {D<-distance(A,weighted=FALSE)
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
  rownames(ecc)<-colnames(A)
  
  return(list(ASPL=aspl,ASPLi=aspli,ecc=ecc,diameter=d))}
  else{print("Weighted not coded.")}
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
#' 
#' @references
#' Blondel, V. D., Guillaume, J. L., Lambiotte, R., & Lefebvre, E. (2008).
#' Fast unfolding of communities in large networks. 
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}(10), P10008.
#'  
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, \emph{52}(3), 1059-1069.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Louvain Community Detection----
louvain <- function (A, gamma = 1, M0 = 1:ncol(A), method = "modularity")
{
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
#' @param iter Number of random networks to generate, which are used to calculate the mean random ASPL and CC
#' @return Returns a value of small-worldness
#' @examples
#' 
#' A<-TMFG(hex)$A
#'
#' swm <- smallworldness(A)
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Small-worldness Measure----
smallworldness <- function (A, iter = 100)
{
    mat<-matrix(0,nrow=nrow(A),ncol=ncol(A)) #Initialize bootstrap matrix
    asamps<-matrix(0,nrow=iter) #Initialize sample matrix
    csamps<-matrix(0,nrow=iter) #Initialize sample matrix
    pb <- txtProgressBar(max=iter, style = 3)
    for(i in 1:iter) #Generate array of bootstrapped samples
    {
        f<-round(runif(i,min=1,max=1000000),0)
        set.seed(f[round(runif(i,min=1,max=length(f)),0)])
        rand<-randnet(ncol(A),sum(ifelse(A!=0,1,0))/2)
        asamps[i,]<-pathlengths(rand)$ASPL
        csamps[i,]<-clustcoeff(rand)$CC
        setTxtProgressBar(pb, i)
    }
    close(pb)
    nodes<-ncol(A)
    edges<-sum(ifelse(A!=0,1,0))/2
    rand<-randnet(nodes,edges)
    #rASPL<-pathlengths(rand)$ASPL
    rASPL<-mean(asamps)
    ASPL<-pathlengths(A)$ASPL
    CC<-clustcoeff(A)$CC
    
    #if(method=="rand")
    #{
    #rCC<-clustcoeff(rand)$CC
    rCC<-mean(csamps)
    swm<-(CC/rCC)/(ASPL/rASPL)
    return(list(S=swm,ASPL=ASPL,randASPL=rASPL,CC=CC,randCC=rCC))
    #}
    
    #if(method=="HG")
    #{lCC<-clustcoeff(rand)$CC
    #swm<-(rASPL/ASPL)-(CC/lCC)
    #return(list(S=swm,ASPL=ASPL,randASPL=rASPL,CC=CC,latCC=lCC))}
}
#----
#' Edge Replication
#' @description Computes the number of edges that replicate between two cross-sectional networks
#' @param A An adjacency matrix of network A
#' @param B An adjacency matrix of network B
#' @return Returns a list of the number of edges that replicate (Replicated), total number of edges (Possible), the percentage of edges that replicate (Percentage), the density of edges (Density), the mean difference between edges that replicate (MeanDifference), the sd of the difference between edges that replicate (SdDifference), and the correlation between the edges that replicate for both networks (Correlation)
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
  count<-0
if(!isSymmetric(A))
{A<-A+t(A)
warning("Adjacency matrix A was made to be symmetric")}
if(!isSymmetric(B))
{B<-B+t(B)
warning("Adjacency matrix B was made to be symmetric")}
  
  for(i in 1:ncol(A))
    for(j in 1:nrow(A))
      if(A[i,j]&&B[i,j]!=0){count<-count+1}
      count<-count/2
  
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
#' @description Computes the average and standard deviation of the
#' @param A An adjacency matrix of network A
#' @return Returns a list of the edge weights (Weights), the mean (Mean), the standard deviation (SD), and the sum of the edge weights (Total) in the network
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
#' Bootstrapped Network Preprocessing
#' @description Bootstraps the sample to identify the most stable correlations
#' @param data A set of data
#' @param method A network filtering method (e.g, "TMFG", "MaST", "ECO", "ECOplusMaST")
#' @param binary Is dataset dichotomous? Defaults to FALSE. Set TRUE if dataset is dichotomous (tetrachoric correlations are computed)
#' @param n Number of people to use in the bootstrap. Defaults to full sample size
#' @param iter Number of bootstrap iterations. Defaults to 1000 iterations
#' @param a Alpha to be used for determining the critical value of correlation coefficients. Defaults to .05
#' @param depend Is network a dependency (or directed) network? Defaults to FALSE. Set TRUE to generate a TMFG-filtered dependency network
#' @return Returns a list that includes the original filtered network (orignet),
#' correlation matrix of the mean bootstrapped network (bootmat),
#' reliabilities of the connections in the original network (netrel),
#' reliabilities of the connections in the bootstrapped network (bootrel),
#' a plot of the bootrel reliability matrix (netrel; upper triangle = actual network reliabilites, bootrel; lower triangle = overall network reliablities),
#' a plot of included correlations on their reliability (ConR)
#' @examples
#' \dontrun{
#' 
#' prepTMFG<-prepboot(hex,method="TMFG")
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
#' @importFrom stats lm na.omit
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#Network Preprocessing Bootstrap----
prepboot <- function (data, method, binary = FALSE, n = nrow(data), iter = 1000, a = .05, depend = FALSE)
{
    if(nrow(data)==ncol(data)){stop("Input must be a dataset")}else
        if(binary){realmat<-psych::tetrachoric(data)$rho}else{realmat<-cor(data)}
    mat<-matrix(0,nrow=n,ncol=ncol(data)) #Initialize bootstrap matrix
    samps<-array(0,c(nrow=nrow(realmat),ncol=ncol(realmat),iter)) #Initialize sample matrix
    dsamps<-array(0,c(nrow=nrow(realmat),ncol=ncol(realmat),iter)) #Initialize dsample matrix
    pb <- txtProgressBar(max=iter, style = 3)
    for(i in 1:iter) #Generate array of bootstrapped samples
    {
        f<-round(runif(i,min=1,max=1000000),0)
        set.seed(f[round(runif(i,min=1,max=length(f)),0)])
        mat<-data[round(runif(n,min=1,max=n),0),]
        if(any(colSums(mat)<=1)){stop("Increase sample size: not enough observations")}
        cormat<-cor(mat)
        if(!depend){
        if(method=="TMFG")
        {samps[,,i]<-TMFG(cormat)$A
        tru<-TMFG(data)$A
        }else if(method=="MaST")
        {samps[,,i]<-MaST(cormat)
        tru<-MaST(data)
        }else if(method=="ECOplusMaST")
        {samps[,,i]<-ECOplusMaST(cormat)
        tru<-ECOplusMaST(data)
        }else if(method=="ECO")
        {samps[,,i]<-ECO(cormat)
        tru<-ECO(data)
        }else stop("Method not available")
        }else{if(method=="TMFG")
        {samps[,,i]<-TMFG(cormat)$A
        dsamps[,,i]<-TMFG(depend(cormat,progBar=FALSE),depend=TRUE)$A
        }else stop("Method not available")
        }
        setTxtProgressBar(pb, i)
    }
    close(pb)
    
    
    tru<-TMFG(data)$A
    
    #Mean matrix
    meanmat<-matrix(0,nrow=nrow(realmat),ncol=ncol(realmat)) #Initialize Mean matrix
    for(j in 1:nrow(realmat))
        for(k in 1:ncol(realmat))
        {meanmat[j,k]<-mean(samps[j,k,])}

    #Remove non-significant edges
        critical.r <- function(iter, a = .05){
            df <- iter - 2
            critical.t <- qt( a/2, df, lower.tail = F )
            cvr <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
            return(cvr)}
    
    for(x in 1:nrow(meanmat))
        for(y in 1:ncol(meanmat))
            if(meanmat[x,y]<=critical.r(iter))
            {meanmat[x,y]<-0}
    
    #return meanmat to bootmat
    bootmat<-meanmat
    colnames(bootmat)<-colnames(bootmat)
    
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
            xo<-na.omit(x)
            yo<-na.omit(y)
            
            mar=c(2,2,2,2)
            cpo<-{plot(xo,yo,pch=16,ylab="Correlation Strength",xlab="Reliability",
                       main="Correlation Strength on Reliability",xlim=c(0,1),ylim=range(yo))
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
                
                #critical value
                dcritical.r <- function(iter, a = .05){
                    df <- iter - (3 + (ncol(realmat)-1))
                    critical.t <- qt( a/2, df, lower.tail = F )
                    cvr <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
                    return(cvr)}
                    
                    for(x in 1:nrow(dmeanmat))
                        for(y in 1:ncol(dmeanmat))
                            if(dmeanmat[x,y]<=dcritical.r(iter))
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
    
    if(!depend)
    {return(list(orignet=tru,bootmat=bootmat,netrel=upp,bootrel=reprel,plotrel=plt,ConR=cpo))
    }else{orignet<-(list(orignet=tru,bootmat=bootmat,bootrel=reprel,plotrel=plt,ConR=cpo))
          depnet<-(list(orignet=dtru,bootmat=dbootmat,bootrel=dreprel,plotrel=dplt,DonR=dcpo))
          return(list(undirected=orignet,directed=depnet))}
}
#----
#' Bootstrapped Walktrap Communities Likelihood
#' @description Bootstraps the sample with replace to compute walktrap reliability (TMFG-filtered networks only)
#' @param data A set of data
#' @param binary Is dataset dichotomous? Defaults to FALSE. Set TRUE if dataset is dichotomous (tetrachoric correlations are computed)
#' @param n Number of people to use in the bootstrap. Defaults to full sample size
#' @param iter Number of bootstrap iterations. Defaults to 1000 iterations
#' @param steps Number of steps to use in the walktrap algorithm. Defaults to 4. Use a larger number of steps for smaller networks
#' @return The factors and their proportion found across bootstrapped samples (i.e., their likelihood)
#' @examples
#' \dontrun{
#' 
#' walkTMFG<-walkboot(hex)
#' }
#' @references
#' Csardi, G., & Nepusz, T. (2006).
#' The igraph software package for complex network research.
#' \emph{InterJournal, Complex Systems}, \emph{1695}(5), 1-9.
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Bootstrapped Walktrap Reliability----
walkboot <- function (data, binary = FALSE, n = nrow(data), iter = 1000, steps = 4)
{
    col<-ncol(data)
    if(nrow(data)==ncol(data)){stop("Input must be a dataset")}else
        if(binary){realmat<-psych::tetrachoric(data)$rho}else{realmat<-cor(data)}
    mat<-matrix(0,nrow=n,ncol=col) #Initialize bootstrap matrix
    walk<-matrix(0,nrow=iter,ncol=1) #Initialize walktrap matrix
    pb <- txtProgressBar(max=iter, style = 3)
    for(i in 1:iter) #Generate array of bootstrapped samples
    {
        f<-round(runif(i,min=1,max=1000),0)
        set.seed(f[round(runif(i,min=1,max=length(f)),0)])
        mat<-data[round(runif(n,min=1,max=n),0),]
        if(any(colSums(mat)<=1)){stop("Increase sample size: not enough observations")}
        cormat<-cor(mat)
        walk[i,]<-max(igraph::walktrap.community(igraph::as.igraph(qgraph::qgraph(TMFG(cormat)$A,DoNotPlot=TRUE)),steps=steps)$membership)
        setTxtProgressBar(pb, i)
    }
    close(pb)
    
    count<-0
    prop<-matrix(0,nrow=length(seq(from=min(walk),to=max(walk))),ncol=1)
    for(i in min(walk):max(walk))
    {
        count<-count+1
        prop[count,]<-length(which(walk==i))
    }
    
    prop<-round(prop/iter,3)
    prop<-cbind(seq(from=min(walk),to=max(walk)),prop)
    colnames(prop)<-c("Factors","Likelihood")
    
    return(prop)
}
#----
#' Dependency Matrix
#' @description Generates a dependency matrix of the data
#' @param data A set of data
#' @param binary Is dataset dichotomous? Defaults to FALSE. Set TRUE if dataset is dichotomous (tetrachoric correlations are computed)
#' @param progBar Should progress bar be displayed? Defaults to TRUE. Set FALSE for no progress bar.
#' @return Returns an adjacency matrix of dependencies
#' @examples
#' D<-depend(hex)
#' 
#' binaryD<-depend(hexb,binary=TRUE)
#' 
#' @references
#' Kenett, D. Y., Tumminello, M., Madi, A., Gur-Gershgoren, G., Mantegna, R. N., & Ben-Jacob, E. (2010).
#' Dominating clasp of the financial sector revealed by partial correlation analysis of the stock market.
#' \emph{PloS one}, \emph{5}(12), e15032.
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Dependency----
depend <- function (data, binary = FALSE, progBar = TRUE)
{
    if(nrow(data)==ncol(data)){cormat<-data}else
        if(binary){cormat<-psych::tetrachoric(data)$rho}else{cormat<-cor(data)}
    
    inter<-((ncol(cormat)*(ncol(cormat)-1)*(ncol(cormat)-2)))
    
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
    
    for(h in 1:j)
    diag(parmat[,,h])<-1

    depmat<-matrix(0,nrow=nrow(parmat),ncol=ncol(parmat))
    for(i in 1:ncol(parmat))
        for(j in 1:ncol(parmat))
        {depmat[j,i]<-mean(parmat[i,-j,j])}
    
    colnames(depmat)<-colnames(data)
    return(depmat)
}
#----
#' Random Network
#' @description Generates a random network
#' @param nodes Number of nodes in random network
#' @param edges Number of edges in random network
#' @return Returns an adjacency matrix of a random network
#' @examples
#' 
#' rand <- randnet(10,27)
#' 
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
#' Binarize Network
#' @description Converts weighted adjacency matrix to a binarized adjacency matrix
#' @param A An adjacency matrix of network data
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
