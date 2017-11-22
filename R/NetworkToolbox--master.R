#  NetworkToolbox Functions
#
#' Triangulated Maximally Filtered Graph
#' @description Applies the Triangulated Maximally Filtered Graph (TMFG) filtering method
#' @param data Can be a dataset or a correlation matrix
#' @param binary Is dataset dichotomous? Defaults to FALSE. Set TRUE if dataset is dichotomous (tetrachoric correlations are computed)
#' @param weighted Should network be weighted? Defaults to TRUE. Set FALSE to produce an unweighted (binary) network
#' @return Returns a list of the adjacency matrix (A) and separators (sep)
#' @examples
#' weighted_TMFGnetwork<-TMFG(data)
#' 
#' weighted_binary_TMFGnetwork<-TMFG(data,binary=TRUE)
#' 
#' unweighted_TMFGnetwork<-TMFG(data,weighted=FALSE)
#' 
#' unweighted_binary_TMFGnetwork<-TMFG(data,binary=TRUE,weighted=FALSE)
#' @references 
#' Massara, G. P., Di Matteo, T., & Aste, T. (2016).
#' Network filtering for big data: Triangulated maximally filtered graph.
#' \emph{Journal of Complex Networks}, 5(2), 161-178.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @importFrom stats cor sd runif qt
#' @export
#TMFG Filtering Method----
TMFG <-function (data,binary=FALSE,weighted=TRUE)
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
  s<-colSums(cormat*(cormat>mean(matrix(unlist(cormat),nrow=1)))*1)
  in_v[1:4]<-order(s,decreasing=TRUE)[1:4]
  ou_v<-setdiff(1:nrow(in_v),in_v)
  #build tetrahedron with the largest strength
  tri[1,]<-in_v[1:3,]
  tri[2,]<-in_v[2:4,]
  tri[3,]<-in_v[c(1,2,4),]
  tri[4,]<-in_v[c(1,3,4),]
  S<-matrix(nrow=(3*nrow(cormat)-6),ncol=3) #initialize sparse matrix
  S[1,]<-c(in_v[1],in_v[2],1)
  S[2,]<-c(in_v[1],in_v[3],1)
  S[3,]<-c(in_v[1],in_v[4],1)
  S[4,]<-c(in_v[2],in_v[3],1)
  S[5,]<-c(in_v[2],in_v[4],1)
  S[6,]<-c(in_v[3],in_v[4],1)
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
  L<-S
  L[,1]<-S[,2]
  L[,2]<-S[,1]
  K<-rbind(S,L)
  x<-as.matrix(Matrix::sparseMatrix(i=K[,1],j=K[,2],x=K[,3]))
  diag(x)<-0
  if(weighted)
   {
    for(r in 1:nrow(x))
      for(z in 1:ncol(x))
      {if(x[r,z]==1){x[r,z]<-cormat[r,z]}
   }
  }else

  x<-as.matrix(x)
  colnames(x)<-colnames(cormat)
  return(list(A=x, sep=separators))
}
#----
#' Maximum Spanning Tree
#' @description Applies the Maximum Spanning Tree (MaST) filtering method
#' @param data Can be a dataset or a correlation matrix
#' @param binary Is dataset dichotomous? Defaults to FALSE. Set TRUE if dataset is dichotomous (tetrachoric correlations are computed)
#' @param weighted Should network be weighted? Defaults to TRUE. Set FALSE to produce an unweighted (binary) network
#' @return A sparse association matrix
#' @examples
#' weighted_MaSTnetwork<-MaST(data)
#' 
#' weighted_binary_MaSTnetwork<-MaST(data,binary=TRUE)
#' 
#' unweighted_MaSTnetwork<-MaST(data,weighted=FALSE)
#' 
#' unweighted_binary_MaSTnetwork<-MaST(data,binary=TRUE,weighted=FALSE)
#' @references 
#' Adapted from: \url{https://www.mathworks.com/matlabcentral/fileexchange/23276-maximum-weight-spanning-tree--undirected}
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Maximum Spanning Tree----
MaST <- function (data,binary=FALSE,weighted=TRUE)
{FIND_PathCompression <- function (temproot=temproot)
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
L[,2]<-S[,3]
L[,3]<-S[,2]
K<-rbind(S,L)
K<-cbind(K,K[,1])
K<-K[,-1]
x<-as.matrix(Matrix::sparseMatrix(i=K[,1],j=K[,2],x=K[,3]))
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
#' weighted_undirected_ECOnetwork<-ECO(data)
#' 
#' unweighted_undirected_ECOnetwork<-ECO(data,weighted=FALSE)
#' 
#' weighted_directed_ECOnetwork<-ECO(data,directed=TRUE)
#' 
#' unweighted_directed_ECOnetwork<-ECO(data,weighted=FALSE,directed=TRUE)
#' 
#' weighted_undirected_binary_ECOnetwork<-ECO(data,binary=TRUE)
#' 
#' unweighted_undirected_binary_ECOnetwork<-ECO(data,weighted=FALSE,binary=TRUE)
#' 
#' weighted_directed_binary_ECOnetwork<-ECO(data,directed=TRUE,binary=TRUE)
#' 
#' unweighted_directed_binary_ECOnetwork<-ECO(data,weighted=FALSE,directed=TRUE,binary=TRUE)
#' @references 
#' Fallani, F. D. V., Latora, V., & Chavez, M. (2017).
#' A topological criterion for filtering information in complex brain networks.
#' \emph{PLoS Computational Biology}, 13(1), e1005305.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#ECO Neural Network Filter----
ECO <- function (data, weighted=TRUE, binary=FALSE, directed=FALSE)
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
#' weighted_ECOplusMaSTnetwork<-ECOplusMaST(data)
#' 
#' weighted_binary_ECOplusMaSTnetwork<-ECOplusMaST(data,binary=TRUE)
#' 
#' unweighted_binary_ECOplusMaSTnetwork<-ECOplusMaST(data,weighted=FALSE)
#' 
#' unweighted_binary_ECOplusMaSTnetwork<-ECOplusMaST(data,binary=TRUE,weighted=FALSE)
#' @references 
#' Fallani, F. D. V., Latora, V., & Chavez, M. (2017).
#' A topological criterion for filtering information in complex brain networks.
#' \emph{PLoS Computational Biology}, 13(1), e1005305.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#ECO Filter + MaST----
ECOplusMaST <- function (data, weighted=TRUE, binary=FALSE)
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
#' Betwenness Centrality
#' @description Computes betweenness centrlaity of each node in a network
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted? Defaults to TRUE. Set to FALSE for unweighted measure of betwenness centrality
#' @return A vector of betweenness centrality values for each node in the network
#' @examples
#' weighted_BC<-betweenness(A)
#' 
#' unweighted_BC<-betweenness(A,weighted=FALSE)
#' @references 
#' Epskamp, S., Cramer, A. O., Waldorp, L. J., Schmittmann, V. D., & Borsboom, D. (2012).
#' qgraph: Network visualizations of relationships in psychometric data.
#' \emph{Journal of Statistical Software}, 48(4), 1-18.
#' 
#' Opsahl, T., Agneessens, F., & Skvoretz, J. (2010).
#' Node centrality in weighted networks: Generalizing degree and shortest paths.
#' \emph{Social Networks}, 32(3), 245-251.
#' 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, 52(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Betweenness Centrality----
betweenness <- function (A,weighted=TRUE)
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
  BC<-BC-(ncol(A)-1)
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
#' rspbc<-rspbc(A, beta=0.01)
#' @references 
#' Kivimaki, I., Lebichot, B., Saramaki, J., & Saerens, M. (2016).
#' Two betweenness centrality measures based on Randomized Shortest Paths.
#' \emph{Scientific Reports}, 6 (19668), 1-15.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Randomized Shortest Paths Betweennesss Centrality----
rspbc <- function (A, beta=0.01)
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
#' \emph{Neuroimage}, 52(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export 
#Closeness Centrality----
closeness <- function (A,weighted=TRUE)
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
  colnames(LC)<-c("LCu")}else{LC<-qgraph::centrality_auto(A)$node.centrality[,2]*100}
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
#' deg<-degree(A)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, 52(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Degree----
degree <- function (A)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  A<-ifelse(A!=0,1,0)
  Deg<-as.data.frame(colSums(A))
  colnames(Deg)<-c("Degree")
  Deg<-as.matrix(Deg)
  return(Deg)
}
#----
#' Node Strength
#' @description Computes strength of each node in a network
#' @param A An adjacency matrix of network data
#' @return A vector of strength values for each node in the network
#' @examples
#' str<-strength(A)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, 52(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Node Strength----
strength <- function (A)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  strength<-round(as.data.frame(colSums(A)),2)
  colnames(strength)<-c("Strength")
  strength<-as.matrix(strength)
  return(strength)
}
#----
#' Eigenvector Centrality
#' @description Computes eigenvector centrality of each node in a network
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted? Defaults to TRUE. Set to FALSE for unweighted measure of eigenvector centrality
#' @return A vector of eigenvector centrality values for each node in the network
#' @examples
#' weighted_EC<-eigenvector(A)
#' 
#' unweighted_EC<-eigenvector(A,weighted=FALSE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, 52(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Eigenvector----
eigenvector <- function (A,weighted=TRUE)
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
#' weighted_lev<-leverage(A)
#'
#' unweighted_lev<-leverage(A, weighted=FALSE)
#' @references 
#' Joyce, K. E., Laurienti, P. J., Burdette, J. H., & Hayasaka, S. (2010).
#' A new measure of centrality for brain networks. 
#' \emph{PLoS One}, 5(8), e12200.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export 
#Leverage Centrality
leverage <- function (A, weighted=TRUE)
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
#' Hybrid Centrality
#' @description Computes hybrid centrality of each node in a network
#' @param A An adjacency matrix of network data
#' @return A vector of hybrid centrality values for each node in the network (higher values are more central, lower values are more peripheral)
#' @examples
#' \dontrun{
#' 
#' HC<-hybrid(A)
#' }
#' @references 
#' Pozzi, F., Di Matteo, T., & Aste, T. (2013).
#' Spread of risk across financial markets: Better to invest in the peripheries. 
#' \emph{Scientific Reports}, 3(1655), 1-7.
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
#' \dontrun{
#' 
#' weighted_centralitylist<-centlist(A)
#' 
#' unweighted_centralitylist<-centlist(A,weighted=FALSE)
#' }
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, 52(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Centrality List----
centlist <- function (A, weighted=TRUE)
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
    return(list(Betweenness=BC,Closeness=CC,Strength=Str,Eigenvector=EC,Leverage=lev))}
}
#----
#' Distance
#' @description Computes distance matrix of the network (Weighted not coded)
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted? Defaults to FALSE. Set to TRUE for weighted measure of distance
#' @return A distance matrix of the network
#' @examples
#' unweighted_D<-distance(A)
#' 
#' weighted_D<-distance(A,weighted=TRUE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage, 52(3)}, 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Distance----
distance<-function (A,weighted=FALSE)
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
#' unweighted_PL<-pathlengths(A)
#' 
#' weighted_PL<-pathlengths(A,weighted=TRUE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, 52(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Path Lengths----
pathlengths <- function (A, weighted=FALSE)
{
  if(nrow(A)!=ncol(A))
  {stop("Input not an adjacency matrix")}
  if(!weighted)
  {D<-distance(A,weighted=FALSE)
  n<-nrow(D)
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
  
  return(list(ASPL=aspl,ASPLi=aspli,Eccentricity=ecc,Diameter=d))}
  else{print("Weighted not coded.")}
}
#----
#' Clustering Coefficient
#' @description Computes global clustering coefficient (CC) and local clustering coefficient (CCi)
#' @param A An adjacency matrix of network data
#' @param weighted Is the network weighted? Defaults to FALSE. Set to TRUE for weighted measures of CC and CCi
#' @return Returns a list of CC and CCi
#' @examples
#' unweighted_CC<-clustcoeff(A)
#' 
#' weighted_CC<-clustcoeff(A,weighted=TRUE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' \emph{Neuroimage}, 52(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Clustering Coefficient----
clustcoeff <- function (A, weighted=FALSE)
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
#' Edge Replication
#' @description Computes the number of edges that replicate between two cross-sectional networks
#' @param A An adjacency matrix of network A
#' @param B An adjacency matrix of network B
#' @return Returns a list of the number of edges that replicate (Replicated), total number of edges (Possible), the percentage of edges that replicate (Percentage), the density of edges (Density), the mean difference between edges that replicate (MeanDifference), the sd of the difference between edges that replicate (SdDifference), and the correlation between the edges that replicate for both networks (Correlation)
#' @examples
#' edges<-edgerep(A,B)
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Edge Replication----
edgerep <- function (A, B)
{
  count<-0
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
  
  return(list(Replicated=count,
              TotalEdgesA=possibleA,TotalEdgesB=possibleB,
              PercentageA=percentA,PercentageB=percentB,
              DensityA=densityA,DensityB=densityB,
              MeanDifference=mvec,SdDifference=svec,Correlation=corr))
}
#----
#' Network Connectivity
#' @description Computes the average and standard deviation of the
#' @param A An adjacency matrix of network A
#' @return Returns a list of the edge weights (Weights), the mean (Mean), the standard deviation (SD), and the sum of the edge weights (Total) in the network
#' @examples
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
    
    return(list(Weights=weights,Mean=mea,SD=s,Total=tot))
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
#' @return Returns a list that includes a correlation matrix of the mean bootstrapped network (bootmat), reliabilities of the connections in the network (bootrel), and a plot of the bootrel reliability matrix (plotrel; upper triangle = actual network reliabilites, lower triangle = overall network reliablities)
#' @examples
#' \dontrun{
#' 
#' prepTMFG<-prepboot(data,method="TMFG")
#' }
#' @references
#' Tumminello, M., Coronnello, C., Lillo, F., Micciche, S., & Mantegna, R. N. (2007).
#' Spanning trees and bootstrap reliability estimation in correlation-based networks.
#' \emph{International Journal of Bifurcation and Chaos}, 17(7), 2319-2329.
#' 
#' Wei, T. & Simko, V.(2017).
#' R package "corrplot": Visualization of a correlation matrix (Version 0.84).
#' Available from \url{https://github.com/taiyun/corrplot}
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Network Preprocessing Bootstrap----
prepboot <- function (data, method, binary = FALSE, n = nrow(data), iter = 1000, a = .05)
{
    if(nrow(data)==ncol(data)){stop("Input must be a dataset")}else
        if(binary){realmat<-psych::tetrachoric(data)$rho}else{realmat<-cor(data)}
    mat<-matrix(0,nrow=n,ncol=ncol(data)) #Initialize bootstrap matrix
    samps<-array(0,c(nrow=nrow(realmat),ncol=ncol(realmat),iter)) #Initialize sample matrix
    for(i in 1:iter) #Generate array of bootstrapped samples
    {
        f<-round(runif(i,min=1,max=1000),0)
        set.seed(f[round(runif(i,min=1,max=length(f)),0)])
        mat<-data[round(runif(n,min=1,max=n),0),]
        if(any(colSums(mat)<=1)){stop("Increase sample size: not enough observations")}
        cormat<-cor(mat)
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
    }
    
    #Mean matrix
    meanmat<-matrix(0,nrow=nrow(realmat),ncol=ncol(realmat)) #Initialize Mean matrix
    for(j in 1:nrow(realmat))
        for(k in 1:ncol(realmat))
        {meanmat[j,k]<-mean(samps[j,k,])}
    
    #Reliability matrix
    samp<-array(0,c(nrow=nrow(realmat),ncol=ncol(realmat),iter))
    rel<-matrix(0,nrow=nrow(realmat),ncol=ncol(realmat))
    for(j in 1:nrow(realmat))
        for(k in 1:ncol(realmat))
            for(l in 1:iter)
            if(samps[j,k,l]!=0)
            {samp[j,k,l]<-1}
    for(j in 1:nrow(realmat))
        for(k in 1:ncol(realmat))
            rel[j,k]<-sum(samp[j,k,])/iter
            colnames(rel)<-colnames(data)
            #reliablity plot
            upp<-matrix(0,nrow=nrow(rel),ncol=ncol(rel))
            for(i in 1:nrow(rel))
                for(j in 1:ncol(rel))
                    if(rel[i,j]!=0&&tru[i,j]!=0)
                    {upp[i,j]<-rel[i,j]}
            reprel<-rel
            colnames(upp)<-colnames(rel)
            rel[upper.tri(rel)]<-upp[upper.tri(upp)]
            row.names(rel)<-colnames(rel)
            plt<-corrplot::corrplot(rel,method="color",
            title="Bootstrapped Correlation Reliabilities",
            mar=c(2,2,2,2),tl.col="black",tl.cex=.75,
            cl.lim = c(0,1),addgrid.col = "grey")
            #reliablity on correlation plot
            x<-matrix(nrow=length(upp))
            y<-matrix(nrow=length(tru))
            wc<-0
            for(i in 1:60)
                for(j in 1:60)
                    if((upp[i,j]&&tru[i,j])!=0)
                    {wc<-wc+1
                    x[wc]<-upp[i,j]
                    y[wc]<-tru[i,j]}
            
            xo<-na.omit(x)
            yo<-na.omit(y)
            
            cpo<-{plot(xo,yo,pch=16,ylab="Correlation Strength",xlab="Reliability",
                 main="Reliability on Correlation Strength",xlim=c(0,1),ylim=range(yo))
            abline(lm(yo~xo))
            text(x=.2,y=(min(yo)+(max(range(yo))/2)),labels = paste("r = ",round(cor(yo,xo),3)))}
    
    j<-meanmat
    critical.r <- function(iter, a = .05){
        df <- iter - 2
        critical.t <- qt( a/2, df, lower.tail = F )
        cvr <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
        return(cvr)
    }
    
    for(x in 1:nrow(meanmat))
        for(y in 1:ncol(meanmat))
            if(meanmat[x,y]<=critical.r(iter))
            {
                meanmat[x,y]<-0
            }
    j<-meanmat
    j<-as.data.frame(j)
    colnames(j)<-colnames(data)
    j<-as.matrix(j)
    
    return(list(bootmat=j,bootrel=reprel,plotrel=plt,netrel=upp))
}
#----
#HEXACO Openness data----
#' HEXACO Openness to Experience Correlation Matrix
#' 
#' Correlation matrix (n = 802) of HEXACO's Openness to Experience
#' from Christensen, Cotter, & Silvia (in preparation).
#' Nomological Network of Openness to Experience.
#' 
#' @docType data
#' 
#' @usage data(data)
#' 
#' @format A 16x16 correlation matrix
#' 
#' @keywords datasets
#' 
#' @references
#' 
#' Christensen, A.P., Cotter, K.N., Silvia, P.J. (in preparation).
#' Nomological network of openness to experience:
#' A network analysis of four openness to experience inventories.
#' 
#' @examples 
#' 
#' data(data)
"data"
#----
#HECACO Openness to Experience TMFG Adjacency matrix----
#' HEXACO Openness to Experience TMFG Adjacency Matrix
#' 
#' TMFG filtered association matrix (n = 802) of HEXACO's Openness to Experience
#' from Christensen, Cotter, & Silvia (in preparation).
#' Nomological Network of Openness to Experience.
#' 
#' @docType data
#' 
#' @usage data(A)
#' 
#' @format A 16x16 TMFG filtered adjacency matrix
#' 
#' @keywords datasets
#' 
#' @references
#' 
#' Christensen, A.P., Cotter, K.N., Silvia, P.J. (in preparation).
#' Nomological network of openness to experience:
#' A network analysis of four openness to experience inventories.
#' 
#' @examples 
#' 
#' data(A)
"A"
#----
#HECACO Openness to Experience MaST Adjacency matrix----
#' HEXACO Openness to Experience MaST Adjacency Matrix
#' 
#' MaST filtered association matrix (n = 802) of HEXACO's Openness to Experience
#' from Christensen, Cotter, & Silvia (in preparation).
#' Nomological Network of Openness to Experience.
#' 
#' @docType data
#' 
#' @usage data(B)
#' 
#' @format A 16x16 TMFG filtered adjacency matrix
#' 
#' @keywords datasets
#' 
#' @references
#' 
#' Christensen, A.P., Cotter, K.N., Silvia, P.J. (in preparation).
#' Nomological network of openness to experience:
#' A network analysis of four openness to experience inventories.
#' 
#' @examples 
#' 
#' data(B)
"B"
#----