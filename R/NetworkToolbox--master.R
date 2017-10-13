#Network Toolbox Functions

#' Triangulated Maximally Filtered Graph
#' @description Applies the Triangulated Maximally Filtered Graph (TMFG) filtering method.
#' @param data Can be a dataset or a correlation matrix
#' @param binary Is dataset dichotomous? Defaults to FALSE. Set TRUE if dataset is dichotomous but could be on a continuous, normal distribution.
#' @param weighted Should network be weighted? Defaults to TRUE. Set FALSE to produce an unweighted (binary) network.
#' @return A sparse association matrix
#' @examples
#' weighted_TMFGnetwork<-TMFG(data)
#' weighted_binarydata_TMFGnetwork<-TMFG(data,binary=TRUE)
#' unweighted_TMFGnetwork<-TMFG(data,weighted=FALSE)
#' unweighted_binarydata_TMFGnetwork<-TMFG(data,binary=TRUE,weighted=FALSE)
#' @references 
#' Massara, G. P., Di Matteo, T., & Aste, T. (2016).
#' Network filtering for big data: Triangulated maximally filtered graph.
#' Journal of Complex Networks, 5(2), 161-178.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#TMFG Filtering Method----
TMFG <-function (data=data,binary=FALSE,weighted=TRUE)
{
  if(nrow(data)==ncol(data)){cormat<-data}else
    if(binary){cormat<-psych::tetrachoric(data)$rho}else{cormat<-cor(data)}
  n<-ncol(cormat)
  if(n<9){print("Matrix is too small")}
  if(any(cormat<0)){print("Matrix has negative elements!")}
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
  s<-colSums(cormat*(cormat>mean(matrix(cormat,nrow=1)))*1)
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
      {
        if(x[r,z]==1)
        {
          x[r,z]<-cormat[r,z]
        }
      }
  }else
    x<-as.data.frame(x)
  colnames(x)<-colnames(cormat)
  round(x,3)
}
#----
#' Maximum Spanning Tree
#' @description Applies the Maximum Spanning Tree (MaST) filtering method.
#' @param data Can be a dataset or a correlation matrix
#' @param binary Is dataset dichotomous? Defaults to FALSE. Set TRUE if dataset is dichotomous but could be on a continuous, normal distribution.
#' @param weighted Should network be weighted? Defaults to TRUE. Set FALSE to produce an unweighted (binary) network.
#' @return A sparse association matrix
#' @examples
#' weighted_MaSTnetwork<-MaST(data)
#' weighted_binarydata_MaSTnetwork<-MaST(data,binary=TRUE)
#' unweighted_MaSTnetwork<-MaST(data,weighted=FALSE)
#' unweighted_binarydata_MaSTnetwork<-MaST(data,binary=TRUE,weighted=FALSE)
#' @references 
#' Adapted from: <https://www.mathworks.com/matlabcentral/fileexchange/23276-maximum-weight-spanning-tree--undirected>
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Maximum Spanning Tree----
MaST <- function (data=data,binary=FALSE,weighted=TRUE)
{FIND_PathCompression <- function (temproot=temproot)
{
  ParentPointer[temproot]
  if(ParentPointer[temproot]!=temproot)
  {ParentPointer[temproot]<-FIND_PathCompression(ParentPointer[temproot])}
  parent<-ParentPointer[temproot]
}
if(nrow(data)==ncol(data)){cormat<-data}else
  if(binary){cormat<-psych::tetrachoric(data)$rho}else{cormat<-cor(data)}
nodeT<-0
nodeF<-0
weights<-0
wc<-0
n<-ncol(cormat)
for (i in 1:n)
  for (j in 1:n)
    if (cormat[i,j]!=0) #Figure out how to remove warning
    {
      wc<- wc+1
      nodeT[wc] <- i
      nodeF[wc] <- j
      weights[wc] <- cormat[i,j]
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
x<-as.data.frame(x)
ifelse(x!=0,cormat,0)
colnames(x)<-colnames(cormat)
round(x,3)
if(!weighted)
{x<-ifelse(x!=0,1,0)}
round(as.matrix(x),3)
}
#----
#' ECO Neural Network Filter
#' @description Applies the ECO neural network filtering method.
#' @param data Can be a dataset or a correlation matrix
#' @param weighted Should network be weighted? Defaults to TRUE. Set FALSE to produce an unweighted (binary) network.
#' @param directed Is the network directed? Defaults to FALSE. Set TRUE if the network is directed.
#' @return A sparse association matrix
#' @examples
#' weighted_undirected_ECOnetwork<-ECO(data)
#' unweighted_undirected_ECOnetwork<-ECO(data,weighted=FALSE)
#' weighted_directed_ECOnetwork<-ECO(data,directed=TRUE)
#' unweighted_directed_ECOnetwork<-ECO(data,weighted=FALSE,directed=TRUE)
#' @references 
#' Fallani, F. D. V., Latora, V., & Chavez, M. (2017).
#' A topological criterion for filtering information in complex brain networks.
#' PLoS Computational Biology, 13(1), e1005305.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#ECO Neural Network Filter
ECO <- function (data=data, weighted=TRUE, directed=FALSE)
{
  if(nrow(data)==ncol(data)){C<-data}else{C<-cor(data)}
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
  {W}else{W<-C+t(C)
  diag(W)<-1}
  J<-S+t(S)
  diag(J)<-1
  if(weighted)
  {
    W<-ifelse(W!=0,J,0)
  }
  round(W,3)
}
#----
#' ECO+MaST Network Filter
#' @description Applies the ECO neural network filtering combined with the MaST filtering method.
#' @param data Can be a dataset or a correlation matrix
#' @param weighted Should network be weighted? Defaults to TRUE. Set FALSE to produce an unweighted (binary) network.
#' @return A sparse association matrix
#' @examples
#' weighted_ECOplusMaSTnetwork<-ECO(data)
#' unweighted_ECOplusMaSTnetwork<-ECO(data,weighted=FALSE)
#' @references 
#' Fallani, F. D. V., Latora, V., & Chavez, M. (2017).
#' A topological criterion for filtering information in complex brain networks.
#' PLoS Computational Biology, 13(1), e1005305.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#ECO Filter + MaST
ECOplusMaST <- function (data=data, weighted=TRUE)
{if(weighted)
{
  a<-MaST(data)
  b<-ECO(data)
  k<-matrix(NA,nrow=nrow(a),ncol(a))
  for(i in 1:nrow(a))
    for(j in 1:ncol(a))
      if(a[i,j]==b[i,j]){k[i,j]<-a[i,j]}else k[i,j]<-b[i,j]
}else{
  a<-MaST(data,weighted=FALSE)
  b<-ECO(data,weighted=FALSE)
  k<-matrix(NA,nrow=nrow(a),ncol(a))
  for(i in 1:nrow(a))
    for(j in 1:ncol(a))
      if(a[i,j]==b[i,j])
      {k[i,j]<-a[i,j]}else k[i,j]<-b[i,j]}
  k
}
#----
#' Semantic Network Cleaner
#' @description An automated cleaning function for semantic network data (still in the alpha stage).
#' @param data A dataset of verbal fluency or linguistic data.
#' @return A binary matrix of responses (rows = participants, cols = responses)
#' @examples
#' responsematrix<-semnetcleaner(data)
#' @references 
#' qdap package: Hornik, K., & Murdoch, D. (2010).
#' Watch Your Spelling!.
#' The R Journal, 3(2), 22-28.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Semantic Network Cleaner----
semnetcleaner<-function(data=data)
{
  #install/load packages
  if (!require("pluralize"))
  {
    install.packages("devtools")
    devtools::install_github("hrbrmstr/pluralize")
  }else library(pluralize)
  
  #perform spell check
  v<-apply(data,c(2),qdap::check_spelling_interactive)
  
  #transform data into a writeable format
  y<-as.data.frame(do.call(cbind,ifelse(v=="NULL",data,v)))
  
  #singularize data
  w<-apply(y,c(2),singularize)
  
  #grab unique responses only
  uni<-rbind(sort(unique(tolower(unlist(apply(w,c(2),unique))))))
  uni[uni==""]<-NA
  uni[uni==" "]<-NA
  uni
  for (i in 1:length(uni))
    if (is.na(uni[i]))
    {
      uni<-uni[-i]
    }

  #attach unique responses to response matrix
  resp<-t(w) #transpose response
  z<-matrix(nrow=nrow(resp),ncol=ncol(uni)) #initialize matrix
  for (i in 1:ncol(resp)) #populate response matrix
  {
    z[,i]<-resp[,i]
  }
  o<-rbind(uni,z) #add unique response to top
  
  #binarize responses
  k<-matrix(nrow=nrow(o),ncol=ncol(o))
  for (i in 1:ncol(o))
    for (j in 2:nrow(o))
      k[j,]<-match(o[1,],o[j,])
  k[is.na(k)]<-0
  for (i in 1:ncol(o))
    for (j in 2:nrow(o))
      if (k[j,i]>0)
      {
        k[j,i]<-1
      }
  k<-k[-1,] #figure out a way so you don't have to do this
  colnames(k)<-o[1,]
  k<-as.data.frame(k)
  
  bin<-k[,-(1:2)] #change as needed
  bin
  write.table(bin,file="filtered.csv",quote=TRUE,row.names=FALSE,col.names=TRUE,sep=",")
}
#----
#' Betwenness Centrality
#' @description Computes betweenness centrlaity of each node in a network (Weighted not coded).
#' @param A An adjacency matrix of network data.
#' @param weighted Is the network weighted? Defaults to TRUE. Set to FALSE for unweighted measure of betwenness centrality.
#' @return A vector of betweenness centrality values for each node in the network.
#' @examples
#' weighted_BC<-Betweenness(A)
#' unweighted_BC<-Betweenness(A,weighted=FALSE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' Neuroimage, 52(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Betweenness Centrality----
Betweenness <- function (A=A,weighted=TRUE)
{
  if(!weighted)
  {A<-ifelse(A!=0,1,0)
  n<-ncol(A)
  I<-diag(60)
  d<-1
  NPd<-A
  NSPd<-NPd
  NSP<-NSPd
  diag(NSP)<-1
  L<-NSPd
  diag(L)<-1
  while (!is.na(which(NSPd!=0)[1]))
  {
    d<-d+1
    NPd<-as.matrix(NPd)%*%as.matrix(A)
    NSPd<-NPd*(L==0)
    NSP<-NSP+NSPd
    L<-L+d*(NSPd!=0)
  }
  L[!L]<-Inf
  diag(L)<-0
  NSP[!NSP]<-1
  At<-t(A)
  DP<-matrix(0,nrow=nrow(A),ncol=ncol(A))
  diam<-d-1
  
  for(d in diam:2)
  {
    DPd1<- (as.matrix(((L==d)*(1+DP)/NSP))%*%as.matrix(At))*((L==(d-1))*NSP)
    DP<-DP+DPd1
  }
  BC<-round(as.data.frame(colSums(DP)),0)
  colnames(BC)<-c("BCu")
  BC}else{print("Weighted not coded--use qgraph")}
}
#----
#' Closeness Centrality
#' @description Computes closeness centrlaity of each node in a network (weighted not coded).
#' @param A An adjacency matrix of network data.
#' @param weighted Is the network weighted? Defaults to TRUE. Set to FALSE for unweighted measure of closeness centrality.
#' @return A vector of closeness centrality values for each node in the network.
#' @examples
#' weighted_LC<-Closeness(A)
#' unweighted_LC<-Closeness(A,weighted=FALSE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' Neuroimage, 52(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export 
#Closeness Centrality----
Closeness <- function (A=A,weighted=TRUE)
{
  if (!weighted)
  {D<-Distance(A,weighted=FALSE)
  C<-matrix(0,ncol=ncol(D))
  for(i in 1:ncol(D))
  {
    C[i]<-1/sum(D[,i])
  }
  LC<-t(as.data.frame(C)*100)
  rownames(LC)<-colnames(A)
  colnames(LC)<-c("LCu")
  LC<-round(as.data.frame(LC),3)
  LC}else{print("Weighted not coded--use qgraph")}
}
#----
#' Degree
#' @description Computes degree of each node in a network.
#' @param A An adjacency matrix of network data.
#' @return A vector of degree values for each node in the network.
#' @examples
#' Deg<-Degree(A)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' Neuroimage, 52(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Degree----
Degree <- function (A=A)
{
  A<-ifelse(A[]!=0,1,0)
  Deg<-as.data.frame(colSums(A))
  colnames(Deg)<-c("Degree")
  Deg
}
#----
#' Node Strength
#' @description Computes strength of each node in a network.
#' @param A An adjacency matrix of network data.
#' @return A vector of strength values for each node in the network.
#' @examples
#' Str<-Strength(A)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' Neuroimage, 52(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Node Strength----
Strength <- function (A=A)
{
  strength<-round(as.data.frame(colSums(A)),2)
  colnames(strength)<-c("Strength")
  strength
}
#----
#' Eigenvector Centrality
#' @description Computes eigenvector centrality of each node in a network.
#' @param A An adjacency matrix of network data.
#' @param weighted Is the network weighted? Defaults to TRUE. Set to FALSE for unweighted measure of eigenvector centrality.
#' @return A vector of eigenvector centrality values for each node in the network.
#' @examples
#' weighted_EC<-Eigenvector(A)
#' unweighted_EC<-Eigenvector(A,weighted=FALSE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' Neuroimage, 52(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Eigenvector

Eigenvector <- function (A=A,weighted=TRUE)
{
  if (!weighted){A<-ifelse(A[]!=0,1,0)
  eigenvector<-eigen(A)
  eigenvector<-round(as.data.frame(abs(eigenvector$vectors[,1])),3)
  colnames(eigenvector)<-c("ECu")}else{eigenvector<-eigen(A)
  eigenvector<-round(as.data.frame(abs(eigenvector$vectors[,1])),3)
  colnames(eigenvector)<-c("ECw")}
  rownames(eigenvector)<-colnames(A)
  eigenvector
}
#----
#' Hybrid Centrality
#' @description Computes hybrid centrality of each node in a network (Weigted not coded).
#' @param A An adjacency matrix of network data.
#' @return A vector of hybrid centrality values for each node in the network (lower values are more central, higher values are more peripheral).
#' @examples
#' HC<-Hybrid(A)
#' @references 
#' Pozzi, F., Di Matteo, T., & Aste, T. (2013).
#' Spread of risk across financial markets: Better to invest in the peripheries. 
#' Scientific Reports, 3(1655), 1-7.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Hybrid Centality----
Hybrid <- function (A=A)
{
  BCu<-Betweenness(A,weighted=FALSE)
  BCw<-Betweenness(A)
  CCu<-Closeness(A,weighted=FALSE)
  CCw<-Closeness(A)
  Deg<-Degree(A)
  Str<-Strength(A)
  ECu<-Eigenvector(A,weighted=FALSE)
  ECw<-Eigenvector(A)
  Eu<-PathLengths(A,weighted=FALSE)$ecc
  Ew<-PathLengths(A)$ecc
  
  hybrid<-((rank(BCu,ties.method="max")+
            rank(BCw,ties.method="max")+
            rank(CCu,ties.method="max")+
            rank(CCw,ties.method="max")+
            rank(Deg,ties.method="max")+
            rank(Str,ties.method="max")+
            rank(ECu,ties.method="max")+
            rank(ECw,ties.method="max")+
            rev(rank(Eu,ties.method="max"))+
            rev(rank(Ew,ties.method="max"))-
            10)/(10*((ncol(A))-10)))
  hybrid<-round(as.data.frame(hybrid),3)
  colnames(hybrid)<-c("HC")
  rownames(hybrid)<-colnames(A)
  hybrid
}
#----
#' List of Centrality Measures
#' @description Computes centrality measures of the network (Weighted not coded).
#' @param A An adjacency matrix of network data.
#' @param weighted Is the network weighted? Defaults to TRUE. Set to FALSE for unweighted list of centrality measures.
#' @return Returns a list of betweenness, closeness, degree (weighted = strength), and eigenvector centralities.
#' @examples
#' weighted_centralitylist<-CentList(A)
#' unweighted_centralitylist<-CentList(A,weighted=FALSE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' Neuroimage, 52(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
CentList <- function (A=A, weighted=TRUE)
{
  if(!weighted){BCu<-Betweenness(A,weighted=FALSE)
  CCu<-Closeness(A,weighted=FALSE)
  Deg<-Degree(A)
  ECu<-Eigenvector(A,weighted=FALSE)}else{print("Weighted not coded.")}
  
  list(Betweenness=BCu,Closeness=CCu,Degree=Deg,Eigenvector=ECu)
}
#' Distance
#' @description Computes distance matrix of the network (Weighted not coded).
#' @param A An adjacency matrix of network data.
#' @param weighted Is the network weighted? Defaults to FALSE. Set to TRUE for weighted measure of distance.
#' @return A distance matrix of the network.
#' @examples
#' unweighted_D<-Distance(A)
#' weighted_D<-Distance(A,weighted=TRUE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' Neuroimage, 52(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Distance:
Distance<-function (A=A,weighted=FALSE)
{
  if(!weighted)
  {B<-ifelse(A[]!=0,1,0)
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
  D}else{print("Weighted not coded.")}
}
#----
#' Path Lengths
#' @description Computes global average shortest path length (ASPL), local average shortest path length (ASPLi), eccentricity (E), and diameter (D) of a network (Weighted not coded).
#' @param A An adjacency matrix of network data.
#' @param weighted Is the network weighted? Defaults to FALSE. Set to TRUE for weighted measures of ASPL, ASPLi, E, and D.
#' @return A list of ASPL, ASPLi, E, and D of a network.
#' @examples
#' unweighted_PL<-PathLengths(A)
#' weighted_PL<-PathLengths(A,weighted=TRUE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' Neuroimage, 52(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Path Lengths
PathLengths <- function (A=A, weighted=FALSE)
{
  if(!weighted)
  {D<-Distance(A,weighted=FALSE)
  n<-nrow(D)
  aspli<-sum(D*(D!=Inf))/(length(which((D!=Inf)!=0)))
  aspl<-sum(sum(D*(D!=Inf))/(length(which((D!=Inf)!=0))))
  Emat<-(D*(D!=Inf))
  ecc<-matrix(nrow=nrow(Emat),ncol=1)
  for(i in 1:nrow(Emat))
  {
    ecc[i,]<-max(Emat[i,])
  }
  d<-max(ecc)
  
  ecc<-as.data.frame(ecc)
  rownames(ecc)<-colnames(A)
  
  list(ASPL=aspl,ASPLi=aspli,Eccentricity=ecc,Diameter=d)}
  else{print("Weighted not coded.")}
}
#----
#' Clustering Coefficient
#' @description Computes global clustering coefficient (CC) and local clustering coefficient (CCi).
#' @param A An adjacency matrix of network data.
#' @param weighted Is the network weighted? Defaults to FALSE. Set to TRUE for weighted measures of CC and CCi.
#' @return A list of CC and CCi.
#' @examples
#' unweighted_CC<-ClustCoeff(A)
#' weighted_CC<-ClustCoeff(A,weighted=TRUE)
#' @references 
#' Rubinov, M., & Sporns, O. (2010). 
#' Complex network measures of brain connectivity: Uses and interpretations. 
#' Neuroimage, 52(3), 1059-1069.
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' @export
#Clustering Coefficient
ClustCoeff <- function (A=A, weighted=FALSE)
{
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
  CC<-mean(C)
  }else{K<-colSums(A!=0)
  m<-A^(1/3)
  cyc<-diag(m%*%m%*%m)
  K[cyc==0]<-Inf
  C<-round(cyc/(K*(K-1)),3)
  CC<-mean(C)}
  list(CC=CC, CCi=C)
}