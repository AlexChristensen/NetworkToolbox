#Network Toolbox Functions
#' @export
#TMFG Filtering Method----

TMFG <-function (data=data)
{
  cormat<-cor(data)
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
  for(r in 1:nrow(x))
    for(z in 1:ncol(x))
    {
      if(x[r,z]==1)
      {
        x[r,z]<-cormat[r,z]
      }
    }
  x<-as.data.frame(x)
  colnames(x)<-colnames(cormat)
  x
}
#----
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
  #read in file
  data<-read.csv("trial.csv",header=FALSE,sep=",",as.is=TRUE)
  
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
  uni
  
  badToken<-c(""," ")
  for (i in badToken)
  {
    uni[uni==i]<-NA
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
#' @export
#Betweenness Centrality: Binary----

Betweenness <- function (A=A)
{
  B<-ifelse(A[]!=0,1,0)
  n<-ncol(B)
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
  
  BC<-round(as.data.frame(colSums(DP)),0)
  colnames(BC)<-c("BC")
  BC
}
#----
#' @export 
#Closeness Centrality: Binary----

Closeness <- function (A=A)
{
  D<-distance_bin(A)
  C<-matrix(0,ncol=ncol(D))
  for(i in 1:ncol(D))
  {
    C[i]<-1/sum(D[,i])
  }
  CC<-t(as.data.frame(C)*100)
  rownames(CC)<-colnames(A)
  colnames(CC)<-c("CC")
  CC<-round(as.data.frame(CC),3)
  CC
}
#----
#' @export
#Degree----

Degree <- function (A=A)
{
  B<-ifelse(A[]!=0,1,0)
  Deg<-as.data.frame(colSums(B))
  colnames(Deg)<-c("Degree")
  Deg
}
#----
#' @export
#Node Strength----

Strength <- function (A=A)
{
  strength<-round(as.data.frame(colSums(A)),2)
  colnames(strength)<-c("Strength")
  strength
}
#----
#' @export
#Eigenvector: Binary----

EigenU <- function (A=A)
{
  B<-ifelse(A[]!=0,1,0)
  eigenvectoru<-eigen(B)
  eigenvectoru<-round(as.data.frame(abs(eigenvectoru$vectors[,1])),3)
  colnames(eigenvectoru)<-c("EigenU")
  rownames(eigenvectoru)<-colnames(A)
  eigenvectoru
}
#----
#' @export
#Eigenvector: Weighted----

EigenW <- function (A=A)
{
  eigenvectorw<-eigen(A)
  eigenvectorw<-round(as.data.frame(abs(eigenvectorw$vectors[,1])),3)
  colnames(eigenvectorw)<-c("EigenW")
  rownames(eigenvectorw)<-colnames(A)
  eigenvectorw
}
#----
#' @export
#Hybrid Centality----

Hybrid <- function (A=A)
{
  BCu<-Betweenness(A)
  CCu<-Closeness(A)
  Deg<-Degree(A)
  Str<-Strength(A)
  ECu<-EigenU(A)
  ECw<-EigenW(A)
  
  hybrid<-((rank(BCu,ties.method="max")+
              rank(CCu,ties.method="max")+
              rank(Deg,ties.method="max")+
              rank(Str,ties.method="max")+
              rank(ECw,ties.method="max")+
              rank(ECu,ties.method="max")+
              -6)/(6*((ncol(A))-6)))
  hybrid<-round(as.data.frame(hybrid),3)
  colnames(hybrid)<-c("HC")
  rownames(hybrid)<-colnames(A)
  hybrid
}
#----
#' @export
#Centrality Table----

CentTable <- function (A)
{
  BCu<-Betweenness(A)
  CCu<-Closeness(A)
  Deg<-Degree(A)
  Str<-Strength(A)
  ECu<-EigenU(A)
  ECw<-EigenW(A)
  HC<-Hybrid(A)
  cent<-cbind(BCu,CCu,Deg,Str,ECu,ECw,HC)
  list(CentrailtyTable=cent)
}
#----
#' @export
#List of Centralities----

AllCents <- function (A)
{
  BCu<-Betweenness(A)
  CCu<-Closeness(A)
  Deg<-Degree(A)
  Str<-Strength(A)
  ECu<-EigenU(A)
  ECw<-EigenW(A)
  HC<-Hybrid(A)
  list(BCu=BCu,CCu=CCu,Deg=Deg,Str=Str,ECu=ECu,ECw=ECw,HC=HC)
}
#----
#' @export
#Distnace: Binary----

distance_bin<-function (A=A)
{
  B<-ifelse(A[]!=0,1,0)
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
  D
}
#----
#' @export
#Average Shortest Path Length: Binary----

PathLengths <- function (A=A)
{
  n<-nrow(D)
  D<-distance_bin(A)
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
  
  list(ASPL=aspl,Eccentricity=ecc,Diameter=d)
}
