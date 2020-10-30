#' Critical Correlation Value
#' 
#' @description Computes critical correlation value given sample size and alpha
#' 
#' @param n Numeric.
#' Sample size
#' 
#' @param alpha Numeric.
#' Significance value.
#' Defaults to \code{.05}
#' 
#' @return Critical correlation value
#' 
#' @examples 
#' critical.r(100, .05)
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @noRd
#' 
# Critical Correlation Value----
# Updated 30.10.2020
critical.r <- function(n, alpha)
{
  df <- n - 2
  critical.t <- qt( alpha/2, df, lower.tail = F )
  cvr <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
  return(cvr)
}

#' cpmIV separate
#' 
#' @noRd
#' 
# cpmIV separate----
# Updated 10.09.2020
cpmIV.separate <- function(neuralarray, bstat, kfolds, covar, thresh = .01,
                           groups = NULL, method = c("mean", "sum"),
                           model = c("linear","quadratic","cubic"),
                           corr = c("pearson","spearman"), nEdges, 
                           standardize = FALSE, cores, progBar = TRUE,
                           plots = TRUE)
{
  bstat<-as.vector(bstat)
  if(standardize)
  {bstat<-scale(bstat)}
  
  #number of subjects
  no_sub<-dim(neuralarray)[3]
  #number of nodes
  no_node<-ncol(neuralarray)
  
  if(missing(kfolds)){
    kfolds <- no_sub
  }
  
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
  
  
  #k-folds
  ##kfold breaks
  folds <- cut(1:no_sub, breaks = kfolds, labels = FALSE)
  ##permutate
  perm.folds <- sample(folds, no_sub)
  
  ##reset no_sub as list
  no_sub <- vector("list", length = kfolds)
  
  for(i in 1:kfolds){
    no_sub[[i]] <- which(perm.folds == i)
  }
  
  #perform leave-out analysis
  if(progBar)
  {pb <- txtProgressBar(max=length(no_sub), style = 3)}
  
  for(k in 1:length(no_sub))
  {
    #leftout
    leftout <- no_sub[[k]]
    
    #get neuralarray
    train_mats<-neuralarray
    train_mats<-train_mats[,,-leftout]
    
    #get parts
    no_train <- dim(train_mats)[3]
    
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
    if(nrow(train_vcts)!=(no_train))
    {train_vcts<-t(train_vcts)}
    
    rmat<-vector(mode="numeric",length=ncol(train_vcts))
    pmat<-vector(mode="numeric",length=ncol(train_vcts))
    
    if(is.list(covar))
    {
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      
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
      cvr<-critical.r((no_train),thresh)
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
    train_sumpos<-matrix(0,nrow=(no_train),ncol=1)
    train_sumneg<-matrix(0,nrow=(no_train),ncol=1)
    
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
    
    if(any(is.na(fit_pos)))
    {fit_pos[which(is.na(fit_pos))] <- 0}
    
    if(any(is.na(fit_neg)))
    {fit_neg[which(is.na(fit_neg))] <- 0}
    
    #run model on TEST sub
    test_mat<-neuralarray[,,leftout]
    
    if(is.na(dim(test_mat)[3])){
      
      if(method=="sum"){
        
        test_sumpos<-sum(test_mat*pos_mask)/2
        test_sumneg<-sum(test_mat*neg_mask)/2
        
      }else if(method=="mean"){
        
        test_sumpos<-mean(test_mat*pos_mask)/2
        test_sumneg<-mean(test_mat*neg_mask)/2
        
      }
      
    }else{
      
      if(method=="sum"){
        
        test_sumpos <- apply(test_mat, 3, function(x){
          sum(x * pos_mask) / 2
        })
        
        test_sumneg <- apply(test_mat, 3, function(x){
          sum(x * neg_mask) / 2
        })
        
      }else if(method=="mean"){
        
        test_sumpos <- apply(test_mat, 3, function(x){
          mean(x * pos_mask) / 2
        })
        
        test_sumneg <- apply(test_mat, 3, function(x){
          mean(x * neg_mask) / 2
        })
      }
      
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
    {setTxtProgressBar(pb, k)}
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
  
  posmask <- ifelse(pos_mat>=nEdges,1,0)
  negmask <- ifelse(neg_mat>=nEdges,1,0)
  
  if(!is.null(colnames(neuralarray)))
  {
    colnames(posmask) <- colnames(neuralarray)
    row.names(posmask) <- colnames(posmask)
    colnames(negmask) <- colnames(neuralarray)
    row.names(negmask) <- colnames(negmask)
  }
  
  R_pos<-cor(behav_pred_pos,bstat,use="pairwise.complete.obs")
  P_pos<-cor.test(behav_pred_pos,bstat)$p.value
  R_neg<-cor(behav_pred_neg,bstat,use="pairwise.complete.obs")
  P_neg<-cor.test(behav_pred_neg,bstat)$p.value
  
  P_pos<-ifelse(round(P_pos,3)!=0,round(P_pos,3),noquote("< .001"))
  P_neg<-ifelse(round(P_neg,3)!=0,round(P_neg,3),noquote("< .001"))
  
  bstat<-as.vector(bstat)
  behav_pred_pos<-as.vector(behav_pred_pos)
  behav_pred_neg<-as.vector(behav_pred_neg)
  
  #error
  perror <- behav_pred_pos - bstat
  nerror <- behav_pred_neg - bstat
  
  #mae
  mae_pos <- mean(abs(perror), na.rm = TRUE)
  mae_neg <- mean(abs(nerror), na.rm = TRUE)
  
  #rmse
  pos_rmse <- sqrt(mean(perror^2, na.rm = TRUE))
  neg_rmse <- sqrt(mean(nerror^2, na.rm = TRUE))
  
  results<-matrix(0,nrow=2,ncol=4)
  
  results[1,1]<-round(R_pos,3)
  results[1,2]<-P_pos
  results[1,3]<-round(mae_pos,3)
  results[1,4]<-round(pos_rmse,3)
  results[2,1]<-round(R_neg,3)
  results[2,2]<-P_neg
  results[2,3]<-round(mae_neg,3)
  results[2,4]<-round(neg_rmse,3)
  
  colnames(results)<-c("r","p","mae","rmse")
  row.names(results)<-c("positive","negative")
  
  #Results list
  res <- list()
  res$results <- results
  res$posMask <- posmask
  res$negMask <- negmask
  res$posArray <- neuralarray
  res$negArray <- neuralarray
  
  for(i in 1:dim(neuralarray)[3])
  {
    res$posArray[,,i] <- posmask * neuralarray[,,i]
    res$negArray[,,i] <- negmask * neuralarray[,,i]
  }
  
  res$behav <- bstat
  res$posPred <- behav_pred_pos
  res$negPred <- behav_pred_neg
  res$groups <- groups
  res$connections <- "separate"
  
  class(res) <- "cpm"
  
  if(plots)
  {plot(res)}
  
  return(res)
}

#' cpmIV overall
#' 
#' @noRd
#' 
# cpmIV overall----
# Updated 10.09.2020
cpmIV.overall <- function(neuralarray, bstat, kfolds, covar, thresh = .01,
                           groups = NULL, method = c("mean", "sum"),
                           model = c("linear","quadratic","cubic"),
                           corr = c("pearson","spearman"), nEdges, 
                           standardize = FALSE, cores, progBar = TRUE,
                           plots = TRUE)
{
  bstat<-as.vector(bstat)
  if(standardize)
  {bstat<-scale(bstat)}
  
  #number of subjects
  no_sub<-dim(neuralarray)[3]
  #number of nodes
  no_node<-ncol(neuralarray)
  
  if(missing(kfolds)){
    kfolds <- no_sub
  }
  
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
  
  #k-folds
  ##kfold breaks
  folds <- cut(1:no_sub, breaks = kfolds, labels = FALSE)
  ##permutate
  perm.folds <- sample(folds, no_sub)
  
  ##reset no_sub as list
  no_sub <- vector("list", length = kfolds)
  
  for(i in 1:kfolds){
    no_sub[[i]] <- which(perm.folds == i)
  }
  
  #perform leave-out analysis
  if(progBar)
  {pb <- txtProgressBar(max=length(no_sub), style = 3)}
  
  for(k in 1:length(no_sub))
  {
    #leftout
    leftout <- no_sub[[k]]
    
    #get neuralarray
    train_mats<-neuralarray
    train_mats<-train_mats[,,-leftout]
    
    #get parts
    no_train <- dim(train_mats)[3]
    
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
    if(nrow(train_vcts)!=(no_train))
    {train_vcts<-t(train_vcts)}
    
    rmat<-vector(mode="numeric",length=ncol(train_vcts))
    pmat<-vector(mode="numeric",length=ncol(train_vcts))
    
    if(is.list(covar))
    {
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      
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
      cvr<-critical.r((no_train),thresh)
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
    train_sumpos<-matrix(0,nrow=(no_train),ncol=1)
    train_sumneg<-matrix(0,nrow=(no_train),ncol=1)
    
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
    {fit_pos<-coef(lm(train_behav~train_sumpos+train_sumneg))
    }else if(model=="quadratic")
    {
      quad_pos<-train_sumpos^2
      quad_neg<-train_sumneg^2
      
      fit_pos<-coef(lm(train_behav~train_sumpos+quad_pos+train_sumneg+quad_neg))
      
    }else if(model=="cubic")
    {
      cube_pos<-train_sumpos^3
      cube_pos<-train_sumneg^3
      
      quad_pos<-train_sumpos^2
      quad_pos<-train_sumneg^2
      
      fit_pos<-coef(lm(train_behav~train_sumpos+quad_pos+cube_pos+train_sumneg+quad_neg+cube_neg))
    }
    
    if(any(is.na(fit_pos)))
    {fit_pos[which(is.na(fit_pos))] <- 0}
    
    #run model on TEST sub
    test_mat<-neuralarray[,,leftout]
    
    if(is.na(dim(test_mat)[3])){
      
      if(method=="sum"){
        
        test_sumpos<-sum(test_mat*pos_mask)/2
        test_sumneg<-sum(test_mat*neg_mask)/2
        
      }else if(method=="mean"){
        
        test_sumpos<-mean(test_mat*pos_mask)/2
        test_sumneg<-mean(test_mat*neg_mask)/2
        
      }
      
    }else{
      
      if(method=="sum"){
        
        test_sumpos <- apply(test_mat, 3, function(x){
          sum(x * pos_mask) / 2
        })
        
        test_sumneg <- apply(test_mat, 3, function(x){
          sum(x * neg_mask) / 2
        })
        
      }else if(method=="mean"){
        
        test_sumpos <- apply(test_mat, 3, function(x){
          mean(x * pos_mask) / 2
        })
        
        test_sumneg <- apply(test_mat, 3, function(x){
          mean(x * neg_mask) / 2
        })
      }
      
    }
    
    if(model=="linear")
    {behav_pred_pos[leftout]<-fit_pos[3]*test_sumneg+fit_pos[2]*test_sumpos+fit_pos[1]
    }else if(model=="quadratic")
    {
      quad_post<-test_sumpos^2
      quad_negt<-test_sumneg^2
      
      behav_pred_pos[leftout]<-fit_pos[5]*quad_negt+fit_pos[4]*test_sumneg+fit_pos[3]*quad_post+fit_pos[2]*test_sumpos+fit_pos[1]
    }else if(model=="cubic")
    {
      cube_post<-test_sumpos^3
      cube_negt<-test_sumneg^3
      
      quad_negt<-test_sumpos^2
      quad_negt<-test_sumneg^2
      
      behav_pred_pos[leftout]<-fit_pos[7]*cube_negt+fit_pos[6]*quad_negt+fit_pos[5]*test_sumneg+fit_pos[4]*cube_post+fit_pos[3]*quad_post+fit_pos[2]*test_sumpos+fit_pos[1]
    }
    
    if(progBar)
    {setTxtProgressBar(pb, k)}
  }
  if(progBar)
  {close(pb)}
  
  pos_mat <- matrix(0, nrow = no_node, ncol = no_node)
  
  for(i in 1:no_node)
    for(j in 1:no_node)
    {pos_mat[i,j] <- sum(pos_array[i,j,])}
  
  posmask <- ifelse(pos_mat>=nEdges,1,0)
  
  if(!is.null(colnames(neuralarray)))
  {
    colnames(posmask) <- colnames(neuralarray)
    row.names(posmask) <- colnames(posmask)
  }
  
  R_pos<-cor(behav_pred_pos,bstat,use="pairwise.complete.obs")
  P_pos<-cor.test(behav_pred_pos,bstat)$p.value
 
  P_pos<-ifelse(round(P_pos,3)!=0,round(P_pos,3),noquote("< .001"))
  
  bstat<-as.vector(bstat)
  behav_pred_pos<-as.vector(behav_pred_pos)
  
  #error
  perror <- behav_pred_pos - bstat
  
  #mae
  mae_pos <- mean(abs(perror), na.rm = TRUE)
  
  #rmse
  pos_rmse <- sqrt(mean(perror^2, na.rm = TRUE))
  
  results<-matrix(0,nrow=1,ncol=4)
  
  results[1,1]<-round(R_pos,3)
  results[1,2]<-P_pos
  results[1,3]<-round(mae_pos,3)
  results[1,4]<-round(pos_rmse,3)
  
  colnames(results)<-c("r","p","mae","rmse")
  row.names(results)<-c("overall")
  
  #Results list
  res <- list()
  res$results <- results
  res$Mask <- posmask
  res$Array <- neuralarray
  
  for(i in 1:dim(neuralarray)[3])
  {res$Array[,,i] <- posmask * neuralarray[,,i]}
  
  res$behav <- bstat
  res$Pred <- behav_pred_pos
  res$groups <- groups
  res$connections <- "overall"
  
  class(res) <- "cpm"
  
  if(plots)
  {plot(res)}
  
  return(res)
}
