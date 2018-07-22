#' Connectome-based Predictive Modeling--External Validation
#' @description Applies the Connectome-based Predictive Modeling approach to neural data.
#' This method predicts a behavioral statistic using neural connectivity from the sample.
#' Results \strong{may} differ from \emph{Matlab} results because of robust GLM methodology.
#' This function is still in its \strong{testing phase}.
#' \strong{Please cite Finn et al., 2015; Rosenberg et al., 2016; Shen et al., 2017}
#' 
#' @param train_na Training dataset
#' (an array from \code{\link{convertConnBrainMat}} function)
#' 
#' @param train_b Behavioral statistic for each participant for the \strong{training} neural data (a vector)
#' 
#' @param valid_na Validation dataset
#' (an array from \code{\link{convertConnBrainMat}} function)
#' 
#' @param valid_b Behavioral statistic for each participant for the \strong{validation} neural data (a vector)
#' 
#' @param thresh Sets an \strong{alpha} threshold for edge weights to be retained.
#' Defaults to .01
#' 
#' @param overlap Should leave-one-out cross-validation be used?
#' Defaults to FALSE (use full dataset, no leave-one-out).
#' Set to TRUE to select edges that appear in every leave-one-out cross-validation network (\emph{time consuming})
#' 
#' @param progBar Should progress bar be displayed?
#' Defaults to TRUE.
#' Set to FALSE for no progress bar
#' 
#' @return Returns a list containing: 
#'
#' \item{results}{A matrix contaning: r coefficient (r), p-value (p-value),
#' mean absolute error (mae), root mean square error (rmse)}
#' 
#' \item{posMask}{Positive connectivity for input in
#' \href{https://bioimagesuiteweb.github.io/webapp/connviewer.html}{BioImage Suite Connectivity Viewer}}
#' 
#' \item{negMask}{Negative connectivity for input in
#' \href{https://bioimagesuiteweb.github.io/webapp/connviewer.html}{BioImage Suite Connectivity Viewer}}
#' 
#' @references 
#' Finn, E. S., Shen, X., Scheinost, D., Rosenberg, M. D., Huang, J., Chun, M. M., Papademetris, X., Constable, R. T. (2015).
#' Functional connectome fingerprinting: Identifying individuals using patterns of brain connectivity.
#' \emph{Nature Neuroscience}, \emph{18}, 1664-1671.
#' doi: \href{https://doi.org/10.1038/nn.4135}{10.1038/nn.4135}
#' 
#' Rosenberg, M. D., Finn, E. S., Scheinost, D., Papademetris, X., Shen, X., Constable, R. T., Chun, M. M. (2016).
#' A neuromarker of sustained attention from whole-brain functional connectivity.
#' \emph{Nature Neuroscience}, \emph{19}, 165-171.
#' doi: \href{https://doi.org/10.1038/nn.4179}{10.1038/nn.4179}
#'
#' Shen, X. Finn, E. S., Scheinost, D., Rosenberg, M. D., Chun, M. M., Papademetris, X., Constable, R. T. (2017).
#' Using connectome-based predictive modeling to predict individual behavior from brain connectivity.
#' \emph{Nature Protocols}, \emph{12}, 506-518.
#' doi: \href{https://doi.org/10.1038/nprot.2016.178}{10.1038/nprot.2016.178}
#' 
#' Wei, T. & Simko, V.(2017).
#' R package "corrplot": Visualization of a correlation matrix (Version 0.84).
#' Available from \url{https://github.com/taiyun/corrplot}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @importFrom MASS psi.bisquare
#' 
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
    
    perror <- vector(mode="numeric",length = length(valid_b))
    nerror <- vector(mode="numeric",length = length(valid_b))
    glmerr <- vector(mode="numeric",length = length(valid_b))
    
    for(i in 1:length(valid_b))
    {
        perror <- pred_pos[i]-valid_b[i]
        nerror <- pred_neg[i]-valid_b[i]
        glmerr <- pred_glm[i]-valid_b[i]
        
        #mae
        mae_pos<-mean(abs(perror))
        mae_neg<-mean(abs(nerror))
        mae_glm<-mean(abs(glmerr))
        
        #rmse
        rmse_pos<-sqrt(mean(perror^2))
        rmse_neg<-sqrt(mean(nerror^2))
        rmse_glm<-sqrt(mean(glmerr^2))
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
    
    results<-matrix(0,nrow=3,ncol=4)
    
    results[1,1]<-round(r_pos,3)
    results[1,2]<-p_pos
    results[1,3]<-round(mae_pos,3)
    results[1,4]<-round(rmse_pos,3)
    results[2,1]<-round(r_neg,3)
    results[2,2]<-p_neg
    results[2,3]<-round(mae_neg,3)
    results[2,4]<-round(rmse_neg,3)
    results[3,1]<-round(r_glm,3)
    results[3,2]<-p_glm
    results[3,3]<-round(mae_glm,3)
    results[3,4]<-round(rmse_pos,3)
    
    colnames(results)<-c("r","p-value","mae","rmse")
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