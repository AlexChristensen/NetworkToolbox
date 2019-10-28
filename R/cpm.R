#' Connectome Predictive Modeling
#' 
#' @name cpm
#'
#' @aliases
#' cpmEV
#' cpmFP
#' cpmFPperm
#' cpmIV
#' cpmPlot
#' 
#' @description Suite of functions for Connectome Predictive Modeling (CPM).
#' \strong{See and cite Finn et al., 2015; Rosenberg et al., 2016; Shen et al., 2017}
#' 
#' \itemize{
#' 
#' \item{\code{cpmIV}}
#'
#' {Internal Validation method (Rosenberg et al., 2016; Shen et al., 2017). Using a leave-one-out approach,
#' this method correlates a behavioral statistic \code{bstat} with each edge of a whole-brain network across
#' participants. Using the significant edges in the network \code{thresh}, a connectome model
#' is built (without the participant's network). A linear regression model is fit, with the behavioral
#' statistic being regressed on the connectome model. The left out participants connectome model is then
#' used with the linear regression weights to compute their predicted behavioral score. This is repeated
#' for every participant. The predicted scores are correlated with their observed score. Significant values
#' suggest that the connectome is related to the behavioral statistic}
#' 
#' \item{\code{cpmEV}}
#' 
#' {UNDER DEVELOPMENT. External Validation method (Beaty et al., 2018). Performs similar function as \code{cpmIV} but uses data
#' to train \code{train_na} the connectome model using a behavioral statistic \code{train_b}.
#' This training connectome model is then used to predict another dataset \code{valid_na},
#' using the same behavioral statistic \code{valid_b}. The full training dataset \code{FALSE} or
#' the leave-one-out \code{overlap = TRUE} approach can be used}
#' 
#' \item{\code{cpmFP}}
#' 
#' {Fingerprinting method (Finn et al., 2015). Uses CPM approach to identify participants across two
#' sessions}
#' 
#' \item{\code{cpmFPperm}}
#' 
#' {Fingerprinting method (Finn et al., 2015). Uses permutation method to estimate the significance of
#' of the \code{cpmFP} results}
#' 
#' \item{\code{cpmPlot}}
#' 
#' {Plots the CPM results}
#' 
#' }
#' 
#' @usage
#' cpmIV(neuralarray, bstat, covar, thresh = .01, method = c("mean", "sum"),
#'       model = c("linear","quadratic","cubic"), corr = c("pearson","spearman"),
#'       nEdges, cores, progBar = TRUE)
#'       
#' cpmEV(train_na, train_b, valid_na, valid_b, thresh = .01,
#'       overlap = FALSE, progBar = TRUE)
#'       
#' cpmFP(session1, session2, progBar = TRUE)
#' 
#' cpmFPperm(session1, session2, iter = 1000, progBar = TRUE)
#' 
#' cpmPlot(cpm.obj, visual.nets = FALSE)
#' 
#' @param neuralarray Array from \code{\link[NetworkToolbox]{convertConnBrainMat}} function
#' 
#' @param bstat Behavioral statistic for each participant with neural data (a vector)
#' 
#' @param covar Covariates to be included in predicting relevant edges (\strong{time consuming}).
#' \strong{Must} be input as a \code{list()} (see examples)
#' 
#' @param thresh Sets an \eqn{\alpha} threshold for edge weights to be retained.
#' Defaults to \code{.01}
#' 
#' @param method Use \code{"mean"} or \code{"sum"} of edge strengths in the positive and negative connectomes.
#' Defaults to \code{"mean"}
#' 
#' @param model Regression model to use for fitting the data.
#' Defaults to \code{"linear"}
#' 
#' @param corr Correlation method for assessing the relationship between the behavioral measure and edges between ROIs.
#' Defaults to \code{"pearson"}.
#' Set to \code{"spearman"} for non-linear or monotonic associations
#' 
#' @param nEdges Number of participants that are required
#' to have an edge to appear in the plots.
#' Defaults to 10 percent of edges in participants
#' 
#' @param cores Number of computer processing cores to use when performing covariate analyses.
#' Defaults to \emph{n} - 1 total number of cores.
#' Set to any number between 1 and maximum amount of cores on your computer
#' 
#' @param progBar Should progress bar be displayed?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} for no progress bar
#' 
#' @param train_na Training dataset
#' (an array from \code{\link[NetworkToolbox]{convertConnBrainMat}} function)
#' 
#' @param train_b Behavioral statistic for each participant for the \strong{training} neural data (a vector)
#' 
#' @param valid_na Validation dataset
#' (an array from \code{\link[NetworkToolbox]{convertConnBrainMat}} function)
#' 
#' @param valid_b Behavioral statistic for each participant for the \strong{validation} neural data (a vector)
#' 
#' @param overlap Should leave-one-out cross-validation be used?
#' Defaults to \code{FALSE} (use full dataset, no leave-one-out).
#' Set to \code{TRUE} to select edges that appear in every leave-one-out cross-validation network (\emph{time consuming})
#' 
#' @param session1 Array from \code{\link[NetworkToolbox]{convertConnBrainMat}} function
#' (first session)
#' 
#' @param session2 Array from \code{\link[NetworkToolbox]{convertConnBrainMat}} function
#' (second session)
#' 
#' @param iter Number of iterations to perform.
#' Defaults to \code{1000}
#' 
#' @param cpm.obj \code{\link[NetworkToolbox]{cpm}} object
#' 
#' @param visual.nets Boolean.
#' Uses \code{\link[qgraph]{qgraph}} to plot connectivity
#' between the networks as a network.
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} to visualize the networks
#' 
#' @return
#' 
#' \code{cpmIV} and \code{cpmEV}:
#' 
#' Returns a list containing:
#'
#' \item{results}{A matrix containing: r coefficient (\code{r}), p-value (\code{p-value}),
#' mean absolute error (\code{mae}), root mean square error (\code{rmse})}
#' 
#' \item{posMask}{Positive connectivity for input in
#' \href{https://bioimagesuiteweb.github.io/webapp/connviewer.html}{BioImage Suite Connectivity Viewer}}
#' 
#' \item{negMask}{Negative connectivity for input in
#' \href{https://bioimagesuiteweb.github.io/webapp/connviewer.html}{BioImage Suite Connectivity Viewer}}
#' 
#' \code{cpmFP}:
#' 
#' Returns a matrix containing the percentage
#' and number of correctly identified subjects for sessions 1 and 2
#' 
#' \code{cpmPlot}:
#' 
#' Returns plot of connectivity differences between the
#' positive and negative masks
#' 
#' @references
#' Beaty, R. E., Kenett, Y. N., Christensen, A. P., Rosenberg, M. D., Benedek, M., Chen, Q.,
#' Fink, A., Qiu, J., Kwapil, T. R., Kane, M. J., & Silvia, P. J. (2018).
#' Robust prediction of individual creative ability from brain functional connectivity.
#' \emph{Proceedings of the National Academy of Sciences}, \emph{115}, 1087-1092.
#' doi:\href{https://doi.org/10.1073/pnas.1713532115}{10.1073/pnas.1713532115}
#' 
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
#' @examples 
#' # Load data
#' behav <- behavOpen
#' 
#' \dontrun{
#' 
#' # Download associated brain data
#' https://drive.google.com/file/d/1ugwi7nRrlHQYuGPzEB4wYzsizFrIMvKR/view
#' 
#' # Load brain data
#' load("restOpen.rda")
#' 
#' # Run cpmIV
#' res <- cpmIV(neuralarray = restOpen, bstat = behav, cores = 4)
#' 
#' # Plot cpmIV results
#' cpmPlot(res)
#' 
#' }
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @importFrom foreach %dopar%
#' @importFrom stats cor.test coef cov2cor lm na.omit
#' @importFrom graphics par abline hist plot text
#' @importFrom grDevices colorRampPalette dev.new
#' @importFrom MASS psi.bisquare
#' 
#' @export
#CPM Functions----
#CPM Internal Validation----
cpmIV <- function (neuralarray, bstat, covar, thresh = .01,
                   method = c("mean", "sum"),
                   model = c("linear","quadratic","cubic"),
                   corr = c("pearson","spearman"), nEdges, 
                   cores, progBar = TRUE)
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
    
    if(missing(nEdges))
    {nEdges<-length(bstat)*.10
    }else{nEdges <- nEdges}
    
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
    
    posmask <- ifelse(pos_mat>=nEdges,1,0)
    negmask <- ifelse(neg_mat>=nEdges,1,0)
    
    
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
    
    bstat<-as.vector(bstat)
    behav_pred_pos<-as.vector(behav_pred_pos)
    behav_pred_neg<-as.vector(behav_pred_neg)
    perror <- vector(mode="numeric",length = length(bstat))
    nerror <- vector(mode="numeric",length = length(bstat))
    
    for(i in 1:length(bstat))
    {
        perror[i] <- behav_pred_pos[i]-bstat[i]
        nerror[i] <- behav_pred_neg[i]-bstat[i]
        
        #mae
        mae_pos<-mean(abs(perror))
        mae_neg<-mean(abs(nerror))
        
        #rmse
        pos_rmse<-sqrt(mean(perror^2))
        neg_rmse<-sqrt(mean(nerror^2))
    }
    
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
    
    class(res) <- "CPM"
    
    return(res)
}
#----
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
    
    #Results list
    res <- list()
    res$results <- results
    res$posMask <- pos_mask
    res$negMask <- neg_mask
    
    class(res) <- "CPM"
    
    return(res)
}
#----
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
#Plots for CPM Results----
cpmPlot <- function (cpm.obj, visual.nets = FALSE)
{
    # Check if CPM object
    if(class(cpm.obj) != "CPM")
    {stop("'cpm.obj' is not a 'CPM' object")}
    
    # Masks
    posmask <- cpm.obj$posMask
    negmask <- cpm.obj$negMask
    
    # Number of nodes
    n <- ncol(posmask)
    
    # Determine plots
    if(n == 268)
    {atlas <- "Shen"
    }else if(n == 300)
    {atlas <- "Schaefer"}
    
    # Shen plots----
    if(atlas == "Shen")
    {
        ####################
        #### Shen Atlas ####
        ####################
        
        # Macroscale Resgions
        rtlobenets <- c(rep("PFC",22),rep("Mot",11),rep("Ins",4),rep("Par",13),
                        rep("Tem",21),rep("Occ",11),rep("Lim",17),rep("Cer",20),
                        rep("Sub",9),rep("Bsm",5))
        
        ltlobenets <- c(rep("PFC",24),rep("Mot",10),rep("Ins",3),rep("Par",14),
                        rep("Tem",18),rep("Occ",14),rep("Lim",19),rep("Cer",21),
                        rep("Sub",8),rep("Bsm",4))
        
        atlasReg <- c(rtlobenets, ltlobenets)
        
        lobe.title <- "Macroscale Regions"
        
        # Canonical Networks
        atlasNet <- c(2,4,3,2,3,3,2,2,2,1,4,1,3,2,4,1,2,4,2,4,2,
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
        
        atlasNet.names <- c("MF","FP","DM","SubC","MT","VI","VII","VA")
        
        for(i in 1:length(unique(atlasNet.names)))
        {atlasNet[which(atlasNet==i)] <- atlasNet.names[i]}
        
        con.title <- "Canonical Networks"
        
    }else if(atlas == "Schaefer")
    {
        ########################
        #### Schaefer Atlas ####
        ########################
        
        # 17-Network Parcellation
        ltlobenets <- c(rep("VisCent",11),rep("VisPeri",9),rep("SomMotA",15),rep("SomMotB",12),
                        rep("DorsAttnA",8),rep("DorsAttnB",8),rep("SalVentAttnA",10),
                        rep("SalVentAttnB",6),rep("Limbic",9),rep("ContA",12),rep("ContB",6),
                        rep("ContC",4),rep("DefaultA",11),rep("DefaultB",18),rep("DefaultC",6),
                        rep("TempPar",5))
        
        rtlobenets <- c(rep("VisCent",11),rep("VisPeri",9),rep("SomMotA",15),rep("SomMotB",9),
                        rep("DorsAttnA",8),rep("DorsAttnB",8),rep("SalVentAttnA",15),
                        rep("SalVentAttnB",8),rep("Limbic",11),rep("ContA",10),rep("ContB",12),
                        rep("ContC",4),rep("DefaultA",12),rep("DefaultB",7),rep("DefaultC",4),
                        rep("TempPar",7))
        
        atlasReg <- c(ltlobenets,rtlobenets)
        
        lobe.title <- "17-Network Parcellation"
        
        # 7-Network Parcellation
        ltnets <- c(rep("Vis",24),rep("SomMot",29),rep("DorsAttn",16),rep("SalVentAttn",16),
                    rep("Limbic",10),rep("Cont",17),rep("Default",38))
        
        rtnets <- c(rep("Vis",23),rep("SomMot",28),rep("DorsAttn",18),rep("SalVentAttn",18),
                    rep("Limbic",10),rep("Cont",23),rep("Default",30))
        
        atlasNet <- c(ltnets,rtnets)
        
        con.title <- "7-Network Parcellation"
    }
    
    ############################
    #### Macroscale Regions ####
    ############################
    
    # Names of regions
    lobes <- unique(atlasReg)
    
    pos_lobe <- posmask
    neg_lobe <- negmask
    colnames(pos_lobe) <- atlasReg
    colnames(neg_lobe) <- atlasReg
    
    poslobemat <- matrix(0, nrow = length(lobes), ncol = length(lobes))
    neglobemat <- matrix(0, nrow = length(lobes), ncol = length(lobes))
    
    for(i in 1:length(lobes))
        for(j in 1:length(lobes))
        {
            poslobemat[i,j] <- sum(pos_lobe[which(colnames(pos_lobe)==lobes[i]),which(colnames(pos_lobe)==lobes[j])])
            neglobemat[i,j] <- sum(neg_lobe[which(colnames(neg_lobe)==lobes[i]),which(colnames(neg_lobe)==lobes[j])])
        }
    
    ldiffmat <- (poslobemat - neglobemat)
    
    colnames(ldiffmat) <- lobes
    row.names(ldiffmat) <- lobes
    colnames(poslobemat) <- lobes
    row.names(poslobemat) <- lobes
    colnames(neglobemat) <- lobes
    row.names(neglobemat) <- lobes
    
    ldiffmat[upper.tri(ldiffmat)] <- 0
    
    llim <- ifelse(abs(min(ldiffmat))>max(ldiffmat),abs(min(ldiffmat)),max(ldiffmat))
    
    colo <- colorRampPalette(c("skyblue2","white","darkorange2"))
    
    dev.new()
    corrplot::corrplot(ldiffmat,is.corr=FALSE,method="color",
                       tl.col="black",col = colo(100),na.label="square",
                       na.label.col = "white",addgrid.col="black",
                       title=paste("Difference in the Number of Edges\nin",lobe.title,sep=" "),
                       mar=c(0,0,4,0),cl.length=3,cl.pos="b",cl.lim=c(-llim,llim))
    
    
    ############################
    #### Canonical Networks ####
    ############################
    
    # Names of regions
    nets <- unique(atlasNet)
    
    pos_nets <- posmask
    neg_nets <- negmask
    
    colnames(pos_nets) <- atlasNet
    colnames(neg_nets) <- atlasNet
    
    posnetmat <- matrix(0,nrow=length(nets),ncol=length(nets))
    negnetmat <- matrix(0,nrow=length(nets),ncol=length(nets))
    
    for(i in 1:length(nets))
        for(j in 1:length(nets))
        {
            posnetmat[i,j] <- sum(pos_nets[which(colnames(pos_nets)==nets[i]),which(colnames(pos_nets)==nets[j])])
            negnetmat[i,j] <- sum(neg_nets[which(colnames(neg_nets)==nets[i]),which(colnames(neg_nets)==nets[j])])
        }
    
    diffmat <- (posnetmat-negnetmat)
    
    colnames(diffmat) <- nets
    row.names(diffmat) <- nets
    colnames(posnetmat) <- nets
    row.names(posnetmat) <- nets
    colnames(negnetmat) <- nets
    row.names(negnetmat) <- nets
    
    diffmat[upper.tri(diffmat)] <- 0
    
    dlim <- ifelse(abs(min(diffmat))>max(diffmat),abs(min(diffmat)),max(diffmat))
    
    colo <- colorRampPalette(c("skyblue2","white","darkorange2"))
    
    dev.new()
    corrplot::corrplot(diffmat,is.corr=FALSE,method="color",
                       tl.col="black",col = colo(100),na.label="square",
                       na.label.col = "white",addgrid.col="black",
                       title=paste("Difference in the Number of Edges\nin",con.title,sep=" "),
                       mar=c(0,0,4,0),cl.length=3,cl.pos="b",cl.lim=c(-dlim,dlim))
    
    if(visual.nets)
    {
        dev.new()
        qgraph::qgraph(posnetmat,title=paste("Positive",con.title,"Connectivity",sep=" "),
                       edge.color="darkorange2")
        dev.new()
        qgraph::qgraph(negnetmat,title=paste("Negative",con.title,"Connectivity",sep=" "),
                       edge.color="skyblue2")
        dev.new()
        diffmat <- (diffmat + t(diffmat))/2
        qgraph::qgraph(diffmat,title=paste("Difference",con.title,"Connectivity",sep=" "),
                       posCol="darkorange2",negCol="skyblue2")
        
        dev.new()
        qgraph::qgraph(poslobemat,title=paste("Positive",lobe.title,"Connectivity",sep=" "),
                       edge.color="darkorange2")
        dev.new()
        qgraph::qgraph(neglobemat,title=paste("Negative",lobe.title,"Connectivity",sep=" "),
                       edge.color="skyblue2")
        dev.new()
        ldiffmat <- (ldiffmat + t(ldiffmat))/2
        qgraph::qgraph(ldiffmat,title=paste("Difference",lobe.title,"Connectivity",sep=" "),
                       posCol="darkorange2",negCol="skyblue2")
    }
}
#----