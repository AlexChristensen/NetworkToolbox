#' Connectome-based Predictive Modeling--Internal Validation
#' @description Applies the Connectome-based Predictive Modeling approach to neural data.
#' This method predicts a behavioral statistic using neural connectivity from the sample.
#' \strong{Please cite Finn et al., 2015; Rosenberg et al., 2016; Shen et al., 2017}
#' 
#' @param neuralarray Array from \code{\link{convertConnBrainMat}} function
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
#' @param corr Correlation method for assessing the relatonship between the behavioral measure and edges between ROIs.
#' Defaults to \code{"pearson"}.
#' Set to \code{"spearman"} for non-linear or monotonic associations
#' 
#' @param shen Are ROIs from Shen et al. 2013 atlas?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} for canonical networks plot
#' 
#' @param cores Number of computer processing cores to use when performing covariate analyses.
#' Defaults to \emph{n} - 1 total number of cores.
#' Set to any number between 1 and maxmimum amount of cores on your computer
#' 
#' @param progBar Should progress bar be displayed?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} for no progress bar
#' 
#' @return Returns a list containing: 
#'
#' \item{results}{A matrix contaning: r coefficient (\code{r}), p-value (\code{p-value}),
#' mean absolute error (\code{mae}), root mean square error (\code{rmse})}
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
#' @importFrom graphics par abline hist plot text
#' @importFrom grDevices colorRampPalette dev.new
#' @importFrom foreach %dopar%
#' @importFrom utils menu
#' @importFrom stats cor.test coef cov2cor lm na.omit
#' 
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
        colnames(poslobemat)<-c("PFC","Mot","Ins","Par","Tem","Occ","Lim","Cer","Sub","Bsm")
        row.names(poslobemat)<-c("PFC","Mot","Ins","Par","Tem","Occ","Lim","Cer","Sub","Bsm")
        colnames(neglobemat)<-c("PFC","Mot","Ins","Par","Tem","Occ","Lim","Cer","Sub","Bsm")
        row.names(neglobemat)<-c("PFC","Mot","Ins","Par","Tem","Occ","Lim","Cer","Sub","Bsm")
        
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
        colnames(posnetmat)<-c("MF","FP","DM","SubC","MT","VI","VII","VA")
        row.names(posnetmat)<-c("MF","FP","DM","SubC","MT","VI","VII","VA")
        colnames(negnetmat)<-c("MF","FP","DM","SubC","MT","VI","VII","VA")
        row.names(negnetmat)<-c("MF","FP","DM","SubC","MT","VI","VII","VA")
        
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
    
    if(shen)
    {
        ans <- menu(c("Yes","No"),title="Visualize canonical and macro-scale regions connectivity?")
        
        if(ans==1)
        {
            dev.new()
            qgraph::qgraph(posnetmat,title="Positive Canonical Connectivity",
                           edge.color="darkorange2")
            dev.new()
            qgraph::qgraph(negnetmat,title="Negative Canonical Connectivity",
                           edge.color="skyblue2")
            dev.new()
            diffmat <- (diffmat + t(diffmat))/2
            qgraph::qgraph(diffmat,title="Difference Canonical Connectivity",
                           posCol="darkorange2",negCol="skyblue2")
            
            dev.new()
            qgraph::qgraph(poslobemat,title="Positive Macro-scale Regions Connectivity",
                           edge.color="darkorange2")
            dev.new()
            qgraph::qgraph(neglobemat,title="Negative Macro-scale Regions Connectivity",
                           edge.color="skyblue2")
            dev.new()
            ldiffmat <- (ldiffmat + t(ldiffmat))/2
            qgraph::qgraph(ldiffmat,title="Difference Macro-scale Regions Connectivity",
                           posCol="darkorange2",negCol="skyblue2")
        }
    }
    
    if(shen)
    {return(list(results=results,posMask=posmask,negMask=negmask,
                 posCanon=posnetmat,negCanon=negnetmat,
                 posMacro=poslobemat,negMacro=neglobemat))
    }else{return(list(results=results,posMask=posmask,negMask=negmask))}
}
#----