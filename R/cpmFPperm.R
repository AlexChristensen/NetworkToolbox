#' Connectome-based Predictive Modeling--Fingerprinting Permutation
#' @description Applies the Connectome-based Predictive Modeling approach to neural data.
#' This method identifies individuals based on their specific connectivity patterns.
#' \strong{Please cite Finn et al., 2015; Rosenberg et al., 2016; Shen et al., 2017}
#' 
#' @param session1 Array from \code{\link{convertConnBrainMat}} function
#' (first session)
#' 
#' @param session2 Array from \code{\link{convertConnBrainMat}} function
#' (second session)
#' 
#' @param iter Number of iterations to perform.
#' Defaults to \code{1000}
#' 
#' @param progBar Should progress bar be displayed?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} for no progress bar
#' 
#' @return Returns a matrix containing the percentage and
#' number of correctly identified subjects for sessions 1 and 2
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
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
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