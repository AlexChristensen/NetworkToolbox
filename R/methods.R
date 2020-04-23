# General plot function for CPM
# Updated 22.04.2020
plot.cpm <- function(x, ...)
{
  bstat <- x$behav
  behav_pred_pos <- x$posPred
  behav_pred_neg <- x$negPred
  groups <- x$groups
  
  R_pos<-cor(behav_pred_pos,bstat,use="pairwise.complete.obs")
  P_pos<-cor.test(behav_pred_pos,bstat)$p.value
  R_neg<-cor(behav_pred_neg,bstat,use="pairwise.complete.obs")
  P_neg<-cor.test(behav_pred_neg,bstat)$p.value
  
  P_pos<-ifelse(round(P_pos,3)!=0,round(P_pos,3),noquote("< .001"))
  P_neg<-ifelse(round(P_neg,3)!=0,round(P_neg,3),noquote("< .001"))
  
  #bstat range
  lower.bstat <- floor(range(bstat))[1]
  upper.bstat <- ceiling(range(bstat))[2]
  lower.pos.pred <- floor(range(behav_pred_pos))[1]
  upper.pos.pred <- ceiling(range(behav_pred_pos))[2]
  lower.neg.pred <- floor(range(behav_pred_neg))[1]
  upper.neg.pred <- ceiling(range(behav_pred_neg))[2]
  
  text.one <- lower.bstat - (lower.bstat * .20)
  
  #set up groups
  if(!is.null(groups))
  {
    #group labels
    labs_groups <- unique(groups)
    
    #number of groups
    n_groups <- length(labs_groups)
    
    #plot positive
    dev.new()
    par(mar=c(5,5,4,2))
    plot(bstat,behav_pred_pos,xlab="Observed Score\n(Z-score)",ylab="Predicted Score\n(Z-score)",
         main="Positive Prediction",xlim=c(lower.bstat,upper.bstat),
         ylim=c(lower.pos.pred,upper.pos.pred),pch=c(rep(16,length(which(groups == labs_groups[1]))),
                                                     rep(1,length(which(groups == labs_groups[2])))),col="darkorange2")
    abline(lm(behav_pred_pos~bstat))
    if(R_pos>=0)
    {
      text.two <- upper.pos.pred - (upper.pos.pred * .20)
      
      text(x=text.one,y=text.two,labels = paste("r = ",round(R_pos,3),"\np = ",P_pos))
      legend("bottomright",legend=labs_groups,col="darkorange2",pch=c(16,1))
    }else if(R_pos<0)
    {
      text.two <- lower.pos.pred - (lower.pos.pred * .20)
      
      text(x=text.one,y=text.two,labels = paste("r = ",round(R_pos,3),"\np = ",P_pos))
      legend("topright",legend=labs_groups,col="darkorange2",pch=c(16,1))
    }
    
    #plot negative
    dev.new()
    par(mar=c(5,5,4,2))
    plot(bstat,behav_pred_neg,xlab="Observed Score\n(Z-score)",ylab="Predicted Score\n(Z-score)",
         main="Negative Prediction",xlim=c(lower.bstat,upper.bstat),
         ylim=c(lower.neg.pred,upper.neg.pred),pch=c(rep(16,length(which(groups == labs_groups[1]))),
                                                     rep(1,length(which(groups == labs_groups[2])))),col="skyblue2")
    abline(lm(behav_pred_neg~bstat))
    if(R_neg>=0)
    {
      text.two <- upper.neg.pred - (upper.neg.pred * .20)
      
      text(x=text.one,y=text.two,labels = paste("r = ",round(R_neg,3),"\np = ",P_neg))
      legend("bottomright",legend=labs_groups,col="skyblue2",pch=c(16,1))
    }else if(R_neg<0)
    {
      text.two <- lower.neg.pred - (lower.neg.pred * .20)
      
      text(x=text.one,y=text.two,labels = paste("r = ",round(R_neg,3),"\np = ",P_neg))
      legend("topright",legend=labs_groups,col="skyblue2",pch=c(16,1))
    }
  }else{
    #plot positive
    dev.new()
    par(mar=c(5,5,4,2))
    plot(bstat,behav_pred_pos,xlab="Observed Score\n(Z-score)",ylab="Predicted Score\n(Z-score)",
         main="Positive Prediction",xlim=c(lower.bstat,upper.bstat),
         ylim=c(lower.pos.pred,upper.pos.pred),pch=16,col="darkorange2")
    abline(lm(behav_pred_pos~bstat))
    if(R_pos>=0)
    {
      text.two <- upper.pos.pred - (upper.pos.pred * .20)
      
      text(x=text.one,y=text.two,labels = paste("r = ",round(R_pos,3),"\np = ",P_pos))
      
    }else if(R_pos<0)
    {
      text.two <- lower.pos.pred - (lower.pos.pred * .20)
      
      text(x=text.one,y=text.two,labels = paste("r = ",round(R_pos,3),"\np = ",P_pos))
    }
    
    #plot negative
    dev.new()
    par(mar=c(5,5,4,2))
    plot(bstat,behav_pred_neg,xlab="Observed Score\n(Z-score)",ylab="Predicted Score\n(Z-score)",
         main="Negative Prediction",xlim=c(lower.bstat,upper.bstat),
         ylim=c(lower.neg.pred,upper.neg.pred),pch=16,col="skyblue2")
    abline(lm(behav_pred_neg~bstat))
    if(R_neg>=0)
    {
      text.two <- upper.neg.pred - (upper.neg.pred * .20)
      
      text(x=text.one,y=text.two,labels = paste("r = ",round(R_neg,3),"\np = ",P_neg))
      
    }else if(R_neg<0)
    {
      text.two <- lower.neg.pred - (lower.neg.pred * .20)
      
      text(x=text.one,y=-text.two,labels = paste("r = ",round(R_neg,3),"\np = ",P_neg))
    }
  }
}