#' Plots CPM results
#' 
#' @param x A \code{\link[NetworkToolbox]{cpm}} object
#' 
#' @param ... Additional arguments for \link{plot}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
# General plot function for CPM
# Updated 10.09.2020
plot.cpm <- function(x, ...)
{
  groups <- x$groups
  
  plot.args <- list(...)
  
  if(length(plot.args) == 0)
  {plot.args$type <- "p"}
  
  if(x$connections == "separate")
  {
    if(!"main" %in% names(plot.args))
    {mains <- c("Positive Prediction", "Negative Prediction")
    }else{mains <- plot.args$main}
  }else{
    if(!"main" %in% names(plot.args))
    {mains <- c("Overall Prediction")
    }else{mains <- plot.args$main}
  }
  
  if(!"xlab" %in% names(plot.args))
  {plot.args$xlab <- "Observed Score\n(Z-score)"}
  
  if(!"ylab" %in% names(plot.args))
  {plot.args$ylab <- "Predicted Score\n(Z-score)"}
  
  if(!"col" %in% names(plot.args))
  {
    cols <- ifelse(is.null(groups), c("darkorange2", "skyblue2"), c("darkorange2", "darkorange2", "skyblue2", "skyblue2"))
    
  }else{
    
    if(length(plot.args$col) == 2)
    {cols <- rep(c(plot.args$col[1], plot.args$col[2]), 2)
    }else{cols <- plot.args$col}
    
  }
  
  if(!"pch"%in% names(plot.args))
  {
    pchs <- ifelse(is.null(groups), c(16,16), c(16,1,16,1))
    
  }else{
    
    if(length(plot.args$pch) == 2)
    {pchs <- rep(c(plot.args$pch[1], plot.args$pch[2]), 2)
    }else{pchs <- plot.args$pch}
    
  }
  
  # Check for missing predictions
  
  if(x$connections == "separate")
  {
    if(any(is.na(x$posPred)))
    {
      pos.na <- which(is.na(x$posPred))
      bstat_pos <- x$behav[-pos.na]
      behav_pred_pos <- x$posPred[-pos.na]
    }else{
      bstat_pos <- x$behav
      behav_pred_pos <- x$posPred
    }
    
    if(any(is.na(x$negPred)))
    {
      neg.na <- which(is.na(x$negPred))
      bstat_neg <- x$behav[-neg.na]
      behav_pred_neg <- x$posPred[-neg.na]
    }else{
      bstat_neg <- x$behav
      behav_pred_neg <- x$negPred
    }
    
    R_pos<-cor(behav_pred_pos,bstat_pos,use="pairwise.complete.obs")
    P_pos<-cor.test(behav_pred_pos,bstat_pos)$p.value
    R_neg<-cor(behav_pred_neg,bstat_neg,use="pairwise.complete.obs")
    P_neg<-cor.test(behav_pred_neg,bstat_neg)$p.value
    
    P_pos<-ifelse(round(P_pos,3)!=0,round(P_pos,3),noquote("< .001"))
    P_neg<-ifelse(round(P_neg,3)!=0,round(P_neg,3),noquote("< .001"))
    
    #bstat range
    lower.bstat_pos <- floor(range(bstat_pos))[1]
    upper.bstat_pos <- ceiling(range(bstat_pos))[2]
    lower.bstat_neg <- floor(range(bstat_neg))[1]
    upper.bstat_neg <- ceiling(range(bstat_neg))[2]
    lower.pos.pred <- floor(range(behav_pred_pos))[1]
    upper.pos.pred <- ceiling(range(behav_pred_pos))[2]
    lower.neg.pred <- floor(range(behav_pred_neg))[1]
    upper.neg.pred <- ceiling(range(behav_pred_neg))[2]
    
    text.one_pos <- lower.bstat_pos - (lower.bstat_pos * .20)
    text.one_neg <- lower.bstat_neg - (lower.bstat_neg * .20)
    
    #set up groups
    if(!is.null(groups))
    {
      #group labels
      labs_groups <- unique(groups)
      
      #number of groups
      n_groups <- length(labs_groups)
      
      #set up colors
      plot.args$col <- c(rep(cols[1],length(which(groups == labs_groups[1]))),
                         rep(cols[2],length(which(groups == labs_groups[2]))))
      
      #set up points
      plot.args$pch <- c(rep(pchs[1],length(which(groups == labs_groups[1]))),
                         rep(pchs[2],length(which(groups == labs_groups[2]))))
      
      #set up x
      plot.args$x <- bstat_pos
      
      #set up y
      plot.args$y <- behav_pred_pos
      
      #set up ylim
      plot.args$ylim <- c(lower.pos.pred, upper.pos.pred)
      
      #set up xlim
      plot.args$xlim <- c(lower.bstat_pos, upper.bstat_pos)
      
      #set up main
      plot.args$main <- mains[1]
      
      #plot positive
      dev.new()
      par(mar=c(5,5,4,2))
      do.call(plot, plot.args)
      
      abline(lm(behav_pred_pos~bstat_pos))
      
      if(R_pos>=0)
      {
        text.two <- upper.pos.pred - (upper.pos.pred * .20)
        
        text(x = text.one_pos, y = text.two, labels = paste("r = ", round(R_pos,3), "\np = ", P_pos))
        
        legend("bottomright", legend = labs_groups, col = c(cols[1], cols[2]), pch = c(pchs[1], pchs[2]))
        
      }else if(R_pos<0)
      {
        text.two <- lower.pos.pred - (lower.pos.pred * .20)
        
        text(x = text.one_pos, y = text.two, labels = paste("r = ", round(R_pos,3), "\np = ", P_pos))
        
        legend("topright", legend = labs_groups, col = c(cols[1], cols[2]), pch = c(pchs[1], pchs[2]))
      }
      
      #set up colors
      plot.args$col <- c(rep(cols[3],length(which(groups == labs_groups[1]))),
                         rep(cols[4],length(which(groups == labs_groups[2]))))
      
      #set up points
      plot.args$pch <- c(rep(pchs[3],length(which(groups == labs_groups[1]))),
                         rep(pchs[4],length(which(groups == labs_groups[2]))))
      
      #set up x
      plot.args$x <- bstat_neg
      
      #set up y
      plot.args$y <- behav_pred_neg
      
      #set up ylim
      plot.args$ylim <- c(lower.neg.pred, upper.neg.pred)
      
      #set up xlim
      plot.args$xlim <- c(lower.bstat_neg, upper.bstat_neg)
      
      #set up main
      plot.args$main <- mains[2]
      
      #plot negative
      dev.new()
      par(mar=c(5,5,4,2))
      do.call(plot, plot.args)
      
      abline(lm(behav_pred_neg~bstat_neg))
      
      if(R_neg>=0)
      {
        text.two <- upper.neg.pred - (upper.neg.pred * .20)
        
        text(x = text.one_neg, y = text.two, labels = paste("r = ", round(R_neg,3), "\np = ", P_neg))
        
        legend("bottomright", legend = labs_groups, col = c(cols[3], cols[4]), pch = c(pchs[3], pchs[4]))
        
      }else if(R_neg<0)
      {
        text.two <- lower.neg.pred - (lower.neg.pred * .20)
        
        text(x = text.one_neg, y = text.two, labels = paste("r = ", round(R_neg,3), "\np = ", P_neg))
        
        legend("topright", legend = labs_groups, col = c(cols[3], cols[4]), pch = c(pchs[3], pchs[4]))
        
      }
      
    }else{
      
      #set up colors
      plot.args$col <- cols[1]
      
      #set up points
      plot.args$pch <- pchs[1]
      
      #set up x
      plot.args$x <- bstat_pos
      
      #set up y
      plot.args$y <- behav_pred_pos
      
      #set up ylim
      plot.args$ylim <- c(lower.pos.pred, upper.pos.pred)
      
      #set up xlim
      plot.args$xlim <- c(lower.bstat_pos, upper.bstat_pos)
      
      #set up main
      plot.args$main <- mains[1]
      
      #plot positive
      dev.new()
      par(mar=c(5,5,4,2))
      do.call(plot, plot.args)
      
      abline(lm(behav_pred_pos~bstat_pos))
      
      if(R_pos>=0)
      {
        text.two <- upper.pos.pred - (upper.pos.pred * .20)
        
        text(x = text.one_pos, y = text.two, labels = paste("r = ", round(R_pos,3), "\np = ", P_pos))
        
      }else if(R_pos<0)
      {
        text.two <- lower.pos.pred - (lower.pos.pred * .20)
        
        text(x = text.one_pos, y = text.two, labels = paste("r = ", round(R_pos,3), "\np = ", P_pos))
      }
      
      #set up colors
      plot.args$col <- cols[2]
      
      #set up points
      plot.args$pch <- pchs[2]
      
      #set up x
      plot.args$x <- bstat_neg
      
      #set up y
      plot.args$y <- behav_pred_neg
      
      #set up ylim
      plot.args$ylim <- c(lower.neg.pred, upper.neg.pred)
      
      #set up xlim
      plot.args$xlim <- c(lower.pos.pred, upper.pos.pred)
      
      #set up main
      plot.args$main <- mains[2]
      
      #plot negative
      dev.new()
      par(mar=c(5,5,4,2))
      do.call(plot, plot.args)
      
      abline(lm(behav_pred_neg~bstat_neg))
      
      if(R_neg>=0)
      {
        text.two <- upper.neg.pred - (upper.neg.pred * .20)
        
        text(x = text.one_neg, y = text.two, labels = paste("r = ", round(R_neg,3), "\np = ", P_neg))
        
      }else if(R_neg<0)
      {
        text.two <- lower.neg.pred - (lower.neg.pred * .20)
        
        text(x = text.one_neg, y= text.two, labels = paste("r = ", round(R_neg,3), "\np = ", P_neg))
      }
    }
  }else{
    
    if(any(is.na(x$Pred)))
    {
      pos.na <- which(is.na(x$Pred))
      bstat_pos <- x$behav[-pos.na]
      behav_pred_pos <- x$Pred[-pos.na]
    }else{
      bstat_pos <- x$behav
      behav_pred_pos <- x$Pred
    }
    
    R_pos<-cor(behav_pred_pos,bstat_pos,use="pairwise.complete.obs")
    P_pos<-cor.test(behav_pred_pos,bstat_pos)$p.value
    
    P_pos<-ifelse(round(P_pos,3)!=0,round(P_pos,3),noquote("< .001"))
    
    #bstat range
    lower.bstat_pos <- floor(range(bstat_pos))[1]
    upper.bstat_pos <- ceiling(range(bstat_pos))[2]
    lower.pos.pred <- floor(range(behav_pred_pos))[1]
    upper.pos.pred <- ceiling(range(behav_pred_pos))[2]
    
    text.one_pos <- lower.bstat_pos - (lower.bstat_pos * .20)
    
    #set up groups
    if(!is.null(groups))
    {
      #group labels
      labs_groups <- unique(groups)
      
      #number of groups
      n_groups <- length(labs_groups)
      
      #set up colors
      plot.args$col <- c(rep(cols[1],length(which(groups == labs_groups[1]))),
                         rep(cols[2],length(which(groups == labs_groups[2]))))
      
      #set up points
      plot.args$pch <- c(rep(pchs[1],length(which(groups == labs_groups[1]))),
                         rep(pchs[2],length(which(groups == labs_groups[2]))))
      
      #set up x
      plot.args$x <- bstat_pos
      
      #set up y
      plot.args$y <- behav_pred_pos
      
      #set up ylim
      plot.args$ylim <- c(lower.pos.pred, upper.pos.pred)
      
      #set up xlim
      plot.args$xlim <- c(lower.bstat_pos, upper.bstat_pos)
      
      #set up main
      plot.args$main <- mains[1]
      
      #plot positive
      dev.new()
      par(mar=c(5,5,4,2))
      do.call(plot, plot.args)
      
      abline(lm(behav_pred_pos~bstat_pos))
      
      if(R_pos>=0)
      {
        text.two <- upper.pos.pred - (upper.pos.pred * .20)
        
        text(x = text.one_pos, y = text.two, labels = paste("r = ", round(R_pos,3), "\np = ", P_pos))
        
        legend("bottomright", legend = labs_groups, col = c(cols[1], cols[2]), pch = c(pchs[1], pchs[2]))
        
      }else if(R_pos<0)
      {
        text.two <- lower.pos.pred - (lower.pos.pred * .20)
        
        text(x = text.one_pos, y = text.two, labels = paste("r = ", round(R_pos,3), "\np = ", P_pos))
        
        legend("topright", legend = labs_groups, col = c(cols[1], cols[2]), pch = c(pchs[1], pchs[2]))
      }
      
    }else{
      
      #set up colors
      plot.args$col <- cols[1]
      
      #set up points
      plot.args$pch <- pchs[1]
      
      #set up x
      plot.args$x <- bstat_pos
      
      #set up y
      plot.args$y <- behav_pred_pos
      
      #set up ylim
      plot.args$ylim <- c(lower.pos.pred, upper.pos.pred)
      
      #set up xlim
      plot.args$xlim <- c(lower.bstat_pos, upper.bstat_pos)
      
      #set up main
      plot.args$main <- mains[1]
      
      #plot positive
      dev.new()
      par(mar=c(5,5,4,2))
      do.call(plot, plot.args)
      
      abline(lm(behav_pred_pos~bstat_pos))
      
      if(R_pos>=0)
      {
        text.two <- upper.pos.pred - (upper.pos.pred * .20)
        
        text(x = text.one_pos, y = text.two, labels = paste("r = ", round(R_pos,3), "\np = ", P_pos))
        
      }else if(R_pos<0)
      {
        text.two <- lower.pos.pred - (lower.pos.pred * .20)
        
        text(x = text.one_pos, y = text.two, labels = paste("r = ", round(R_pos,3), "\np = ", P_pos))
      }
    }
  }
}
