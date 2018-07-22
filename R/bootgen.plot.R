#' Bootstrapped Network Generalization Plots
#' @description Generates reliability plots from the \link[NetworkToolbox]{bootgen} function
#' 
#' @param object An output from the \link[NetworkToolbox]{bootgen} function
#' 
#' @param bootmat Should the bootstrap generalization matrix be included?
#' Defaults to FALSE.
#' Set to TRUE to plot the original and bootmat networks' reliablities
#' 
#' @param breaks If bootmat = TRUE, then reliability comparison histogram is produced.
#' The breaks may not appear as desired, so the researcher can adjust as needed.
#' Defaults to 20 (try 10 if 20 isn't equivalent breaks)
#' 
#' @return Returns the following plots:
#' 
#' \item{orignet}{A reliability matrix for the original network
#' ("Original Network Correlation Reliabilities"
#' upper triangle = network reliabilites from the filtered network,
#' lower triangle = full, unfiltered matrix reliablities)}
#' 
#' \item{orignet retained}{A plot of retained correlations on their reliability
#' ("Original Network Correlation Strength on Reliability")}
#' 
#' \item{bootmat}{A reliability matrix for the bootgen network
#' ("bootgen Network Correlation Reliabilities"
#' upper triangle = network reliabilites from the bootgen filtered network,
#' lower triangle = full, unfiltered matrix reliablities)}
#' 
#' \item{bootmat retained}{A plot of retained correlations on their reliability
#' ("bootgen Network Correlation Strength on Reliability")}
#' 
#' \item{histogram}{A frequency of the reliability of edges
#' in both the original and bootgen network}
#' 
#' @examples
#' \dontrun{
#' bootTMFG <- bootgen(neoOpen)
#' 
#' bootPlot <- bootPlot(bootTMFG)
#' }
#' @references
#' 
#' Wei, T. & Simko, V.(2017).
#' R package "corrplot": Visualization of a correlation matrix (Version 0.84).
#' Available from \url{https://github.com/taiyun/corrplot}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @importFrom graphics abline plot text hist legend
#' 
#' @importFrom grDevices rgb
#' 
#' @importFrom stats lm na.omit
#' 
#' @export
#Bootstrap Network Generalization Plots----
bootgen.plot <- function (object, bootmat = FALSE, breaks = 20)
{
    tru<-object$orignet
    rel<-object$netrel
    if(bootmat==TRUE)
    {
        boot<-object$bootmat
        brel<-rel
    }
    
    #create tru reliability plot matrix
    row.names(rel)<-colnames(rel)
    diag(rel)<-1
    upp<-matrix(0,nrow=nrow(rel),ncol=ncol(rel))
    for(i in 1:nrow(rel))
        for(j in 1:ncol(rel))
            if(rel[i,j]!=0&&tru[i,j]!=0)
            {upp[i,j]<-rel[i,j]}
    colnames(upp)<-colnames(rel)
    rel[upper.tri(rel)]<-upp[upper.tri(upp)]
    
    #reliablity on correlation plot
    x<-matrix(nrow=length(upp))
    y<-matrix(nrow=length(tru))
    diag(upp)<-0
    diag(tru)<-0
    wc<-0
    for(i in 1:nrow(tru))
        for(j in 1:ncol(tru))
            if((upp[i,j]!=0&&tru[i,j])!=0)
            {wc<-wc+1
            x[wc]<-upp[i,j]
            y[wc]<-tru[i,j]}
    xo<-na.omit(abs(x))
    yo<-na.omit(abs(y))
    
    dev.new()
    mar=c(2,2,2,2)
    cpo<-{plot(xo,yo,pch=16,ylab="Edge Weight Strength",xlab="Reliability",
               main="Original Network Edge Weight Strength on Reliability",xlim=c(0,1),ylim=range(yo))
        abline(lm(yo~xo))
        text(x=.05,y=max(yo-.05),labels = paste("r = ",round(cor(yo,xo),3)))}
    
    
    #plot reliability matrix
    dev.new()
    if(ncol(tru)<=20)
    {plt<-corrplot::corrplot(rel,method="color",
                             title="Original Network Edge Weight Reliabilities",
                             mar=c(2,2,2,2),tl.col="black",tl.cex=.75,
                             cl.lim = c(0,1),addgrid.col = "grey",addCoef.col = "black")
    }else if(ncol(tru)>20){
        plt<-corrplot::corrplot(rel,method="color",
                                title="Original Network Edge Weight Reliabilities",
                                mar=c(2,2,2,2),tl.col="black",tl.cex=.75,
                                cl.lim = c(0,1),addgrid.col = "grey")}
    
    
    ############################################################################################
    #create boot reliability plot matrix
    if(bootmat==TRUE)
    {
        row.names(brel)<-colnames(brel)
        diag(brel)<-1
        upp<-matrix(0,nrow=nrow(brel),ncol=ncol(brel))
        for(i in 1:nrow(brel))
            for(j in 1:ncol(brel))
                if(brel[i,j]!=0&&boot[i,j]!=0)
                {upp[i,j]<-brel[i,j]}
        colnames(upp)<-colnames(brel)
        brel[upper.tri(brel)]<-upp[upper.tri(upp)]
        
        #reliablity on correlation plot
        x<-matrix(nrow=length(upp))
        y<-matrix(nrow=length(boot))
        diag(upp)<-0
        diag(boot)<-0
        wc<-0
        for(i in 1:nrow(boot))
            for(j in 1:ncol(boot))
                if((upp[i,j]!=0&&boot[i,j])!=0)
                {wc<-wc+1
                x[wc]<-upp[i,j]
                y[wc]<-boot[i,j]}
        xo<-na.omit(abs(x))
        yo<-na.omit(abs(y))
        
        dev.new()
        mar=c(2,2,2,2)
        cpo<-{plot(xo,yo,pch=16,ylab="Edge Weight Strength",xlab="Reliability",
                   main="bootgen Network Edge Weight Strength on Reliability",xlim=c(0,1),ylim=range(yo))
            abline(lm(yo~xo))
            text(x=.05,y=max(yo-.05),labels = paste("r = ",round(cor(yo,xo),3)))}
        
        
        #plot reliability matrix
        dev.new()
        if(ncol(boot)<=20)
        {plt<-corrplot::corrplot(brel,method="color",
                                 title="bootgen Network Edge Weight Reliabilities",
                                 mar=c(2,2,2,2),tl.col="black",tl.cex=.75,
                                 cl.lim = c(0,1),addgrid.col = "grey",addCoef.col = "black")
        }else if(ncol(boot)>20){
            plt<-corrplot::corrplot(brel,method="color",
                                    title="bootgen Network Edge Weight Reliabilities",
                                    mar=c(2,2,2,2),tl.col="black",tl.cex=.75,
                                    cl.lim = c(0,1),addgrid.col = "grey")}
        
        dev.new()
        hist(na.omit(as.vector(ifelse(tru!=0,rel,NA))),
             breaks=breaks,xlab="Reliability",ylab="Number of Edges",main="Reliability Frequency",xlim=c(0,1),
             col = rgb(0,0,1,0.5))
        hist(na.omit(as.vector(ifelse(boot!=0,brel,NA))),breaks=breaks,
             add=TRUE, col = rgb(1,0,0,0.5))
        legend("topleft",c("orignet","bootgen","overlap"),fill=c(rgb(0,0,1,0.5),rgb(1,0,0,0.5),"maroon3"))
        
        
    }
}
#----