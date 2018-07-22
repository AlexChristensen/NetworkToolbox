#' Core Items
#' @description Automatically determines core, intermediary, and peripheral items in the network.
#' The entire network or within-community gradations can be determined
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param comm A vector or matrix corresponding to the community each node belongs to
#' 
#' @param by Should the core items be defined by network or communities?
#' Defaults to "network".
#' Set to "communities" to define core items within communities
#' 
#' @return Returns a list containing:
#' 
#' \item{core}{Core items for each community}
#' 
#' \item{inter}{Intermediate items for each community}
#' 
#' \item{peri}{Peripheral items for each community}
#' 
#' @examples
#' #network
#' A <- TMFG(neoOpen)$A
#' 
#' #core items by network
#' coreBYnetwork <- core.items(A, by = "network")
#' 
#' #theoretical factors
#' comm <- c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,8),rep(6,8))
#' 
#' #core items by communities
#' coreBYcomm <- core.items(A, comm, by = "communities")
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Core Items----
core.items <- function (A, comm, by = c("network","communities"))
{
    if(missing(comm))
    {by <- "network"
    }else{comm <- comm}
    
    if(missing(by))
    {by <- "communities"
    }else{by <- match.arg(by)}
    
    gradlist <- list()
    
    hc <- hybrid(A, BC = "random")
    
    left <- length(hc)
    
    name <- names(hc)
    
    if(by == "network")
    {
        dat <- hc
        
        corenum <- floor(left/3)
    
        cst <- 1:corenum
    
        left <- left - corenum
    
        internum <- floor(left/2)
    
        ist <- (corenum+1):(corenum+internum)
    
        left <- left - internum
    
        perinum <- left
    
        pst <- (corenum+internum+1):(corenum+internum+perinum)
    
        ord <- dat[order(dat,decreasing=TRUE)]
        
        gradlist$core <- ord[cst]
        gradlist$inter <- ord[ist]
        gradlist$peri <- ord[pst]
    }else if(by == "communities")
    {
        uniq <- unique(comm)
        if(is.numeric(uniq))
        {uniq <- as.character(uniq)}
        
        len <- length(uniq)
        
        for(i in 1:len)
        {
            left <- length(which(uniq[i]==comm))
            
            corenum <- floor(left/3)
            
            cst <- 1:corenum
            
            left <- left - corenum
            
            internum <- floor(left/2)
            
            ist <- (corenum+1):(corenum+internum)
            
            left <- left - internum
            
            perinum <- left
            
            pst <- (corenum+internum+1):(corenum+internum+perinum)
            
            items <- which(uniq[i]==comm)
            
            #core
            core <- hc[items][order(hc[items],decreasing = TRUE)][cst]
            gradlist[[uniq[i]]]$core <- core
            
            #inter
            inter <- hc[items][order(hc[items],decreasing = TRUE)][ist]
            gradlist[[uniq[i]]]$inter <- inter
            
            #peri
            peri <- hc[items][order(hc[items],decreasing = TRUE)][pst]
            gradlist[[uniq[i]]]$peri <- peri
        }
    }
    
    return(gradlist)
}
#----