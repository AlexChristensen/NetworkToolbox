#' Determines significant values for \code{\link[NetworkToolbox]{node.overlap}}
#' 
#' @description Identifies a signficant value for nodes that have a topological
#' overlap over a certain standard deviation above the
#' mean topological overlap (TO) value.
#' The function bins TO values by rounding to the nearest tenth. Then, the
#' local maximas are derived. Typically, there will be two. First, there
#' are 
#' 
#' @param vec Numeric vector.
#' Bootstrapped data from \code{\link[NetworkToolbox]{node.overlap}} function
#' 
#' @param sd.above Numeric.
#' Standard deviations above the mean that should be used to determine
#' the significance value.
#' Defaults to \code{2}
#' 
#' @return Returns a significance value determined from the bootstrapped
#' values of weighted topological overlap
#' 
#' @examples
#' 
#' #bimodal distribution
#' nn <- 1e4
#' betas <- rbeta(nn,2,2)
#' sims <- c(betas[1:(nn/20)]*2+3, betas[(nn/20+1):nn]*2+1)
#' 
#' #distribution akin to bootstrapped weighted topological overlap
#' wto.dist <- na.omit(ifelse(sims<=2,NA,sims))*.10
#' hist(wto.dist)
#' 
#' #significance value of the right distribution
#' result <- node.overlap.wrapper(wto.dist)
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Node Overlap Function
node.overlap.wrapper <- function (vec, sd.above = 2)
{
    #Local maxima function
    localMaxima <- function(x) {
        # Use -Inf instead if x is numeric (non-integer)
        y <- diff(c(-.Machine$integer.max, x)) > 0L
        rle(y)$lengths
        y <- cumsum(rle(y)$lengths)
        y <- y[seq.int(1L, length(y), 2L)]
        if (x[[1]] == x[[2]]) {
            y <- y[-1]
        }
        y
    }
    
    #round values to place into discrete bins
    round.val <- round(vec, 1)
    
    #unique values
    uniq.val <- unique(round.val)
    
    #ordered unique values
    uniq.ord <- uniq.val[order(uniq.val)]
    
    #initialize frequency counts vector
    counts <- vector("numeric", length(uniq.ord))
    
    #compute counts
    for(i in 1:length(counts))
    {counts[i] <- length(which(round.val==uniq.ord[i]))}
    
    #transfer ordered unique values to counts
    names(counts) <- uniq.ord
    
    #identify lowest minimum local maximum
    vals <- min(as.numeric(names(localMaxima(counts))))
    
    #determine cut point to identify normal distribution
    cut.point <- as.numeric(names(counts[which(names(counts)==vals)+1]))
    
    #omit values below cut point
    clean.val <- na.omit(ifelse(vec<=cut.point,NA,vec))
    
    #identify significance value
    sig.val <- mean(clean.val) + (sd.above*sd(clean.val))
    
    return(sig.val)
}
