#' Detects Node Crossings in a Network
#' 
#' @description UNDER DEVELOPMENT.
#' Computes \code{\link[NetworkToolbox]{rspbc}} for connections
#' between dimensions in a network. Multidimensional nodes can be detected 
#' 
#' @param A An adjacency matrix of network data
#' 
#' @param comm A vector or matrix corresponding to the
#' community each node belongs to
#' 
#' @param plot Should a plot be produced?
#' 
#' @return Produces a list containing:
#' 
#' @examples
#' # Uses Pearson's correlation 
#' tmfg <- TMFG(neoOpen)$A
#' 
#' \dontrun{
#' # Better to use polychoric correlations with this dataset
#' ega.glasso <- EGAnet::EGA(neoOpen)
#' 
#' result <- node.multidimensional(A = ega.glasso$network, comm = ega.glasso$wc, plot = FALSE)
#' 
#' }
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Node Mutildimensional Function
node.multidimensional <- function(A, comm, plot = FALSE)
{
    #Compute randomized shortest paths betweenness
    #with communities
    bc <- rspbc(A = A,comm = comm)
    
    #Compute between-dimension strength for each node and dimension
    cc <- comcat(A,comm,cent="strength",metric="each")
    
    #Compute RSPBC that accounts for strength of connections
    #between dimensions
    sbc <- (cc / rowSums(cc,na.rm=TRUE)) * bc
    
    #Set NA values to zero and round to zero decimal places
    res <- round(ifelse(is.na(sbc),0,sbc))
    
    #Plot?
    if(plot)
    {
        #Absolute
        res <- abs(res)
        
        #Standardize by maximum rspbc
        std.res <- res / max(rowSums(res))
        
        #Ensure that pie value is not greater than 1
        std.res <- std.res - .001
        std.res <- ifelse(std.res==-.001,0,std.res)
        
        #Split results to list for each node
        pies <- split(std.res, rep(1:nrow(std.res)))
        
        #Plot
        suppressWarnings(
            qgraph::qgraph(A, layout = "spring",
                       groups = as.factor(comm),
                       label.prop = 1.5,
                       pie = pies,
                       negDashed = TRUE,
                       vTrans = 200)
        )
    }
    
    return(res)
} 
#----