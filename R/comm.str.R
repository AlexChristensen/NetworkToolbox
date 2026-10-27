#' Community Strength/Degree Centrality
#' @description Computes the community
#' \code{\link[NetworkToolbox]{strength}}/\code{\link[NetworkToolbox]{degree}}
#' centrality measure of each community in a network or computes the
#' \code{\link[NetworkToolbox]{strength}}/\code{\link[NetworkToolbox]{degree}}
#' centrality measure of each community's connections to the other communities
#'
#' @param A An adjacency matrix of network data
#'
#' @param comm A vector corresponding to the
#' community each node belongs to
#'
#' @param weighted Is the network weighted?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} for weighted measures
#'
#' @param measure Type of measure to compute:
#'
#' \itemize{
#'
#' \item \code{"within"} --- Computes the community strength or degree of nodes within its own community
#'
#' \item \code{"between"} --- Computes the community strength or degree of nodes outside of its own community
#'
#' }
#'
#' @return A vector of community strength/degree centrality values for each specified
#' community in the network
#' (larger values suggest more central positioning)
#'
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(neoOpen, normal = FALSE)$A
#'
#' comm <- igraph::walktrap.community(convert2igraph(abs(A)))$membership
#'
#' #Strength
#' within.ns <- comm.str(A, comm, measure = "within")
#' between.ns <- comm.str(A, comm, measure = "between")
#'
#' #Degree
#' within.deg <- comm.str(A, comm, weighted = FALSE, measure = "within")
#' between.deg <- comm.str(A, comm, weighted = FALSE, measure = "between")
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#Community Strength/Degree Centrality
# Updated 26.11.2022
comm.str <- function (A, comm, weighted = TRUE,
                      measure = c("within","between"))
{
    if(is.null(comm))
    {stop("comm must be input")}

    if(missing(measure))
    {measure <- "between"
    }else{measure <- match.arg(measure)}

    comm <- as.vector(comm)

    if(ncol(A)!=length(comm))
    {stop("length of comm does not match nodes in matrix")}

    uniq <- unique(comm)
    len <- length(uniq)

    remove <- matrix(0,nrow=len,ncol=1)

    if(!weighted)
    {A <- binarize(A)}

    for(j in 1:len)
    {
        rem <- which(comm==uniq[j])

        if(measure == "within")
        {
            # a single node has no within-community connections
            if(length(rem)!=1)
            {remove[j,] <- sum(A[rem,rem,drop=FALSE])
            }else{remove[j,] <- 0}
        }else if(measure == "between")
        {remove[j,] <- sum(A[-rem,rem,drop=FALSE])}
    }

    norm <- as.vector(round(remove,3))

    names(norm) <- uniq

    # order labels numerically when possible
    num <- suppressWarnings(as.numeric(names(norm)))

    if(!anyNA(num))
    {res <- norm[order(num)]
    }else{res <- norm[sort(names(norm))]}

    return(res)
}
#----