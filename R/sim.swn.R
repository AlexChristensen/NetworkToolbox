#' Simulate Small-world Network
#' @description Simulates a small-world network based on specified topological properties.
#' Data will also be simulated based on the true network structure
#'
#' @param nodes Number of nodes in the simulated network
#'
#' @param n Number of cases in the simulated dataset
#'
#' @param pos Proportion of positive correlations in the simulated network
#'
#' @param ran Range of correlations in the simulated network
#'
#' @param nei Adjusts the number of connections each node has to
#' neighboring nodes (see \code{\link[igraph]{sample_smallworld}})
#'
#' @param p Adjusts the rewiring probability (default is .5).
#' p > .5 rewires the simulated network closer to a random network.
#' p < .5 rewires the simulated network closer to a lattice network
#'
#' @param corr Should the simulated network be a correlation network?
#' Defaults to FALSE.
#' Set to TRUE for a simulated correlation network
#'
#' @param replace If noise > 0, then should participants be sampled with replacement?
#' Defaults to TRUE.
#' Set to FALSE to not allow the potential for participants to be consecutively entered
#' into the simulated dataset.
#'
#' @param ordinal Should simulated continuous data be converted to ordinal?
#' Defaults to FALSE.
#' Set to TRUE for simulated ordinal data
#'
#' @param ordLevels If ordinal = TRUE, then how many levels should be used?
#' Defaults to NULL.
#' Set to desired number of intervals (defaults to 5)
#'
#' @return Returns a list containing:
#'
#' \item{simNetwork}{Adjacency matrix of the simulated network}
#'
#' \item{simData}{Simulated data from sim.correlation in the \code{psych} package
#' based on the simulated network}
#'
#' \item{simRho}{Simulated correlation from sim.correlation in the \code{psych} package}
#'
#' @examples
#' #Continuous data
#' sim.Norm <- sim.swn(25, 500, nei = 3)
#'
#' #Ordinal data
#' sim.Likert <- sim.swn(25, 500, nei = 3, replace = TRUE, ordinal = TRUE, ordLevels = 5)
#'
#' #Dichotomous data
#' sim.Binary <- sim.swn(25, 500, nei = 3, replace = TRUE, ordinal = TRUE, ordLevels = 2)
#'
#' @references
#' Csardi, G., & Nepusz, T. (2006).
#' The \emph{igraph} software package for complex network research.
#' \emph{InterJournal, Complex Systems}, \emph{1695}, 1-9.
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#Simulate small-world network----
sim.swn <- function (nodes, n, pos = .80, ran = c(.3,.7),
                     nei = 1, p = .5, corr = FALSE,
                     replace = NULL,
                     ordinal = FALSE, ordLevels = NULL)
{
    posdefnet <- function (net, pos, ran)
    {
        #from bootnet
        net[upper.tri(net)] <- net[upper.tri(net)] *
            sample(c(-1,1),sum(upper.tri(net)),TRUE,prob=c(pos,1-pos)) *
            runif(sum(upper.tri(net)), min(ran),max(ran))

        net[lower.tri(net)]<-net[upper.tri(net)]*net[lower.tri(net)]

        diag(net) <- rowSums(abs(net)) * 1
        diag(net) <- ifelse(diag(net)==0,1,diag(net))
        net <- net/diag(net)[row(net)]
        net <- (net + t(net)) / 2

        cornet<--cov2cor(net)
        diag(cornet)<-0

        return(cornet)
    }

    adj<-igraph::as_adj(igraph::sample_smallworld(1,nodes,nei=nei,p=p),sparse=FALSE)

    net<-adj

    graph <- posdefnet(net,pos,ran)

    #compute covariance matrix
    I<-diag(1, dim(graph)[1])
    covmat<-solve(I-graph)%*%t(solve(I-graph))

    #from bootnet
    #covmat<-solve(diag(nodes) - ngraph)

    #compute correlation matrix
    rho<-cov2cor(covmat)

    while(any(is.na(rho)))
    {
        adj<-igraph::as_adj(igraph::sample_smallworld(1,nodes,nei=nei,p=p),sparse=FALSE)

        net<-adj+t(adj)

        graph<-posdefnet(net,pos,ran)

        #compute covariance matrix
        I<-diag(1, dim(graph)[1])
        covmat<-solve(I-graph)%*%t(solve(I-graph))

        #from bootnet
        #covmat<-solve(diag(nodes) - ngraph)

        #compute correlation matrix
        rho<-cov2cor(covmat)
    }

    #generate data
    dat <- psych::sim.correlation(R=rho,n=n,data=TRUE)

    if(ordinal)
    {
        ordData <- function(data,ordLevels)
        {
            if(is.null(ordLevels))  #set number of ordinal levels
            {
                ordLevels<-5
                message("Default of 5 ordinal levels was used")
            }else{ordLevels<-ordLevels}

            for (i in 1:ncol(data))
            {data[,i] <- as.numeric(cut(data[,i],ordLevels))}

            return(data)
        }

        dat <- ordData(dat,ordLevels)
    }

    if(corr)
    {graph<-ifelse(abs(rho)>=min(ran)&abs(rho)<=max(ran),rho,0)}

    return(list(simNetwork=graph,simData=dat,simRho=rho))
}
#----