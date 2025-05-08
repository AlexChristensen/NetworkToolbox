#' Simulate Chordal Network
#' @description Simulates a chordal network based on number of nodes.
#' Data will also be simulated based on the true network structure
#'
#' @param nodes Numeric.
#' Number of nodes in the simulated network
#'
#' @param inverse Character.
#' Method to produce inverse covariance matrix.
#'
#' \itemize{
#'
#' \item \code{"cases"} --- Estimates inverse covariance matrix
#' based on \code{n} number of cases and \code{nodes} number of
#' variables, which are drawn from a random normal distribution
#' \code{rnorm}. Data generated will be continuous unless
#' \code{ordinal} is set to \code{TRUE}
#'
#' \item \code{"matrix"} --- Estimates inverse covariance matrix
#' based on sigma
#'
#' }
#'
#' @param n Numeric.
#' Number of cases in the simulated dataset
#'
#' @param ordinal Boolean.
#' Should simulated continuous data be converted to ordinal?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} for simulated ordinal data
#'
#' @param ordLevels Numeric.
#' If \code{ordinal = TRUE}, then how many levels should be used?
#' Defaults to \code{5}.
#' Set to desired number of intervals
#'
#' @param idio Numeric.
#' DESCRIPTION.
#' Defaults to \code{0.10}
#'
#' @param eps Numeric.
#' DESCRIPTION.
#' Defaults to \code{2}
#'
#' @return Returns a list containing:
#'
#' \item{cliques}{The cliques in the network}
#'
#' \item{separators}{The separators in the network}
#'
#' \item{inverse}{Simulated inverse covariance matrix of the network}
#'
#' \item{data}{Simulated data from sim.correlation in the \code{psych}
#' package based on the simulated network}
#'
#' @examples
#' #Continuous data
#' sim.Norm <- sim.chordal(nodes = 20, inverse = "cases", n = 1000)
#'
#' #Ordinal data
#' sim.Likert <- sim.chordal(nodes = 20, inverse = "cases", n = 1000, ordinal = TRUE)
#'
#' #Dichotomous data
#' sim.Binary <- sim.chordal(nodes = 20, inverse = "cases", n = 1000, ordinal = TRUE, ordLevels = 5)
#'
#' @references
#' Massara, G. P. & Aste, T. (2019).
#' Learning clique forests.
#' \emph{ArXiv}.
#'
#' @author Guido Previde Massara <gprevide@gmail.com>
#'
#' @importFrom stats rnorm
#'
#' @export
#Simulate random chordal network
sim.chordal <- function(nodes, inverse = c("cases","matrix"),
                        n = NULL, ordinal = FALSE, ordLevels = NULL,
                        idio = NULL, eps = NULL)
{
    # initialize cliques and separators list
    clq <- list()
    sep <- list()

    # initialize number of cliques and separators
    num_clq <- 0
    num_sep <- 0

    # initiliaze nodes
    vertices <- 1:nodes

    # first vertex is random
    v <- sample(x = vertices, size = 1)
    # remove from oustanding vertices
    vertices <- setdiff(vertices, v)
    # create first clique with first vetrex
    num_clq <- num_clq + 1
    clq[[num_clq]] <- v

    while (length(vertices) > 0)
    {
        # pick a random vertex
        if (length(vertices) == 1)
        {v <- vertices
        }else
        {v <- sample(x = vertices, size = 1)}
        # (and remove from outstanding)
        vertices <- setdiff(vertices, v)

        # pick a clique, if zero means create a new clique
        clique_to_expand <- sample(x = 0:num_clq, size = 1)

        if (clique_to_expand == 0)
        {
            # new clique
            num_clq <- num_clq + 1
            clq[[num_clq]] <- v
            # no separator
        }else
        {
            clq_len <- length(clq[[clique_to_expand]])
            sep_size <- sample(x=1:clq_len, size = 1)

            if(length(clq[[clique_to_expand]]) == 1)
            {tsep <- clq[[clique_to_expand]]
            }else
            {tsep <- sample(x = clq[[clique_to_expand]], size = sep_size)}

            if(sep_size == clq_len)
            {
                # expand existing clique
                clq[[clique_to_expand]] <- c(clq[[clique_to_expand]], v)
                # no separator
            }else
            {
                # new clique
                num_clq <- num_clq + 1
                clq[[num_clq]] <- c(tsep, v)
                # add separator  in this case
                num_sep <- num_sep + 1
                sep[[num_sep]] <- tsep
            }
        }
    }

    res <- list(cliques = clq, separators = sep)

    if(inverse == "cases")
    {
        # missing arguments
        if(is.null(n))
        {
            n <- 1000
            message("Default sample size of 1000 was used")
        }else{n <- n}

        if(is.null(idio))
        {
            idio <- 0.10
            message("Default idio of 0.10 was used")
        }else{idio <- idio}

        if(is.null(ordLevels))
        {
            ordLevels <- 5
            message("Default of 5 ordinal levels was used")
        }else{ordLevels <- ordLevels}

        # initialize cases and inverse correlation matrix
        cases <- matrix(0, nrow = n, ncol = nodes)
        J <- matrix(0, ncol = nodes, nrow = nodes)

        # generate data for cases
        for (cl in clq)
        {cases[, cl] <- cases[, cl] + (1+rnorm(n)) + idio * matrix(data=rnorm(length(cl)*n), nrow = n, ncol = length(cl))}

        # generate inverse correlations for cliques
        for (cl in clq)
        {
            j <- solve(cor(cases[, cl], cases[, cl]))
            J[cl, cl] <- J[cl, cl] + j
        }

        # substract inverse correlations from cliques given separators
        for (sp in sep)
        {
            j <- solve(cor(cases[, sp], cases[, sp]))
            J[sp, sp] <- J[sp, sp] - j
        }

        if(ordinal)
        {
            ordData <- function(data,ordLevels)
            {
                if(is.null(ordLevels))  #set number of ordinal levels
                {ordLevels <- 5
                }else{ordLevels <- ordLevels}

                for (i in 1:ncol(data))
                {data[,i] <- as.numeric(cut(data[,i], ordLevels))}

                return(data)
            }

            cases <- ordData(cases, ordLevels)
        }

        res$inverse <- J
        res$data <- cases
    }else if(inverse == "matrix")
    {
        # missing arguments
        if(is.null(eps))
        {eps <- 2
        }else{eps <- eps}

        # initialize inverse correlation matrix
        J <- matrix(0, ncol = nodes, nrow = nodes)

        # generate inverse correlations for cliques
        for (cl in clq)
        {
            j <- solve((1-eps) * diag(length(cl)) + matrix(eps, ncol = length(cl), nrow = length(cl)))
            J[cl, cl] <- J[cl, cl] + j
        }

        # substract inverse correlations from cliques given separators
        for (sp in sep)
        {
            j <- solve((1-eps) * diag(length(sp)) + matrix(eps, ncol = length(sp), nrow = length(sp)))
            J[sp, sp] <- J[sp, sp] - j
        }

        res$inverse <- J
    }

    return(res)

}
#----