#' Permutation Test for Network Measures
#' 
#' @description Computes a permutation test to determine whether
#' there are difference in centrality and global network measures
#' 
#' @param sample1 Matrix or data frame.
#' Sample to be compared with \code{sample2}
#' 
#' @param sample2 Matrix or data frame.
#' Sample to be compared with \code{sample1}
#' 
#' @param iter Numeric.
#' Number of iterations to perform.
#' Defaults to \code{1000}
#' 
#' @param network Character.
#' Network estimation method to apply to the datasets.
#' Defaults to \code{"glasso"}
#' 
#' @param measure Character.
#' Network measure to be compared in the permutation test
#' 
#' @param alternative Character.
#' Alternative hypothesis test to perform.
#' Defaults to \code{"two.tailed"}
#' 
#' @param ncores Numeric.
#' Number of computer processing cores to use for bootstrapping samples.
#' Defaults to \emph{n} - 1 total number of cores.
#' Set to any number between 1 and maximum amount of cores on your computer
#' (see \code{parellel::detectCores()})
#' 
#' @param prev.perm \code{network.permutation} class object.
#' An object of previously performed permutation test. The
#' networks generated in the previous permutation will be
#' used to compute other network measures. This saves time
#' when computing multiple permutation tests
#' 
#' @return Returns a list containing two objects:
#' 
#' \item{result}{The results of the permutation test. For centrality measures,
#' this is a matrix where the rows represent each node and the columns are
#' the observed values of the centrality measure for \code{sample1}, \code{sample2},
#' and the \emph{p}-value from the permutation test. For global network measures,
#' this is a vector with the observed values of the global network measure for
#' \code{sample1}, \code{sample2}, and the \emph{p}-value from the permutation test.}
#' 
#' \item{networks}{A list containing two lists: \code{network1} and \code{network2}.
#' The network lists correspond to the networks generated in the permutation test
#' for \code{sample1} and \code{sample2}, respectively. This output is used primarily
#' for the computation of other network measures using the same datasets
#' (see \code{prev.perm} explanation)}
#'
#' @examples
#' # Split data (only for example)
#' split1 <- neoOpen[c(1:401),]
#' split2 <- neoOpen[c(402:802),]
#' 
#' \donttest{
#' # Perform permutation test
#' perm.str <- network.permutation(split1, split2, iter = 1000, network = "glasso",
#' measure = "strength", alternative = "two.tailed", ncores = 2)
#' 
#' # Check results
#' perm.str$result
#' 
#' # Permutation to check other measures (using networks from previous result)
#' perm.aspl <- network.permutation(prev.perm = perm.str, measure = "ASPL", ncores = 2)
#' 
#' # Check results
#' perm.aspl$result}
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Network Permutation Test----
#Updated 24.05.2021
network.permutation <- function(sample1 = NULL, sample2 = NULL, iter,
                               network = c("glasso", "ising", "TMFG", "LoGo"),
                               measure = c("betweenness", "closeness", "strength",
                                           "eigenvector", "rspbc", "hybrid", "ASPL",
                                           "CC", "S", "Q"),
                               alternative = c("less", "greater", "two.tailed"),
                               ncores, prev.perm = NULL)
{
  #### Argument check ####
  
  if(missing(iter))
  {iter <- 1000
  }else{iter <- iter}
  
  if(missing(measure))
  {stop("Argument 'measure' must be defined")
  }else{measure <- match.arg(measure)}
  
  if(is.null(prev.perm))
  {
    if(missing(network))
    {stop("Argument 'network' must be defined")
    }else{network <- match.arg(network)}
    
    if(is.null(sample1) || is.null(sample2))
    {
      stop("Arguments 'sample1' and 'sample2' must input or
           argument 'prev.perm' must have a previous result")
    }
  }else if(class(prev.perm) != "network.permutation")
  {
    stop("Object input into argument 'prev.perm' is not 
         a 'network.permutation' class object")
  }
  
  if(missing(alternative))
  {alternative <- "two.tailed"
  }else{alternative <- match.arg(alternative)}
  
  if(missing(ncores))
  {ncores <- parallel::detectCores() - 1}
  
  #### Argument check ####
  
  #### Network function ####
  
  get_network <- function(data, network)
  {
    net <- switch(
      network,
      glasso = qgraph::EBICglasso(qgraph::cor_auto(data), n = nrow(data)),
      ising = IsingFit::IsingFit(data, plot = FALSE)$weiadj,
      TMFG = NetworkToolbox::TMFG(data)$A,
      LoGo = NetworkToolbox::LoGo(data)
      )
    
    return(net)
  }
  
  #### Network function ####
  
  #### Network measure function ####
  
  get_measure <- function(net, measure)
  {
    meas <- switch(
      measure,
      betweenness = NetworkToolbox::betweenness(net),
      closeness = NetworkToolbox::closeness(net),
      strength = NetworkToolbox::strength(net),
      eigenvector = NetworkToolbox::eigenvector(net),
      rspbc = NetworkToolbox::rspbc(net),
      hybrid = NetworkToolbox::hybrid(net),
      ASPL = NetworkToolbox::pathlengths(net)$ASPL,
      CC = NetworkToolbox::clustcoeff(net)$CC,
      S = NetworkToolbox::smallworldness(net, method = "HG")$swm,
      Q = NetworkToolbox::louvain(net)$Q
    )
    
    return(meas)
  }
  
  #### Network measure function ####
  
  #### Statistic function ####
  
  get_statistic <- function(meas1, meas2) {meas1 - meas2}
  
  #### Statistic function ####
  
  # Check for previous permutation test
  if(!is.null(prev.perm))
  {
    # Obtain networks
    net.list1 <- prev.perm$networks$network1
    net.list2 <- prev.perm$networks$network2
    
    # Obtain names of samples
    sample1.name <- colnames(prev.perm$result)[1]
    sample2.name <- colnames(prev.perm$result)[2]
    
  }else{
    
    # Obtain names of samples
    sample1.name <- paste(deparse(substitute(sample1)))
    sample2.name <- paste(deparse(substitute(sample2)))
    
    # Combine samples
    comb.sample <- rbind(sample1, sample2)
    
    # Initialize data permutation list
    data.list1 <- vector("list", length = iter)
    data.list2 <- data.list1
    
    data.list1[[1]] <- sample1
    data.list2[[1]] <- sample2
    
    # Message to user
    message("Generating permutated samples...", appendLF = FALSE)
    
    # Generate permutated samples
    for(i in 2:iter)
    {
      # Randomly draw sample 1
      k <- sample(1:nrow(comb.sample), nrow(sample1))
      
      # New sample 1
      data.list1[[i]] <- comb.sample[k,]
      
      # New sample 2
      data.list2[[i]] <- comb.sample[-k,]
    }
    
    # Message to user
    message("done")
    
    #### Estimate networks ####
    message("Estimating networks...\n", appendLF = FALSE)
    
    # Parallel processing
    cl <- parallel::makeCluster(ncores)
    
    # Export data lists
    parallel::clusterExport(cl = cl,
                            varlist = c("data.list1", "data.list2",
                                        "get_network", "network"),
                            envir = environment())
    
    
    # Compute networks
    net.list1 <- pbapply::pblapply(data.list1, cl = cl,
                                   FUN = get_network,
                                   network = network)
    
    # Compute networks
    net.list2 <- pbapply::pblapply(data.list2, cl = cl,
                                   FUN = get_network,
                                   network = network)
    
    # Stop cluster
    parallel::stopCluster(cl)
    
    #### Estimate networks ####
    
  }
  
  #### Compute measures ####
  message("Computing measures...\n", appendLF = FALSE)
  
  # Parallel processing
  cl <- parallel::makeCluster(ncores)
  
  # Export network lists
  parallel::clusterExport(cl = cl,
                          varlist = c("net.list1", "net.list2",
                                      "get_measure", "measure"),
                          envir = environment())
  
  
  #Compute measures
  stat.list1 <- pbapply::pbsapply(net.list1, cl = cl,
                                 FUN = get_measure,
                                 measure = measure)
  
  #Compute measures
  stat.list2 <- pbapply::pbsapply(net.list2, cl = cl,
                                 FUN = get_measure,
                                 measure = measure)
  
  #Stop cluster
  parallel::stopCluster(cl)
  
  #### Compute measures ####
  
  # Compute statistic
  if(measure %in% c("betweenness", "closeness", "strength",
                    "eigenvector", "rspbc", "hybrid"))
  {
    diff0 <- stat.list1[,1] - stat.list2[,1]
    diff <- stat.list1 - stat.list2
    
    # Alternative test
    if(alternative == "greater")
    {p <- rowMeans(diff >= diff0, na.rm = TRUE)
    }else if(alternative == "less")
    {p <- rowMeans(diff <= diff0, na.rm = TRUE)
    }else{p <- rowMeans(abs(diff) >= abs(diff0), na.rm = TRUE)}
    
    res <- round(cbind(stat.list1[,1], stat.list2[,1], p),4)
    
    colnames(res) <- c(sample1.name, sample2.name, "p-value")
    
  }else{
    diff0 <- stat.list1[1] - stat.list2[1]
    diff <- stat.list1 - stat.list2
    
    # Alternative test
    if(alternative == "greater")
    {p <- mean(diff >= diff0, na.rm = TRUE)
    }else if(alternative == "less")
    {p <- mean(diff <= diff0, na.rm = TRUE)
    }else{p <- mean(abs(diff) >= abs(diff0), na.rm = TRUE)}
    
    res <- round(c(stat.list1[1], stat.list2[1], p),4)
    
    names(res) <- c(sample1.name, sample2.name, "p-value")
  }
  
  # Result list
  res.list <- list()
  res.list$result <- res
  res.list$networks <- list(network1 = net.list1, network2 = net.list2)
  class(res.list) <- "network.permutation"
  
  return(res.list)
}
#----