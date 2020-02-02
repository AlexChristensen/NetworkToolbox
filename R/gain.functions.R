#' MFCF Gain Functions
#' 
#' @name gain.functions
#'
#' @aliases
#' gfcnv_logdet
#' gfcnv_logdet_val
#' gdcnv_lmfit
#' 
#' @description These functions maximize a gain criterion
#' for adding a node to a clique (and the larger network).
#' The flexibility of \code{\link[NetworkToolbox]{MFCF}}
#' allows for any multivariate function to be used as a
#' scoring function.
#' 
#' \itemize{
#' 
#' \item{\code{"logLik"}}
#' 
#' {The log determinant of the matrix restricted to 
#' the separator minus the log determinant of the matrix restricted
#' to the clique.}
#' 
#' \item{\code{"logLik.val"}}
#' 
#' {\code{"logLik"} with a further validation based on
#' the likelihood ratio. If the increase in gain is not significant
#' the routine stops adding nodes to the separator.}
#' 
#' \item{\code{"rSquared.val"}}
#' 
#' {The R squared from the regression of the node against the clique. Only
#' the clique nodes with a regression coefficient significantly different
#' from zero are added to the separator / new clique. The gain is different from 
#' zero only if the F-values is significant, It assumed that the \code{data} 
#' matrix is a dataset of realizations (i.e., \code{p}
#' variables and \code{n} observations).} 
#' 
#' } 
#' 
#' @usage
#' "logLik"
#' gfcnv_logdet(data, clique_id, cl, excl_nodes, ctreeControl)
#' 
#' "logLik.val"
#' gfcnv_logdet_val(data, clique_id, cl, excl_nodes, ctreeControl)
#' 
#' "rSquared.val"
#' gdcnv_lmfit(data, clique_id, cl, excl_nodes, ctreeControl)
#' 
#' @param data Matrix or data frame.
#' Can be a dataset or a correlation matrix
#' 
#' @param clique_id Numeric.
#' Number corresponding to clique to add another node to
#' 
#' @param cl List.
#' List of cliques already assembled in the network
#' 
#' @param excl_nodes Numeric vector.
#' A vector of numbers corresponding to nodes not already
#' included in the network
#' 
#' @param ctreeControl List (length = 5).
#' A list containing several parameters for controlling
#' the clique tree sizes:
#' 
#' \itemize{
#' 
#' \item{\code{min_size}}
#' {Numeric. Minimum number of nodes allowed per
#' clique. Defaults to \code{1}}
#' 
#' \item{\code{max_size}}
#' {Numeric. Maximum number of nodes allowed per
#' clique. Defaults to \code{8}}
#' 
#' \item{\code{pval}}
#' {Numeric. \emph{p}-value used to determine cut-offs for nodes
#' to include in a clique. Defaults to \code{.05}}
#' 
#' \item{\code{pen}}
#' {Numeric. Multiplies the number of edges added to penalize 
#' complex models. Similar to the penalty term in AIC}
#' 
#' \item{\code{drop_sep}}
#' {Boolean. This parameter influences the MFCF only. If TRUE any
#' separator can be used only once, as in the TMFG.}
#' 
#' \item{\code{use_returns}}
#' {Boolean. Only used in rSquared.val. If set to TRUE the regression is 
#' performed on log-returns. Defaults to \code{FALSE}}
#' 
#' }
#' 
#' @return Returns the value with the maximum gain
#' 
#' @references
#' Massara, G. P. & Aste, T. (2019).
#' Learning clique forests.
#' \emph{ArXiv}.
#' 
#' @author Guido Previde Massara <gprevide@gmail.com> and Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @importFrom stats pchisq lm.fit complete.cases pf pt
#'
#' @export
#Gain Functions----
####gfcnv_logdet####
gfcnv_logdet <- function(data, clique_id, cl, excl_nodes, ctreeControl)
{
    # small accessory function for calculation of the log of determinant
    logdet <- function(m) 
    {sum(log(pmax(svd(m)$d,0)))}
    
    # The general idea behind these gain functions is that all the nodes in a 
    # clique are sorted in order of decreasing affinity (correlation in these cases)
    # with an outstanding node. The nodes in the clique are added one by one in
    # this order until (a) the gain is zero (gfcnv_logdet) or statistically not
    # different from zero (gfcnv_logdet_val) or (b) the new clique made of the 
    # oustanding node and the nodes already added reaches ctreeControl$max_size.
    # In case ctreeControl$min_size is > 1 the gain function may add nodes even if 
    # the gain is (not significantly) greater than 0 until the new clique has
    # ctreeControl$min_size elements. The new clique is made of the oustanding node
    # and the nodes added from the previous (parent) clique. The nodes added from 
    # parent clique constitute the separator.
    # The fact that we set the order of the nodes to add in advance is not ideal as
    # it does not guarantee an optimal solution, but it provides satisfactory 
    # performance.
    
    # data is the matrix used to calculate the gains. In these gain functions
    # it is a square matrix (similarity / correlation between nodes). In other 
    # gain functions it could be a data (observations x variable) matrix.
    # clique_id is the identifier of the clique to extend.
    # cl is the vector of nodes that constitutes the clique
    # excl_nodes are the nodes not yet added to the clique forest.
    # the output is dataframe with as many rows as the size of excl_nodes
    # structure score in the dataframe is for future uses.
    retval <- data.frame(clique_id = rep(clique_id, times = length(excl_nodes)),
                         nodes = excl_nodes,
                         structure_score = 0,
                         gain = -Inf)
    
    retval$clique <- list(cl)
    
    # build a matrix with the nodes to be added in the first column and then
    # the nodes of the clique to be extended. Ideally we want to proceed from
    # left to right adding nodes from the clique to the outstanding nodes.
    x <- data[excl_nodes, cl]
    
    # due to different behaviour of some functions when the work on arrays of 
    # length one we have to specify explicitly the case when len(excl_nodes) == 1
    if(length(excl_nodes) == 1){
        # o is the oustanding node followed by the nodes in the clique in order 
        # of decreasing affinity to the oustanding_node.
        o <- c(excl_nodes, cl[order(abs(x), decreasing = TRUE)])
        # max size of final clique, either ctreeControl$max_size or parent clique + 1
        max_sz <- min(ctreeControl$max_size, length(cl)+1)
        # in this case rl is a vector that will contain the running likelihood as
        # we add nodes from the parent clique.
        rl <- matrix(data = rep(0, length(excl_nodes)*(max_sz-1)),nrow = length(excl_nodes))
        # calculate the running likelihood by increasing clique and separator
        # special case for the first case where the separator is empty
        for(sz in 2:(max_sz)){
            vs <- o[2:sz]  # running separator
            v <- o[1:sz]   # running clique
            rl[,sz-1] <- ifelse(length(vs)>1,logdet(data[vs,vs])-logdet(data[v,v]), -logdet(data[v,v]))
        }
        # apply penalty to gain
        g <- rl - ctreeControl$pen * 2:max_sz
        # in this case we have only one gain
        gains <- max(g)
        # next two lines try and get the size that gives the maximum gain or
        # at least the prescribed min clique size.
        best_sz <- which.max(g) + 1
        best_sz <- pmin(length(cl) + 1, pmax(ctreeControl$min_size, best_sz))
        # future use
        retval$structure_score <- best_sz
        # populate return value
        retval$gain <- gains
        retval$separator <- list(o[2:best_sz]) # the separator does not include the oustanding node
        retval$new_clique <- list(o[1:best_sz])
        # separator key is used by MFCF to delete records with a given separator
        # in case ctreeControl$drop_sep = TRUE
        retval$separator_key <- paste0(retval$separator, collapse = ",")
    } else {
      # This is essentially the same code as above, but it is vectorised 
      # as len(excl_nodes) > 1
      
        # order nodes in the parent clique by affinity to each oustanding node
        o <- base::t(apply(X=x, MARGIN = 1, FUN = function(x, cl) {cl[order(abs(x), decreasing = TRUE)]}  , cl = cl))
        # o1 cn11 cn12 cn13 ...
        # o2 cn21 cn22 cn23 ...
        # ...
        # ox cnx1 cnx2 cnx3 ...
        o <- cbind(excl_nodes, o)
        # same as above 
        max_sz <- min(ctreeControl$max_size, length(cl)+1)
        # running likelihood as above
        rl <- matrix(data = rep(0, length(excl_nodes)*(max_sz-1)),nrow = length(excl_nodes))
        for(sz in 2:(max_sz)){
          # apply gain function
            rl[,sz-1] <- apply(X = o[,1:sz], MARGIN = 1, 
                               FUN = function(v, M){vs=v[-1]; ifelse(length(vs)>1,logdet(M[vs,vs])-logdet(M[v,v]), -logdet(M[v,v]))}, 
                               M = data)
        }
        
        # find best number of edges consistent with the constraints
        # particular case when max_size = 2
        if (max_sz == 2){
            g <- rl-ctreeControl$pen * 1;
            retval$gain <- g
            retval$structure_score <- 2
            ol <- lapply(1:nrow(o), function(i) o[i,])
            
            retval$separator <- lapply(ol, FUN = function(v) {v[2]})
            retval$new_clique <- lapply(ol,  FUN = function(v) {v[1:2]})
            retval$separator_key <- as.character(retval$separator)
            
        } else {
            g <- base::t(apply(X = rl, MARGIN = 1, FUN = function(v) {v - ctreeControl$pen * 2:max_sz}))
            gains <- apply(X = g, MARGIN = 1, FUN = max)
            best_sz <- apply(X = g, MARGIN = 1, FUN = which.max)
            best_sz <- best_sz + 1
            idx <- which(gains <= 0)
            best_sz[idx] <- 0
            best_sz <- pmax(ctreeControl$min_size, best_sz)
            retval$structure_score <- best_sz
            retval$gain <- gains
            colnames(o) <- NULL
            
            ol <- lapply(1:nrow(o), function(i) o[i,])
            
            retval$separator <- mapply(ol, as.list(best_sz), FUN = function(v, b) {v[2:b]}, SIMPLIFY = FALSE)
            retval$new_clique <- mapply(ol, as.list(best_sz), FUN = function(v, b) {v[1:b]}, SIMPLIFY = FALSE)
            retval$separator_key <- lapply(retval$separator, FUN = function(x) paste0(sort(x), collapse=","))
        }
    }
    
    return(retval)
}
#----
####gfcnv_logdet_val####
gfcnv_logdet_val <- function(data, clique_id, cl, excl_nodes, ctreeControl){
  

  # small accessory function for calculation of the log of determinant
  logdet <- function(m) 
  {sum(log(pmax(svd(m)$d,0)))}
    
  # The general idea behind these gain functions is that all the nodes in a 
  # clique are sorted in order of decreasing affinity (correlation in these cases)
  # with an outstanding node. The nodes in the clique are added one by one in
  # this order until (a) the gain is zero (gfcnv_logdet) or statistically not
  # different from zero (gfcnv_logdet_val) or (b) the new clique made of the 
  # oustanding node and the nodes already added reaches ctreeControl$max_size.
  # In case ctreeControl$min_size is > 1 the gain function may add nodes even if 
  # the gain is (not significantly) greater than 0 until the new clique has
  # ctreeControl$min_size elements. The new clique is made of the oustanding node
  # and the nodes added from the previous (parent) clique. The nodes added from 
  # parent clique constitute the separator.
  # The fact that we set the order of the nodes to add in advance is not ideal as
  # it does not guarantee an optimal solution, but it provides satisfactory 
  # performance.
  
  # Main function
  gfcn_logdet_val <- function(X, clique_id, cl, node, ctreeControl){
      n <- ctreeControl$n
      min_size <- ctreeControl$min_size
      max_size <- ctreeControl$max_size
      clique <- cl
      candidates <- c(node, clique[order(X[node, clique]^2, decreasing = TRUE, na.last = TRUE)])
      separator <-c()
      gain = 0
      for (it in 2:min(length(candidates), ctreeControl$max_size)){
          # please see Rencher chapter 7 for this calculation. it is in essence
          # the likelihood ratio test. At every stage we build an hypothesis test
          # where H0 assumes that there is no correlation betwen the next node to be 
          # added and the nodes already added to the new clique. If we reject H0
          # we add the new node to the growing clique.
          p <- it
          m <- p - 1
          v <- candidates[1:p]
          pv <- candidates[1:(p-1)]
          s0 <- X[v,v]
          s0[1:m, p] <- 0
          s0[p, 1:m] <- 0
          j0 <- solve(s0)
          tm <- j0 %*% X[v,v]
          aqt <- n * sum(diag(tm)) - n * logdet(tm) - n * p
          qt2 <- small_sample_correction(n-1,p) * aqt
          pc2 <- pchisq(q = qt2, df = m)
          
          clique    <- v
          separator <- v[-1]
          
          if ( pc2 < 1 - ctreeControl$pval ){
              clique    <- pv
              separator <- pv[-1]
              break
          } else if (length(v) > ctreeControl$max_size){
              clique    <- v
              separator <- v[-1]
              break
          }
          gain <- logdet(X[v[-1],v[-1]]) - logdet(X[v, v])
      }
      #cat("Gain: ", gain, "\n", "Sep: ", v, "\n")
      list(gain=gain, separator=separator, clique=clique)
  }
  
  ####gfcnv_mvn_logdet_val####
  # Functions for loglikelihood ratio gain function -----
  # See Alvin Rencher - Methods of Multivariate Analysis
  # Wiley Interscience, Ch. 7 ("Tests on covariance matrices")
  # Degrees of freedom nu = n - 1
  
  # small sample correction for the statistics u (eq 7.2 Rencher)
  small_sample_correction <- function(nu,p)
  {1-1/(6*nu-1)*(2*p +1-2/(p+1))}
    
  # data is the matrix used to calculate the gains. In these gain functions
  # it is a square matrix (similarity / correlation between nodes). In other 
  # gain functions it could be a data (observations x variable) matrix.
  # clique_id is the identifier of the clique to extend.
  # cl is the vector of nodes that constitutes the clique
  # excl_nodes are the nodes not yet added to the clique forest.
  # the output is dataframe with as many rows as the size of excl_nodes
  # structure score in the dataframe is for future uses.
  retval <- data.frame(clique_id = rep(clique_id, times = length(excl_nodes)),
                       node = excl_nodes,
                       structure_score = 0,
                       gain = -Inf)
  retval$clique <- list(cl)
  
  for (it in 1:nrow(retval)) {
    # differently from the previous function we calculate the gain for a single
    # row of the return value dataframe. It is less peforming but cleaner.
    # the workhorse function is \code{gfcn_logdet_val}.
    res <- gfcn_logdet_val(data, retval[it, "clique_id"], cl, retval[it, "node"], ctreeControl)
    retval[it, "structure_score"] <- length(res$clique)
    if (length(res$separator)>0){
      retval[it, "separator"] <- list(list(res$separator))
    } else {
      retval[it, "separator"] <- list(list(NA))
    }
    retval[it, "new_clique"] <- list(list(res$clique))
    retval[it, "separator_key"] <-  paste0(sort(res$separator), collapse = ",")
    retval[it, "gain"] <- res$gain
  }
  retval
}
#----
####gdcnv_lmfit####
gdcnv_lmfit <- function(data, clique_id, cl, excl_nodes, ctreeControl){
    
    #################################
    ###### BEGIN MAIN FUNCTION ######
    #################################
    
    gdcn_lmfit <- function(X, clique_id, cl, node, ctreeControl) {
        # Bare bones linear regression using lm.fit
        # lm.fit performs much better than \code{lm} because there is no nedd to 
        # convert to / from dataframe to matrix.
        clique <- cl
        modelvars <- c(clique, node)
        
        # Take returns if required and standardise
        n <- sum(complete.cases(X[, modelvars]))
        if (n < 5) {
            # if less than 5 points then gain 0
            sep <- as.numeric()
            the_clique <- node
            gain <- 0
            
            return(list(gain = gain, #summary(the_model)$adj.r.squared, #r-squared for lm
                        clique = the_clique,
                        separator = sep))
        }
        X <- X[complete.cases(X[, modelvars]), modelvars]
        if (ctreeControl$use_returns == TRUE) {
            X <- log(X[2:n,])-log(X[1:(n-1),])
        }
        n <- n-1
        # standardise time series
        X <- scale(X)
        attr(X,"scaled:center")<-NULL 
        attr(X,"scaled:scale")<-NULL 
    
        p <- length(clique)
        cliquevars <- 1:p
        nodevar    <- p+1
        the_model <- lm.fit(as.matrix(X[,cliquevars]), X[,nodevar])
        rsq <- stats::var(the_model$fitted.values)/stats::var(X[, nodevar])
        sigma_sq <- sum(the_model$residuals^2)/the_model$df.residual
        R <- chol2inv(the_model$qr$qr[1:p, 1:p])
        se <- sqrt(diag(R)*sigma_sq)
        est <- the_model$coefficients
        tval <- est/se
        pval <- 2 * pt(abs(tval), the_model$df.residual, lower.tail = FALSE)
        # IMPORTANT. the model might change the order of the coefficient 
        # for pivoting, take that into account
        clique <- clique[the_model$qr$pivot]
        good_links <- sum(pval <= ctreeControl$pval, na.rm = TRUE)
        num_links <- max(min(good_links, ctreeControl$max_size - 1), ctreeControl$min_size -1)
        clique <- clique[order(pval)]
        if (num_links>0) {
            vars <- clique[1:num_links]
            sep <- sort(vars)
            the_clique = sort(c(node, vars))
        } else {
            sep <- as.numeric()
            the_clique <- node
        }
        
        # R-squared p-value function
        r2pval <- function(rsq, p, n){
            retval <- 1 - pf((rsq/(1-rsq))*((n-p-2)/p), df1 = p, df2 = n-p-1)
            retval
        }
        
        list(#gain = -log(r2pval(rsq, p, n)), #summary(the_model)$adj.r.squared, #r-squared for lm
            gain = rsq,
            clique = the_clique,
            separator = sep)
    }
    #################################
    ######  END MAIN FUNCTION  ######
    #################################
    
  retval <- data.frame(clique_id = rep(clique_id, times = length(excl_nodes)),
                       node = excl_nodes,
                       structure_score = 0,
                       gain = -Inf)
  retval$clique <- list(cl)
  
  for (it in 1:nrow(retval)) {
    res <- gdcn_lmfit(data, retval[it, "clique_id"], cl, retval[it, "node"], ctreeControl)
    retval[it, "structure_score"] <- length(res$clique)
    if (length(res$separator)>0){
      retval[it, "separator"] <- list(list(res$separator))
    } else {
      retval[it, "separator"] <- list(list(NA))
    }
    retval[it, "new_clique"] <- list(list(res$clique))
    retval[it, "separator_key"] <-  paste0(sort(res$separator), collapse = ",")
    retval[it, "gain"] <- res$gain
  }
  retval
}
#----