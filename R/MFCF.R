#' Maximally Filtered Clique Forest
#' 
#' @description Applies the Maximally Filtered Clique Forest (MFCF) filtering method
#' (\strong{Please see and cite Massara & Aste}).
#' 
#' @param data Matrix (n \code{x} n or p \code{x} n) or data frame.
#' Can be a dataset or a correlation matrix
#' 
#' @param cases Numeric. If \code{data} is a (partial) correlation
#' matrix, then number of cases must be input.
#' Defaults to \code{NULL}
#' 
#' @param na.data Character.
#' How should missing data be handled?
#' 
#' \itemize{
#' 
#' \item{\code{"listwise"}}
#' {Removes case if \strong{any} missing data exists.
#' Applies \code{\link{na.omit}}}
#' 
#' \item{\code{"pairwise"}}
#' {Estimates correlations using the available data
#' for each variable}
#' 
#' \item{\code{"fiml"}}
#' {Estimates correlations using the Full Information
#' Maximum Likelihood. Recommended and most robust but time consuming}
#' 
#' \item{\code{"none"}}
#' {Default. No missing data or missing data has been
#' handled by the user}
#' 
#' }
#' 
#' @param time.series Boolean.
#' Is \code{data} a time-series dataset?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} to handle time-series data (n \code{x} p)
#' 
#' @param gain.fxn Character.
#' Gain function to be used for inclusion of nodes in cliques.
#' There are several options available
#' (see \code{\link[NetworkToolbox]{gain.functions}} for more details):
#' \code{"logLik"}, \code{"logLik.val"}, \code{"rSquared.val"}.
#' Defaults to \code{"rSquared.val"}
#' 
#' @param min_size Numeric. Minimum number of nodes allowed per
#' clique. Defaults to \code{0}
#' 
#' @param max_size Numeric. Maximum number of nodes allowed per
#' clique. Defaults to \code{8}
#' 
#' @param pval Numeric. \emph{p}-value used to determine cut-offs for nodes
#' to include in a clique
#' 
#' @param pen Numeric. Multiplies the number of edges added to penalise 
#' complex models. Similar to the penalty term in AIC
#' 
#' @param drop_sep Boolean. This parameter influences the MFCF only.
#' Defaults to \code{FALSE}.
#' If \code{TRUE}, then any separator can be used only once (similar
#' to the \code{\link[NetworkToolbox]{TMFG}})
#' 
#' @param use_returns Boolean. Only used in \code{"gain.fxn = rSquared.val"}.
#' If set to \code{TRUE} the regression is 
#' performed on log-returns.
#' Defaults to \code{FALSE}
#' 
#' @return Returns a list containing:
#' 
#' \item{A}{MFCF filtered partial correlation network (adjacency matrix)}
#' 
#' \item{J}{MFCF filtered inverse covariance matrix (precision matrix)}
#' 
#' \item{cliques}{Cliques in the network
#' (output for \code{\link[NetworkToolbox]{LoGo}})}
#' 
#' \item{separators}{Separators in the network
#' (output for \code{\link[NetworkToolbox]{LoGo}})}
#' 
#' @examples
#' # Load data
#' data <- neoOpen
#' 
#' \dontrun{ 
#' # Use polychoric correlations and R-squared method
#' MFCF.net <- MFCF(qgraph::cor_auto(data), cases = nrow(neoOpen))$A
#' 
#' }
#' 
#' @references
#' Massara, G. P. & Aste, T. (2019).
#' Learning clique forests.
#' \href{https://arxiv.org/abs/1905.02266}{https://arxiv.org/abs/1905.02266}
#' 
#' @author Guido Previde Massara <gprevide@gmail.com> and Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @importFrom utils combn
#' 
#' @export
#MFCF Filtering Method----
MFCF <- function(data, cases = NULL,
                 na.data = c("pairwise","listwise","fiml","none"),
                 time.series = FALSE,
                 gain.fxn = c("logLik","logLik.val","rSquared.val"),
                 min_size = 0,    
                 max_size = 8,    
                 pval = .05,      
                 pen = 0.0,       
                 drop_sep = FALSE,
                 use_returns = FALSE
                 )
{
    #gain function
    if(missing(gain.fxn))
    {gain.fun <- "rSquared.val"
    }else{gain.fun <- match.arg(gain.fxn)}
    
    if(!time.series)
    {
        #missing data handling
        if(nrow(data)!=ncol(data))
        {
            #number of cases (n)
            cases <- nrow(data)
            
            if(any(is.na(data)))
            {
                if(missing(na.data))
                {stop("Missing values were detected! Set 'na.data' argument")
                }else{data <- qgraph::cor_auto(data, missing = na.data)}
            }else{data <- qgraph::cor_auto(data)}
        }else{
            #check that 'n' is given
            if(is.null(cases))
            {stop("Number of cases (cases) must be input")}
        }
    }else{cases <- nrow(data)}
    
    #number of cases (n) variables (p)
    n <- nrow(data)
    p <- ncol(data)
  
    #initialize the clique-tree control structure from arguments, keeps everything packed together
    ctreeControl = list(n = cases,
                        min_size = min_size,    
                        max_size = max_size,     
                        pval = pval,       
                        pen = pen,      
                        drop_sep = drop_sep,
                        use_returns = use_returns
                        )
    
    #nodes yet to be included (excluded nodes)
    outstanding_nodes <- 1:p
    
    #initialize lists for cliques and separators
    the_cliques <- list()
    the_separators <- list()
    
    #initialize number of cliques and separators
    clique_no <- 0
    sep_no <- 0
    
    #initialize gain function
    if(gain.fun == "logLik")
    {gain.fun <- match.fun(gfcnv_logdet)
    }else if(gain.fun == "logLik.val")
    {gain.fun <- match.fun(gfcnv_logdet_val)
    }else if(gain.fun == "rSquared.val")
    {gain.fun <- match.fun(gdcnv_lmfit)}
    
    #insert first clique
    if(n == p)
    {sums <- rowSums(data * (data > rowMeans(data)))
    }else
    {
        #make kernel function
        make_kernel <- function(X, gain.fun, ct_control)
        {
          gf <- function(i,j) gain.fun(X, 1, i, j, ct_control)$gain
          ct_control$n <- ctreeControl$n
          ct_control$min_size <- 2
          ct_control$max_size <- 2
          ct_control$pval <- 1
          ct_control$pen <- ctreeControl$pen
          ct_control$drop_sep <- ctreeControl$drop_sep
          ct_control$use_returns <- ctreeControl$use_returns
          p <- ncol(X)
          C <- combn(p,2)
          family <- family
          gains <- mapply(gf, C[1,], C[2,])
          K <- matrix(0, nrow = p, ncol = p)
          K[1:(p-1), ] <- as.matrix(Matrix::sparseMatrix(i = C[1,], j = C[2,], x = gains))
          K <- K + t(K) 
          diag(K) <- 1
        }
        
        K <- make_kernel(data,gain.fun,ctreeControl)
        
        sums <- rowSums(K * (K > rowMeans(K)))
        
        K <- NULL
    }
    
    #updated cliques with first clique
    first_clique <- order(sums, decreasing = TRUE)[1:max(2,ctreeControl$min_size)]
    clique_no <- 1
    the_cliques[[clique_no]] <- first_clique
    outstanding_nodes <- setdiff(outstanding_nodes, the_cliques[[clique_no]])
    
    #initialize gain table
    gain_table <- gain.fun(data,clique_no,the_cliques[[clique_no]],outstanding_nodes,ctreeControl)
    
    #MCMF algorithm
    while(length(outstanding_nodes)>0){
        # get most favourable entry in gain table
        # max_rec <- gain_table[gain_table$structure_score == max(gain_table$structure_score),]
        # max_rec <- max_rec[max_rec$gain == max(max_rec$gain),]
        idx <- which.max(gain_table$gain)
        max_rec <- gain_table[idx[1],] # keep only the first match
        if (max_rec$gain <= 0) {
            new_clique <- max_rec$node
            parent_clique <- NA
            parent_clique_id <- NA
            the_node <- max_rec$node
            the_sep <- list(NA)
            clique_no <- clique_no + 1
            clique_to_update <- clique_no
            #new separator
            sep_no <- sep_no + 1
            the_separators[[sep_no]] <- the_sep[[1]]
            the_cliques[[clique_to_update]] <- new_clique
            outstanding_nodes <- outstanding_nodes[outstanding_nodes != the_node]
            gain_table <- gain_table[gain_table$clique_id != clique_no & gain_table$node != the_node,]
            #cat(length(outstanding_nodes), "\n")
            next
        } else {
            new_clique <- max_rec$new_clique[[1]]
            parent_clique <- max_rec$clique[[1]]
            parent_clique_id <- max_rec$clique_id
            the_node <- max_rec$node
            the_sep <- max_rec$separator
        }
        outstanding_nodes <- outstanding_nodes[outstanding_nodes != the_node]
        
        #case new clique
        if (length(new_clique) <= length(parent_clique)){
            clique_no <- clique_no + 1
            clique_to_update <- clique_no
            #new separator
            sep_no <- sep_no + 1
            the_separators[[sep_no]] <- the_sep[[1]]
            the_cliques[[clique_to_update]] <- new_clique
        } else { # extend old clique
            clique_to_update <- parent_clique_id
            the_cliques[[clique_to_update]] <- new_clique
        }
        # don't update at the last round
        if (length(outstanding_nodes) == 0) {
            break
        }
        gain_table <- gain_table[gain_table$clique_id != clique_no & gain_table$node != the_node,]
        #gain_table <- gain_table[gain_table$node != the_node,]
        if(ctreeControl$drop_sep == TRUE) {
            the_sep_key <- paste0(sort(the_sep[[1]]), collapse = ",")
            gain_table <- gain_table[gain_table$separator_key != the_sep_key,]
        }
        gain_table <- rbind(gain_table, 
                            gain.fun(data, clique_to_update, new_clique, outstanding_nodes, ctreeControl))
        # df <- gain_func_p(clique_to_update, new_clique, outstanding_nodes)
        # gain_table <- rbind(gain_table, df)
        #cat(length(outstanding_nodes), "\n")
    }
    
    
    # remove cliques with no nodes
    if(any(is.na(the_cliques)))
    {the_cliques[[which(is.na(the_cliques))]] <- NULL}
    
    # remove separators with no nodes
    if(any(is.na(the_separators)))
    {the_separators[[which(is.na(the_separators))]] <- NULL}
    
    
    #results list
    res <- list()
    res$cliques <- the_cliques
    res$separators <- the_separators
    res$A <- LoGo(data, cliques = the_cliques, separators = the_separators, partial = TRUE)
    res$J <- LoGo(data, cliques = the_cliques, separators = the_separators, partial = FALSE)
    
    return(res)
}
