#' Repeated Responses Check
#' @description Screens data to identify potential cases of repeated responding.
#' The function is based on two criteria: no variance (i.e., a standard
#' deviation of zero for given responses)and frequency proportion of the 
#' response values (which is set by \code{freq.prop}). Note that these
#' criteria are highly related. Additional criteria will be added in 
#' the future.
#' 
#' @param data A dataset
#' 
#' @param scale.lens The number of items for each scale in the data.
#' A vector indicating the length for each scale to be checked in the data
#' 
#' @param max.val Maximum value for data (or scales).
#' If scales have different maximum values, then a vector must be
#' input with each scale's maximum value (see examples)
#' 
#' @param reverse Reverse scored responses.
#' If responses have not yet reversed, then do not reverse them.
#' If responses have been reversed, then a vector indicating
#' which responses have been reverse-scored should be input (see examples).
#' Can be TRUE/FALSE or 1/0 (reversed/not reversed)
#' 
#' @param freq.prop Frequency proportion of the response values.
#' Allows the researcher to determine the maximum frequency proportion
#' of a certain response value is suspicious.
#' The default is set to .80 (or 80 percent responses are a single value)
#' 
#' @return Returns a matrix when \code{scale.lens = NULL} and a
#' list with elements corresponding to the order of scales. In general,
#' the output contains potential bad cases that should be further
#' inspected by the researcher. 
#' 
#' @details If a case is returned, then it does not mean that it is a bad case.
#' Researchers should thoroughly inspect each case that is returned.
#' A general guideline is that if a participant responded with all middle
#' values (e.g., all 3's on a 5-point Likert scale), then they should be
#' dropped. Note that a participant who responds with all maximum or
#' minimum values may be a real case or a bad case. It is up to the
#' researcher to decide and justify why or why not a case is kept.
#' 
#' @examples
#' #Re-reverse responses
#' rev.vec <- c(TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,FALSE,
#' TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,TRUE,FALSE,TRUE,
#' FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,TRUE,FALSE,TRUE,
#' FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE)
#' 
#' #Maximum value (5-point Likert scale)
#' mv.vec <- 5
#' 
#' #Repeated responses check
#' resp.rep(neoOpen, reverse = rev.vec, max.val = mv.vec)
#' 
#' #Example with multiple scales
#' 
#' #Facet scale lengths of NEO-PI-3 Openness to Experience
#' s.len <- c(8, 8, 8, 8, 8, 8)
#' 
#' #Maximum values
#' mv.vec <- c(5, 5, 5, 5, 5, 5)
#' 
#' #Re-reverse responses
#' rev.vec <- c(TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,FALSE,
#' TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,TRUE,FALSE,TRUE,
#' FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,TRUE,FALSE,TRUE,
#' FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE)
#' 
#' #Repeated responses check
#' resp.rep(neoOpen, scale.lens = s.len, max.val = mv.vec, reverse = rev.vec)
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#Repeated Responses Check----
resp.rep <- function (data,
                       scale.lens = NULL,
                       max.val,
                       reverse = NULL,
                       freq.prop = .80)
{
    #as a matrix
    data <- as.matrix(data)
    #number of variables
    n <- ncol(data)
    #reverse items
    if(!is.null(reverse))
    {reverse <- as.matrix(reverse)}
    
    #maximum values
    if(!is.null(scale.lens))
    {
        if(length(scale.lens)!=length(max.val))
        {stop("scale.lens and max.val lengths do not match")
        }else{
        max.vec <- list()
        
        for(i in 1:length(scale.lens))
        {max.vec[[i]] <- rep(max.val[i],scale.lens[i])}
        
        max.vec <- unlist(max.vec)
        }
    }else{max.vec <- NULL}
    
    if(is.null(max.vec))
    {
        if(!is.null(max.val))
        {
            if(length(max.val)!=n)
            {
                if(length(max.val)==1)
                {max.vec <- rep(max.val,n)
                }else{stop("max.val is not a single value or match the number of variables in the data")}
            }else if(length(max.val)==n)
            {max.vec <- max.val
            }else{stop("max.val must be a single value or match the number of variables in the data")}
        }
    }
    
    #reverse code if applicable
    if(!is.null(reverse))
    {
        #convert to vector and TRUE/FALSE to numeric values
        reverse <- as.numeric(as.vector(reverse))
        #identify which variables to reverse code
        rev.vars <- which(reverse==1)
        #length of reverse variables
        rev.len <- length(rev.vars)
        
        for(i in 1:rev.len)
        {
            #maximum value for variable (plus one)
            rev.val <- max.vec[rev.vars[i]] + 1
            
            #data reversal
            data[,rev.vars[i]] <- rev.val - data[,rev.vars[i]]
        }
    }
    
    ##########################
    #POTENTIAL CHECK FUNCTION#
    ##########################
    
    check <- function (data, max.val, freq.prop)
    {
        dat.means <- rowMeans(data)
        dat.sds <- apply(data,1,sd)
    
        if(any(dat.sds==0))
        {pot.chk <- which(dat.sds==0)
        }else{pot.chk <- NULL}
    
        freq.table <- matrix(0, nrow = nrow(data), ncol = max.val)
    
        colnames(freq.table) <- c(min(data):max.val)
    
        for(i in 1:nrow(data))
        {
            vec <- data[i,]
        
            uniq <- unique(vec)
        
            uniq.len <- length(uniq)
        
            for(j in 1:uniq.len)
            {freq.table[i,paste(uniq[j])] <- length(which(vec==uniq[j]))}
        }
    
        prop.table <- freq.table/ncol(data)
    
        if(any(prop.table>=freq.prop))
        {pot.chk2 <- which(prop.table>=freq.prop,arr.ind=TRUE)[,1]
        }else{pot.chk2 <- NULL}
        
        if(!is.null(pot.chk)|!is.null(pot.chk2))
        {pot.check <- unique(c(pot.chk,pot.chk2))
        }else{pot.check <- NULL}
        
        if(!is.null(pot.check))
        {
            mat.check <- data[pot.check,]
            
            if(is.vector(mat.check))
            {
                mat.check <- as.matrix(mat.check)
                colnames(mat.check) <- pot.check
                mat.check <- t(mat.check)
            }else if(is.matrix(mat.check))
            {row.names(mat.check) <- pot.check}
        }else{mat.check <- NULL}
    
        return(mat.check)
    }
    ##########################
    
    #check means for each variable
    if(is.null(scale.lens))
    {
        return(check(data,max.val,freq.prop))
    }else{
        
        nscales <- length(scale.lens)
        
        ends <- cumsum(scale.lens)
        
        starts <- (ends - scale.lens) + 1
        
        list.check <- list()
        
        for(i in 1:nscales)
        {list.check[[as.character(i)]] <- check(data[,c(starts[i]:ends[i])],max.val[i],freq.prop)}
        
        pot <- unique(unlist(lapply(list.check,row.names)))
        
        list.check[["full"]] <- data[as.numeric(pot),]
        row.names(list.check[["full"]]) <- pot
        
        return(list.check)
    }
    
}
