#' Import CONN Toolbox Brain Matrices to R format
#' 
#' @description Converts a Matlab brain z-score connectivity array (n x n x m)
#' where \strong{n} is the n x n connectivity matrices and \strong{m} is the participant.
#' If you would like to simply import a connectivity array from Matlab, then see the examples
#' 
#' @param MatlabData Input for Matlab data file.
#' Defaults to interactive file choice
#' 
#' @param progBar Should progress bar be displayed?
#' Defaults to \code{TRUE}.
#' Set \code{FALSE} for no progress bar
#' 
#' @return Returns a list containing:
#' 
#' \item{rmat}{Correlation matrices for each participant (m) in an array (n x n x m)}
#' 
#' \item{zmat}{Z-score matrices for each participant (m) in an array (n x n x m)}
#' 
#' @examples
#' \dontrun{
#' neuralarray <- convertConnBrainMat()
#' 
#' #Import correlation connectivity array from Matlab
#' library(R.matlab)
#' neuralarray<-readMat(file.choose())
#' }
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
#Convert CONN Toolbox Brain Matrices----
#Updated 07.03.2020
convertConnBrainMat <- function (MatlabData, progBar = TRUE)
{
    ###########################
    #### MISSING ARGUMENTS ####
    ###########################
    
    if(missing(MatlabData))
    {mat <- R.matlab::readMat(file.choose())
    }else{mat <- R.matlab::readMat(MatlabData)}
    
    #######################
    #### MAIN FUNCTION ####
    #######################
    
    # If data imported is not in a list, then return
    if(!is.list(mat))
    {return(mat)
    }else{ # If data is in a list
        
        # Determine structure of data
        if("data" %in% names(mat)) # Time series data
        {
            # Grab data
            dat <- unlist(mat$data, recursive = FALSE)
            
            # Grab names
            dat.names <- unlist(mat$names)
            
            # Get names of ROIs
            names(dat) <- dat.names
            
            # Check for CSF, grey and white matter
            if(any(c("Grey Matter", "White Matter", "CSF") %in% names(dat)))
            {dat <- dat[-which(names(dat) == c("Grey Matter", "White Matter", "CSF"))]}
            
            return(dat)
            
        }else if("Z" %in% names(mat)) # Average time series data
        {
            #read in matlab data
            n1<-nrow(mat$Z) #determine number of rows
            n2<-ncol(mat$Z) #determine number of columns
            if(nrow(mat$Z)!=ncol(mat$Z))
            {warning("Row length does not match column length")}
            m<-length(mat$Z)/n1/n2 #determine number of participants
            
            #change row and column names
            coln1<-matrix(0,nrow=n1) #get row names
            for(i in 1:n1)
            {coln1[i,]<-mat$names[[i]][[1]][1,1]}
            
            coln2<-matrix(0,nrow=n2) #get column names
            for(i in 1:n2)
            {coln2[i,]<-mat$names2[[i]][[1]][1,1]}
            
            dat<-mat$Z
            if(progBar)
            {pb <- txtProgressBar(max=m, style = 3)}
            
            for(i in 1:m) #populate array
            {
                dat[,,i]<-psych::fisherz2r(mat$Z[,,i])
                for(j in 1:n1)
                    for(k in 1:n2)
                        if(is.na(dat[j,k,i]))
                        {dat[j,k,i]<-0}
                if(progBar){setTxtProgressBar(pb, i)}
            }
            if(progBar){close(pb)}
            
            colnames(dat)<-coln2
            row.names(dat)<-coln1
            
            return(list(rmat=dat,zmat=mat$Z))
        }
        
    }
}
#----