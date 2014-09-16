#*/#########################################################################
#  Cache Matrix - Functions to support a matrix that caches its inverse
#
#  Author: Myat Khine Oo
#  Created: 16/09/2014
#  Supported by: Myat Khine Oo
#
#
#  Functions:
#    makeCacheMatrix - Creates a special matrix "object" that caches its inverse
#    cacheSolve      - Wrapper for solve to run on special matrix object to use its cache functionality
#
#*/#########################################################################



#*/#########################################################################
# Creates a special matrix "object" that caches its inverse.
# Matricies are maintained in a list of matrices and the inverses.
# 
# Args:
#     x:          Initial matrix to initialise object
#     size:       Number of matricies/inverse pairs currently maintained in the lists
#     autoExtend: Flag to extend size of lists if the lists are full
# 
# Returns:
#     List of functions to manipulate the matrix object.
#*/#########################################################################
makeCacheMatrix <- function(x = matrix(), size = 1, autoExtend = T) {
    # Exception handling if size argument is not positive number
    if (size < 1) {
        message("Size set to 0 or negative number. Setting to 1")
        size <- 1
    }
    
    # Maintained list to manage matrix/inverse pairs
    Matrix <- as.list(replicate(size, NULL))
    Inverse <- as.list(replicate(size, NULL))
    
    # How much of the lists are full
    currentFill <- 1
    
    # Indext to current matrix being worked on.
    # First value indicates value is in Matrix, else in Inverse list
    # Only ONE value can be not 0
    currentValue <- c(1,0)
    
    # set the initial matrix
    Matrix[[1]] <- x
    
    #*/#########################################################################
    # INTERNAL FUNCTION: Determines if two matricies are equal
    # Reference: https://stat.ethz.ch/pipermail/r-help/2012-June/315408.html
    #
    # Args:
    #     x: A matrix
    #     y: A matrix
    #
    # Returns:
    #    True if equal else false
    #*/#########################################################################
    matequal <- function(x, y) {
        is.matrix(x) && is.matrix(y) && length(dim(x)) == length(dim(y)) &&
        all(dim(x) == dim(y)) && all(x == y) 
    }
    
    #*/#########################################################################
    # INTERNAL FUNCTION: Search for a given matrix in the Matrix and Inverse 
    # lists and return a position if found.
    #
    # Args:
    #     y: A matrix to be searched.
    #
    # Returns:
    #    Returns position similar to currentValue format. If not found returns
    #    c(0,0)
    #*/#########################################################################
    search <- function(y) {
        searchMatrix <- which((lapply(Matrix,matequal,y)) == T)
        searchInverse <- which((lapply(Inverse,matequal,y)) == T)
        
        ans <- c(0,0)
        if (length(searchMatrix) == 0 && length(searchMatrix) == 0){
            ans <- c(0,0)
        } else if (length(searchMatrix) > 0) {
            ans <- c(searchMatrix[1],0)  
        } else {  
            ans <- c(0, searchInverse[1]) 
        }
        
        return(ans)
    }
    
    #*/#########################################################################
    # Extend the size of the Matrix and Inverse lists
    #
    # Args:
    #     rowsize: Number of rows to be extended by
    #
    # Returns:
    #    Returns the new size of the matrix
    #*/#########################################################################
    extend <- function (rowsize = NULL) {
        # Exception Handling
        if (is.null(rowsize)) { 
            stop("Row size must be set. Exit without extending") 
        } else if (rowsize < 1) { 
            stop("Invalid row size. Exit without extending") 
        }
        
        # Extend size of lists
        size <<- size + rowsize 
        length(Matrix) <<- size
        length(Inverse) <<- size
        
        return (size)
    }

    #*/#########################################################################
    # Set auto extend functionality for the size of the Matrix and Inverse lists
    #
    # Args:
    #     extend: Boolean value.
    #
    # Returns:
    #    Returns the new setting of autoExtend.
    #*/#########################################################################
    setAutoExtend <- function (extend = NULL) {
        # Exception Handling
        if (is.null(extend)){ 
            stop("Extend parameter not set. Function not doing anything.") 
        } else if (is.na(as.logical(extend))) {
            stop("Invalid extend parameter.")  
        }
        
        autoExtend <<- as.logical(extend)
        return(autoExtend)
    }
    
    #*/#########################################################################
    # Get the current value based on the currentValue index
    #
    # Returns:
    #    Returns the currently set matrix.
    #*/#########################################################################
    get <- function() {    
        if (currentValue[1] > 0) {
            ans <- Matrix[[currentValue[1]]]
        } else {
            ans <- Inverse[[currentValue[2]]]
        }
        
        return(ans)
    }
    
    #*/#########################################################################
    # Set a new matrix as the current Matrix.
    # If the matrix has been set in the past or an inverse has been calculated
    # then set back to that.
    #
    # Args:
    #     y: New matrix to be set as current.
    #*/#########################################################################
    set <- function (y) {    
        # Search if matrix has been encountered in the past
        index <- search(y)
        
        # If not found
        if (sum(index) == 0) {
            if (currentFill == size && autoExtend == T) {
                # If lists are full and auto Extend has been set
                # Extend lists and move index down
                extend(size)
                currentFill <<- currentFill + 1
            } else if(currentFill == size && size == 1) {
                # If lists are full and only one element in the lists
                # Will overwrite that value later
                # Clear out the Inverse list
                currentFill <<- currentFill
                Inverse [[currentFill]] <<- as.list(replicate(size, NULL))
            } else if (currentFill == size) {
                # If full but list size is more than 1
                # Pop off the oldest entry in list
                Matrix <<- Matrix[2:size]
                Inverse <<- Inverse[2:size]
                # Leave space for new value
                length(Matrix) <<- size 
                length(Inverse) <<- size        
            } else { 
                # Lists are not yet full so just go to next value
                currentFill <<- currentFill + 1 
            }
            
            # Add new matrix to the lists
            Matrix [[currentFill]] <<- y
            # Reset the index
            currentValue <<- c(currentFill,0)
        } else { # If found 
            # Set to current location.
            currentValue <<- index
        }
    }
    
    #*/#########################################################################
    # Set a value to be the inverse of the current matrix.
    # If the current matrix is not "square" matrix or the current matrix matrix
    # multiplied (%*%) is not the identity matrix then error.
    # Otherwise, cache the value.
    #
    # Args:
    #     inverse: Inverse of current matrix
    #*/#########################################################################
    setinverse <- function (inverse) {
        # Current matrix
        matrix <- get()
        
        # Exception handling
        # Not a square matrix
        if(nrow(matrix) != ncol(matrix)) 
            stop("Current matrix cannot have inverse")
        
        # Does not produce identity matrix
        identity <- diag(nrow(matrix))
        if(!(matequal(matrix %*% inverse, identity)))
            stop("Input parameter is not inverse of matrix")
        
        # Cache the value
        if (currentValue[1] > 0) { 
            Inverse[[sum(currentValue)]] <<- inverse 
        } else { 
            Matrix[[sum(currentValue)]] <<- inverse 
        }
    }
    
    #*/#########################################################################
    # Return the value of the inverse of the current matrix.
    # If not cached, then solves the inverse and caches the matrix.
    # If solving does not produce the inverse then the value is returned
    # without caching the value.
    #
    # Returns:
    #     Value of solve(x) where x is the current matrix.
    #*/#########################################################################
    getinverse <- function () {
        if (currentValue[1] > 0) { 
            ans <- Inverse[[sum(currentValue)]]  
        } else { 
            ans <- Matrix[[sum(currentValue)]] 
        }
        
        if (is.null(ans)) { 
            ans <- solve(get())
            tryCatch(setinverse(ans), 
                     error = function(e) {
                         warning("Current matrix does not have an inverse. Result is not cached")
                     }
            )
        } else {
            message("Retrieved value from cache")
        }
        
        return(ans)
    }
    
    #*/#########################################################################
    # Returns current fill statistics of the lists.
    #
    # Returns:
    #     A list including currentFill and TotalSize.
    #*/#########################################################################
    getFillStatistics <- function() {
        list(CurrentFill = currentFill, TotalSize = size)
    }
    
    # Functions available to interact with this object
    list(set               = set, 
         get               = get, 
         setinverse        = setinverse, 
         getinverse        = getinverse, 
         extend            = extend,   
         setAutoExtend     = setAutoExtend,
         getFillStatistics = getFillStatistics
    )
}


#*/#########################################################################
# Wrapper for solve() to work with special matrix "object" created by 
# makeCacheMatrix (). 
# Will use the getInverse() function from makeCacheMatrix() if the solve is used
# to calculate the inverse. Otherwise, it will run solve() as required.
#
# 
# Args:
#     x:   Special matrix object that caches its inverse. Must be square matrix
#     b:   Optional matrix to solve y for x %*% y == b
#     ...: Other parameters to be passed to solve
# 
# Returns:
#     Result of the solve() function applied to x
#*/#########################################################################
cacheSolve <- function(x, b, ...) {    
    mat <- x$get()
    
    # Exception handling
    # x must be a square matrix as per solve() documentation
    if (nrow(mat) != ncol(mat))
        stop("'x' (", nrow(mat), " x ", ncol(mat), ") must be square")
    
    tstB <- F
    # If b has been set, check if be is the identity matrix of x
    if (!missing(b)) {
        tstID <- diag(nrow(mat))
        tstB <- is.matrix(b) && dim(b) == dim(tstID) && all(b == tstID)
    }
    
    # If b has not been set or is the identity matrix of x then get inverse.
    # Otherwise, run the solve function
    if (missing(b) || tstB) { 
        ans <- x$getinverse() 
    } else { 
        ans <- solve(mat, b, ...) 
    }
    
    return(ans)
}
