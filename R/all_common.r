
##--- All the general, common functions that don't belong to any one specific class in the MTVGARCH Package ---##

library(Matrix)

## ---  Enums:

## TV
tvshape <- list(delta0only=0,single=1,double=2,double1loc=3)
speedopt <- list(none=0,gamma=1,gamma_std=2,eta=3,lamda2_inv=4)
## GARCH
garchtype <- list(noGarch=0,general=1,gjr=2)
## COR
corrtype <- list(CCC=1,STCC1=2,STCC2=3,STEC1=4,STEC2=5)
corrshape <- list(single=1,double=2,double1loc=3)
corrspeedopt <- list(gamma=1,gamma_std=2,eta=3)

## -- vecl -- ####
setGeneric(name="vecL",
           valueClass = "numeric",
           signature = c("sqrMatrix"),
           def = function(sqrMatrix){
             ## Returns the lower triangle of a square matrix in vector format.
             ## Note: This operation can be reversed using unVecL()
             idx = 0
             N <- ncol(sqrMatrix)
             vM <- matrix(0,N*(N-1)/2,1)
             for (idxCol in seq(1,N-1)){
               for (idxRow in seq(idxCol+1,N)){
                 idx <- idx+1
                 vM[idx,1] <- sqrMatrix[idxRow,idxCol]
               }
             }
             return(as.vector(vM))
           }
)
## -- unVecl -- ####
setGeneric(name="unVecL",
           valueClass = "matrix",
           signature = c("lowerTri"),
           def = function(lowerTri){
             ## Returns a square matrix constructed using a vector of it's lower triangle.
             ## Note: This operation can be reversed using: vecL(Matrix)
             k <- length(lowerTri)
             N <- (1+sqrt(1+8*k))/2
             M <- matrix(0,N,N)
             M[lower.tri(M)] <- lowerTri
             return(M + t(M) + diag(N))
           }

)
## -- eigVec.EC -- ####
setGeneric(name=".eigVec.EC",
           valueClass = "matrix",
           signature = c("N"),
           def = function(N){
             # compute matrix of eigenvectors (columns of Q) for an EQUI-correlation model
             # N = number of series
             Q <- matrix(0,N,N)
             for (i in 1:N)
             {
               if (i==1) Q[,i] <- (1/sqrt(N))*matrix(1,N,1)
               else
               {
                 tmp <- 1/(sqrt((N-i+2)*(N-i+1)))
                 for (j in seq((i-1),N))
                 {
                   if (j==(i-1)) Q[j,i] <- tmp*(N-i+1)
                   else Q[j,i] <- tmp*(-1)
                 }
               }
             }
             return(Q) # NxN matrix
           }
)
## -- eigVal.EC -- ####
setGeneric(name=".eigVal.EC",
           valueClass = "numeric",
           signature = c("N","rho"),
           def = function(N,rho){
             # compute matrix of eigenvvalues for an EQUI-correlation model
             # N = number of series
             # rho = equicorrelation parameter (scalar)
             L <- rep((1-rho),N)
             L[1] <- L[1]+rho*N
             return(L) # Nx1
           }
)
## -- ar1.Filter -- ####
setGeneric(name=".ar1.Filter",
           valueClass = "matrix",
           signature = c("mX","vB"),
           def = function(mX,vB){
             # mX -- T x s matrix
             # vB -- 1 x s vector of coefficients
             # does AR type filtering with lag 1 only
             # output -- mY -- T x s matrix
             # mY[1,] = 0...0
             # mY[t,s] = mX[t,s]+vB[s]*mY[t-1,s]
             s <- length(vB)
             Tobs <- NROW(mX)
             if (NCOL(mX) != s) stop("Error in ar1.Filter: Number of columns in mX must equal the length of vB")
             mY <- matrix(mX[1,],1,s)
             for (t in 2:Tobs){
               mY <- rbind(mY,mX[t,]+vB*mY[(t-1),])
             }
             return(mY)
           }
)
## -- vec -- ####
setGeneric(name=".vec",
           valueClass = "matrix",
           signature = c("mat"),
           def = function(mat){
             # Convert a matrix to a single-column matrix
             return(matrix(as.vector(mat),ncol=1))
           }
)

## -- lag0 -- ####
setGeneric(name="lag0",
           valueClass = "matrix",
           signature = c("vX","lagRange"),
           def = function(vX,lagRange){
             lags <- matrix(0,nrow=NROW(vX),ncol=length(lagRange))
             #TODO: Replace with apply()
             for (i in 1:length(lagRange)){
               lags[,i] <- c(0*c(1:lagRange[i]),vX[(1:(NROW(vX)-lagRange[i]))])
             }
             return(lags)
           }
)

## -- sqrt_mat1 -- ####
setGeneric(name="sqrt_mat1",
           valueClass = "matrix",
           signature = c("m"),
           def = function(m){
             ## Note: This works for positive semi-definite (or definite) matrices that are diagonalizable (no normal Jordan forms, etc.)
             ##       The 'Matrix' package is used to try to compute the nearest P-D matrix, if required, and a warning is returned.

             N <- NROW(m)
             # Handle the special case where we have a 1 x 1 matrix:
             if(isTRUE(N==1)){
               m.sqrt <- matrix(sqrt(m[1,1]),nrow=1,ncol=1)
             } else {
               m.eig <- eigen(m)
               if (isTRUE(min(m.eig$values) <= 0) ) {
                 nPD <- NULL
                 nPD <- tryCatch(nearPD(m, base.matrix=TRUE, do2eigen = FALSE, conv.tol = 1e-15, maxit = 350),
                                 error = function(e){ message(e) }
                 )
                 if(!is.null(nPD)){ if(isTRUE(nPD$converged)){ m.eig <- eigen(nPD$mat) }
                                    else{
                                        warning("Matrix Inversion failed to converge")
                                        return(matrix())
                                    }
                                  }
               }
               m.sqrt <- m.eig$vectors %*% diag(sqrt(m.eig$values)) %*% solve(m.eig$vectors)
             }
             return(m.sqrt)
           }
)

## -- sqrt_mat2 -- ####
setGeneric(name="sqrt_mat2",
           valueClass = "list",
           signature = c("mat"),
           def = function(mat){
             # Using the Denman-Beavers algorithm:
             maxit <- 50
             stopifnot(nrow(mat) == ncol(mat))
             niter <- 0
             y <- mat
             z <- diag(rep(1,nrow(mat)))
             for (niter in 1:maxit) {
               y.temp <- 0.5*(y+solve(z))
               z <- 0.5*(z+solve(y))
               y <- y.temp
             }
             return(list(sqrt=y,sqrt.inv=z))
           }
)

invertMatrix_NPD = function(mat){

  # Try to invert 'mat'
  mat.inv <- NULL
  mat.inv <- try(solve(mat))
  if(!is.null(mat.inv)) { return(mat.inv) }
  else{
    # Create nearest PD approximation of 'mat'
    mat.nPD <- NULL
    mat.nPD <- try(nearPD(mat,corr = FALSE, base.matrix=TRUE, do2eigen = FALSE, conv.tol = 1e-18, maxit = 500))
      if(!is.null(mat.nPD)){
        if(isTRUE(mat.nPD$converged)){
          # Try to invert the nPD approximation
          mat.inv <- tryCatch(solve(mat.nPD$mat),
                              error = function(e){
                                        msg = cat("Attempt to invert the nearest PD approximation failed\n" ,e)
                                        warning(msg)
                                        return(NULL)
                                      }
                              )
          return(mat.inv)
        }
        else{ return(NULL) }  }
  }

}

