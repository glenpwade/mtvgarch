
##--- All the general, common functions that don't belong to any one specific class in the MTVGARCH Package ---##

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
             ## Note: this only works for positive semi-definite (or definite) matrices that are diagonalizable
             ## (no normal Jordan forms, etc.)

             # Handle the special case where we have a 1 x 1 matrix:
             if(NROW(m) == 1){
               m.sqrt <- matrix(sqrt(m[1,1]),nrow=1,ncol=1)
             } else {
               m.eig <- eigen(m)
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

## -- generateRefData -- ####
setGeneric(name="generateRefData",
           valueClass = "matrix",
           signature = c("nr.series","nr.obs","tvObj","garchObj","corrObj","noiseDist"),
           def =  function(nr.series,nr.obs,tvObj,garchObj,corrObj,noiseDist)
           {
             ## TODO: 1. Override function to handle TV & GARCH as optional params
             ## TODO: 2. Add more validation / better defaults (and messages) for noiseDist

             ## Note:
             ## noiseDist is a named-list describing the error-distribution and parameters
             ## e.g. noiseDist$name = 'Normal'     noiseDist$mean = 0  noiseDist$sd = 1
             ## or   noiseDist$name = 'Student-t'  noiseDist$df = 6    noiseDist$ncp = 0

             # Generate Noise Data:
             if(identical(toupper(substr(trimws(noiseDist$name),1,1)),"S") ){
               # Standardised Student-t Error/Noise Distribution
               if(!is.null(noiseDist$df)) {
                 df <- noiseDist$df
               } else {
                 df <- 6
                 message("param noiseDist did not contain a valid $df value.  df set = 6 ")
               }
               u <- matrix(rt(nr.obs * nr.series,df),nrow=nr.obs, ncol=nr.series)
               u <- u/sqrt(df/(df-2))
             }
             else
             {
               # Normal Error/Noise Distribution (Default is Standard-Normal)
               u <- matrix(rnorm(nr.obs * nr.series),nrow=nr.obs, ncol=nr.series)
             }

             # Step1: Create Correlation

             # e = correlated data
             e <- u

             corrType <- class(corrObj)

             # - - - CCC - - -
             if (corrType[1] == "ccc_class"){
               if(is.null(corrObj$Estimated)) {
                 P <- corrObj$P
               } else P <- corrObj$Estimated$P
               P.sqrt <-  sqrt_mat1(P)
               e <- t(P.sqrt %*% t(u))
             }

             if (corrType[1] == "stcc1_class"){
               if(is.null(corrObj$Estimated)) {
                 corrObj$Estimated$P1 <- corrObj$P1
                 corrObj$Estimated$P2 <- corrObj$P2
                 corrObj$Estimated$pars <- corrObj$pars
               }
               Pt <- .calc.Pt(corrObj)
               for (t in 1:corrObj@Tobs){
                 mPt <- unVecL(Pt[t,,drop=FALSE])
                 mPt.sqrt <- sqrt_mat1(mPt)
                 e[t,] <- t( mPt.sqrt %*% t(u[t,,drop=FALSE]) )
               }
             }

             #End: Generate Correlated Data

             # Step2: Inject GARCH into Data

             garchObj$pars["omega",1] <- ( 1 - garchObj$pars["alpha",1] - garchObj$pars["beta",1] )
             discardObs <- 1500
             discardData <- matrix(rnorm(discardObs * nr.series),nrow=discardObs, ncol=nr.series)
             e <- rbind(discardData,e)
             endRow <- discardObs + nr.obs

             for (b in 1:nr.series){
               w <- z <- e[,b]
               ht_1 <- 1
               w[1] <- z[1]
               for (t in 2:endRow) {
                 ht <- garchObj$pars["omega",1] + garchObj$pars["alpha",1]*(w[t-1])^2 + garchObj$pars["beta",1]*ht_1
                 if(garchObj$type == garchtype$gjr) { ht <- ht + garchObj$pars["gamma",1]*(min(w[t-1],0)^2) }
                 ht_1 <- ht
                 w[t] <- sqrt(ht)*z[t]
               }
               e[,b] <- as.numeric(w)
             }

             # Discard the first 2000
             startRow <- discardObs + 1
             e <- e[(startRow:endRow), ]

             # Step3: Inject TV into Data
             gt <- get_g(tvObj)
             e <- e*sqrt(gt)

             #Return:
             e

           }
)



