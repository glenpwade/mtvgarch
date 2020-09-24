## --- Correlation Structures --- ####
##
## -- We have a number of Correlation objects, with many common properties
## -- This class file maintains the Structure for STCC1 (STCC with One Transition)


##===============  Anna's Matrix Functions  ==================##

setGeneric(name=".Vecl",
           valueClass = "numeric",
           signature = c("sqrMatrix"),
           def = function(sqrMatrix){
             ## Returns the lower triangle of a square matrix in vector format.
             ## Note: This operation can be reversed using myUnVecl()
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
setGeneric(name=".unVecl",
           valueClass = "matrix",
           signature = c("lowerTri"),
           def = function(lowerTri){
             ## Returns a square matrix constructed using a vector of it's lower triangle.
             ## Note: This operation can be reversed using: Matrix[lower.tri(Matrix)]
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

## ===============  End: Anna's Tricky Functions  ================== ##


## --- stcc1_class Definition --- ####

corrtype <- list(CCC=1,CEC=2,STCC1=3,STEC1=4)
corrshape <- list(single=1,double=2,double1loc=3)
corrspeedopt <- list(gamma=1,gamma_std=2,eta=3)

stcc1 <- setClass(Class = "stcc1_class",
               slots = c(st="numeric",nr.covPars="integer",nr.trPars="integer",Tobs="integer",N="integer"),
               contains = c("namedList")
               )

## --- Initialise --- ####
setMethod("initialize","stcc1_class",
          function(.Object,...){
            .Object <- callNextMethod(.Object,...)
            # By definition of STCC1:
            .Object$type <- corrtype$STCC1

            # Default initial values
            .Object@N <- as.integer(0)
            .Object$shape <- corrshape$single
            .Object$speedopt <- corrspeedopt$eta
            .Object$optimcontrol <- list(fnscale = -1, maxit = 1000, reltol = 1e-5)

            # Return:
            .Object
          })

## -- Constructor:stcc1 -- ####
setGeneric(name="stcc1",
           valueClass = "stcc1_class",
           signature = c("mtvgarchObj"),
           def = function(mtvgarchObj){
             this <- new("stcc1_class")

             ## -- Do validation checks -- ####
             objType <- class(mtvgarchObj)
             if(objType[1] != "mtvgarch_class"){
               warning("a valid instance of the mtvgarch_class is required to create an STCC1 model")
               return(this)
             }
             # End validation

             # Set Default Values:
             this@N <- mtvgarchObj@N
             this@st <- 1:mtvgarchObj@Tobs/mtvgarchObj@Tobs
             this@Tobs <- mtvgarchObj@Tobs

             if(this$shape==corrshape$double) {
               this@nr.trPars <- as.integer(3)
               this$pars <- c(2.5,0.33,0.66)
             }else {
               this@nr.trPars <- as.integer(2)
               this$pars <- c(2.5,0.5,NA)
             }
             names(this$pars) <- c("speed","loc1","loc2")

             N <- this@N
             this@nr.covPars <- as.integer((N^2-N)/2)
             this$P1 <- matrix(0.3,N,N)
             diag(this$P1) <- 1
             this$P2 <- matrix(0.7,N,N)
             diag(this$P2) <- 1

             return(this)
           }
)
## -- calc.st_c -- ####
setGeneric(name=".calc.st_c",
           valueClass = "numeric",
           signature = c("stcc1Obj"),
           def = function(stcc1Obj){
             this <- stcc1Obj
             loc1 <- this$Estimated$pars["loc1"]
             loc2 <- this$Estimated$pars["loc2"]

             st_c <- 0
             if(this$shape == corrshape$single) { st_c <- this@st - loc1 }
             if(this$shape == corrshape$double) { st_c <- (this@st - loc1)*(this@st - loc2) }
             if(this$shape == corrshape$double1loc) { st_c <- (this@st - loc1)^2 }
             return(st_c)
           }
)

## -- calc.Gt -- ####
setGeneric(name=".calc.Gt",
           valueClass = "matrix",
           signature = c("stcc1Obj","st_c"),
           def = function(stcc1Obj,st_c){
             this <- stcc1Obj
             speed <- this$Estimated$pars["speed"]

             G <- 0
             if(this$speedopt == corrspeedopt$gamma) { G <- 1/(1+exp(-speed*st_c)) }
             if(this$speedopt == corrspeedopt$gamma_std) { G <- 1/(1+exp(-speed/sd(this@st)*st_c)) }
             if(this$speedopt == corrspeedopt$eta) { G <- 1/(1+exp(-exp(speed)*st_c)) }

             return(matrix(G,nrow = this@Tobs,ncol = 1))
           }
)

## -- calc.Pt -- ####
setGeneric(name=".calc.Pt",
           valueClass = "matrix",
           signature = c("stcc1Obj"),
           def =   function(stcc1Obj){

             this <- stcc1Obj
             mP1 <- this$Estimated$P1
             vP1 <- mP1[lower.tri(mP1)]
             mP2 <- this$Estimated$P2
             vP2 <- mP2[lower.tri(mP2)]

             st_c <- .calc.st_c(this)
             Gt <- .calc.Gt(this,st_c)
             Pt <- t(apply(Gt,MARGIN = 1,FUN = function(X,P1,P2) ((1-X)*P1 + X*P2), P1=vP1, P2=vP2))

             return(Pt)

           }
)

## -- loglik.stcc1() --####
setGeneric(name=".loglik.stcc1",
           valueClass = "numeric",
           signature = c("optimpars","z","stcc1Obj"),
           def = function(optimpars,z,stcc1Obj){

             err_output <- -1e10
             this <- stcc1Obj

             #### ======== constraint checks ======== ####

             # # Check 1: Confirm we have a valid shape


             # # Check 2: Check the boundary values for speed params:
             # #speedoption -- 1=gamma, 2=gamma/std(st), 3=exp(eta), 4=1/lambda^2
             # maxSpeed <- switch(speedoption,1000,(1000/sd(st)),7.0,0.30)
             # if (speed > maxSpeed) return(err_output)
             # if (speed < 0) return(err_output)
             #
             # # Check 3: Check the locations fall within min-max values of st
             # # We must have at least 1 loc1 to be inside this shape..loop, so no need to check if loc1 contains a valid value:
             # if (loc1 < min(st)) return(err_output)
             # if (loc1 > max(st)) return(err_output)
             # if (numTRpars==3) {
             #   if (loc2 < min(st)) return(err_output)
             #   if (loc2 > max(st)) return(err_output)
             # }


             #### ======== calculate loglikelihood ======== ####

             this$Estimated$pars <- tail(optimpars,this@nr.trPars)

             tmp.par <- optimpars

             vP1 <- tmp.par[1:this@nr.covPars]
             mP <- .unVecl(vP1)
             eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
             # Check for SPD - positive-definite check:
             if (min(eig$values) <= 0) return(err_output)

             #Remove the P1 covPars, then extract the P2 covPars
             tmp.par <- tail(tmp.par,-this@nr.covPars)
             vP2 <- tmp.par[1:this@nr.covPars]
             mP <- .unVecl(vP2)
             eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
             # Check for SPD - positive-definite check:
             if (min(eig$values) <= 0) return(err_output)

             st_c <- .calc.st_c(this)
             Gt <- .calc.Gt(this,st_c)

             Pt <- t(apply(Gt,MARGIN = 1,FUN = function(X,P1,P2) ((1-X)*P1 + X*P2), P1=vP1, P2=vP2))

             llt <- vector("numeric")
             for(t in 1:this@Tobs) {
               mPt <- .unVecl(Pt[t,])
               llt[t] <- -0.5*log(det(mPt)) -0.5*( t(z[t,])%*%(qr.solve(mPt))%*%z[t,])
             }
             # Return:
             return(sum(llt))

           }
)
## --- estimateSTCC1 --- ####
setGeneric(name="estimateSTCC1",
           valueClass = "stcc1_class",
           signature = c("z","stcc1Obj","estimationCtrl"),
           def = function(z,stcc1Obj,estimationCtrl){
             this <- stcc1Obj

             calcSE <- estimationCtrl$calcSE
             verbose <- estimationCtrl$verbose

             this$Estimated <- list()

             optimpars <- c( this$P1[lower.tri(this$P1)], this$P2[lower.tri(this$P2)], this$pars )
             optimpars <- optimpars[!is.na(optimpars)]

             ### ---  Call optim to calculate the estimate --- ###
             if (verbose) this$optimcontrol$trace <- 10

             tmp <- NULL
             try(tmp <- optim(optimpars,.loglik.stcc1,z,this,gr=NULL,method="BFGS",control=this$optimcontrol,hessian=calcSE))

             ### ---  Interpret the response from optim --- ###
             # An unhandled error could result in a NULL being returned by optim()
             if (is.null(tmp)) {
               this$Estimated$value <- -1e10
               this$Estimated$error <- TRUE
               return(this)
             }

             #Optim converged successfully => we expect tmp$par to have good estimates!
             if (tmp$convergence==0) {
               this$Estimated$error <- FALSE
               this$Estimated$value <- tmp$value

               tmp.par <- tmp$par
               this$Estimated$P1 <- .unVecl(tmp.par[1:this@nr.covPars])
               tmp.par <- tail(tmp.par,-this@nr.covPars)
               this$Estimated$P2 <- .unVecl(tmp.par[1:this@nr.covPars])
               tmp.par <- tail(tmp.par,-this@nr.covPars)
               this$Estimated$pars <- tmp.par
               if(this$shape != corrshape$double) this$Estimated$pars <- c(this$Estimated$pars,NA)
               names(this$Estimated$pars) <- names(this$pars)

               if (calcSE) {
                 this$Estimated$hessian <- tmp$hessian
                 vecSE <- vector("numeric")
                 try(vecSE <- sqrt(-diag(qr.solve(tmp$hessian))))

                 if(length(vecSE) > 0) {
                   this$Estimated$P1.se <- .unVecl(vecSE[1:this@nr.covPars])
                   vecSE <- tail(vecSE,-this@nr.covPars)
                   this$Estimated$P2.se <- .unVecl(vecSE[1:this@nr.covPars])
                   vecSE <- tail(vecSE,-this@nr.covPars)

                   this$Estimated$pars.se <- vecSE
                   if(this$shape != corrshape$double) this$Estimated$pars.se <- c(this$Estimated$pars.se,NA)
                   names(this$Estimated$pars.se) <- names(this$pars)
                 }
               }
               this$Estimated$Pt <- .calc.Pt(this)

             } else {
               #Failed to converge
               this$Estimated$error <- TRUE
               this$Estimated$value <- -1e10
               this$Estimated$optimoutput <- tmp
             }

             if (verbose) this$Estimated$optimoutput <- tmp
             #Return:
             return(this)
           }
)
setMethod("estimateSTCC1",signature = c("matrix","stcc1_class","missing"),
          function(z,stcc1Obj){
            estimationControl <- list(calcSE = TRUE,verbose = TRUE)
            estimateSTCC1(z,stcc1Obj,estimationControl)
          })

## --- PRIVATE METHODS --- ####
