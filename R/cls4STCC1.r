## --- Correlation Structures --- ####
##
## -- We have a number of Correlation objects, with many common properties
## -- This class file maintains Structures for STCC1 (STCC with One transition)


CORRshape = list(single=1,double=2,double1loc=3)
CORRspeedopt = list(gamma=1,gamma_std=2,eta=3)

## --- Utiltiy Functions & Constants --- ####
#### ===============  Anna's Matrix Functions  ================== ####


myUnVecl <- function(vM){
  ## Returns a square matrix constructed using a vector of it's lower triangle.
  ## Note: This operation can be reversed using myVecl()
  k <- length(vM)
  N <- (1+sqrt(1+8*k))/2
  M <- matrix(0,N,N)
  M[lower.tri(M)] <- vM
  M <- M + t(M) + diag(N)
  #Return:
  M
}
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

myQ.EC <- function(N){
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
  #Return
  Q # NXN
}  #End: myQ.EC(N)


myL.EC <- function(N,rho){
  # compute matrix of eigenvvalues for an EQUI-correlation model
  # N = number of series
  # rho = equicorrelation parameter (scalar)
  L <- rep((1-rho),N)
  L[1] <- L[1]+rho*N
  #Return vector:
  L # Nx1
}

myFilter <- function(mX,vB){
  # mX -- T x s matrix
  # vB -- 1 x s vector of coefficients
  # does AR type filtering with lag 1 only
  # output -- mY -- T x s matrix
  # mY[1,] = 0...0
  # mY[t,s] = mX[t,s]+vB[s]*mY[t-1,s]
  s <- length(vB)
  Tobs <- NROW(mX)
  # GLEN: add error check, NCOL(mX)=s
  mY <- matrix(mX[1,],1,s)
  for (t in 2:Tobs){
    mY <- rbind(mY,mX[t,]+vB*mY[(t-1),])
  }
  # Return:
  mY
}

myVecd <- function(M){
  # vectorises main diagonal and off diagonal
  N <- NCOL(M)
  if (N==3){
    vM <- matrix(0,nrow=N*(N+1)/2,ncol=1)
    vM[1:3] <- diag(M)^2
    vM[4] <- sqrt(2)*M[3,2]
    vM[5] <- sqrt(2)*M[3,1]
    vM[6] <- sqrt(2)*M[2,1]
    return(vM)
  }
  else return(NULL)

}

#### ===============  End: Anna's Tricky Functions  ================== ###


## --- stcc1_class Definition --- ####
stcc1 <- setClass(Class = "stcc1_class",
               slots = c(N="integer",st="numeric",shape="integer",nr.covPars="integer",nr.trPars="integer",Tobs="integer"),
               contains = c("namedList")
               )

## Initialise with no params
setMethod("initialize","stcc1_class",
          function(.Object,...){
            .Object <- callNextMethod(.Object,...)
            .Object@N <- as.integer(2)
            # Return:
            .Object
          })

setGeneric(name="stcc1",
           valueClass = "stcc1_class",
           signature = c("N","st","shape"),
           def = function(N,st,shape){
             this <- new("stcc1_class")
             this@N <- as.integer(N)
             this@nr.covPars <- as.integer(N + (N^2-N)/2)
             this@shape <- as.integer(shape)
             if(shape==CORRshape$double) this@nr.trPars <- as.integer(3) else this@nr.trPars <- as.integer(2)
             this@st <- st
             this@Tobs <- as.integer(NROW(st))

             # Do validation checks:

             # End validation

             this@st <- st
             this@shape <- shape
             # Set Default Values:
             this$speedopt <- CORRspeedopt$eta
             if(shape==CORRshape$double) this$pars <- c(2.5,0.3,0.7) else this$pars <- c(2.5,0.5,NA)
             names(this$pars) <- c("speed","loc1","loc2")
             this$P1 <- matrix(0.2,N,N)
             diag(this$P1) <- 1
             this$P2 <- matrix(0.7,N,N)
             diag(this$P2) <- 1
             this$optimcontrol <- list(fnscale = -1, maxit = 1000, reltol = 1e-5)

             return(this)
           }
)

setGeneric(name=".calc.st_c",
           valueClass = "numeric",
           signature = c("stcc1Obj"),
           def = function(stcc1Obj){
             this <- stcc1Obj
             st_c <- 0
             if(this@shape == CORRshape$single) { st_c <- this@st-this$pars["loc1"] }
             if(this@shape == CORRshape$double) { st_c <- (this@st-this$pars["loc1"])*(this@st-this$pars["loc2"]) }
             if(this@shape == CORRshape$double1loc) { st_c <- (this@st-this$pars["loc1"])^2 }
             return(st_c)
           }
)

setGeneric(name=".calc.Gt",
           valueClass = "numeric",
           signature = c("stcc1Obj","st_c"),
           def = function(stcc1Obj,st_c){
             this <- stcc1Obj
             G <- 0
             if(this$speedopt == CORRspeedopt$gamma) { G <- 1/(1+exp(-this@speed*st_c)) }
             if(this$speedopt == CORRspeedopt$gamma_std) { G <- 1/(1+exp(-this@speed/sd(this@st)*st_c)) }
             if(this$speedopt == CORRspeedopt$eta) { G <- 1/(1+exp(-exp(this@speed)*st_c)) }

             return(matrix(G,nrow = this@Tobs,ncol = 1))
           }
)

## === .calc_Pt1() === ####
setGeneric(name=".calc_Pt1",
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

## --- Public Methods --- ####

## -- loglik.stcc1() ####
setGeneric(name="loglik.stcc1",
           valueClass = "numeric",
           signature = c("optimpars","z","stcc1Obj"),
           def = function(optimpars,z,stcc1Obj){

             # input: optimpars   -- c(speed,loc1,[loc2])
             #        z           -- volatility standardised returns (matrix TxN)
             #        stcc        -- object containing all the other parameters

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

             tmp.par <- optimpars

             vP1 <- tmp.par[1:this@nr.covPars]
             mP <- myUnVecl(vP1)
             eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
             # Check for SPD - positive-definite check:
             if (min(eig$values) <= 0) return(err_output)

             #Remove the P1 covPars, then extract the P2 covPars
             tmp.par <- tail(tmp.par,-this@nr.covPars)
             vP2 <- tmp.par[1:this@nr.covPars]
             mP <- myUnVecl(vP2)
             eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
             # Check for SPD - positive-definite check:
             if (min(eig$values) <= 0) return(err_output)

             st_c <- .calc.st_c(this)
             Gt <- .calc.Gt(this,st_c)

             Pt <- t(apply(Gt,MARGIN = 1,FUN = function(X,P1,P2) ((1-X)*P1 + X*P2), P1=vP1, P2=vP2))

             llt <- vector("numeric")
             for(t in 1:this@Tobs) {
               mPt <- myUnVecl(Pt[t,])
               llt[t] <- -0.5*log(det(mPt)) -0.5*( t(z[t,])%*%(qr.solve(mPt))%*%z[t,])
             }
             # Return:
             return(sum(llt))

           }
)

setGeneric(name="estimateSTCC1",
           valueClass = "stcc1_class",
           signature = c("z","stcc1Obj","estimationCtrl"),
           def = function(z,stcc1Obj,estimationCtrl){
             this <- stcc1Obj
             calcSE <- estimationCtrl$calcSE
             verbose <- estimationCtrl$verbose

             this$Estimated <- list()

             optimpars <- c( this$P1[lower.tri(this$P1)], this$P2[lower.tri(this$P2)], this$pars )
             optimpars <- na.omit(optimpars)

             ### ---  Call optim to calculate the estimate --- ###
             if (verbose) this$optimcontrol$trace <- 10

             tmp <- NULL
             try(tmp <- optim(optimpars,loglik.stcc1,z,this,gr=NULL,method="BFGS",control=this$optimcontrol,hessian=calcSE))

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
               this$Estimated$P1 <- myUnVecl(tmp.par[1:this@nr.covPars])
               tmp.par <- tail(tmp.par,-this@nr.covPars)
               this$Estimated$P2 <- myUnVecl(tmp.par[1:this@nr.covPars])
               tmp.par <- tail(tmp.par,-this@nr.covPars)
               this$Estimated$pars <- tmp.par
               if(this@shape != CORRshape$double) this$Estimated$pars <- c(this$Estimated$pars,NA)
               names(this$Estimated$pars) <- names(this$pars)

               if (calcSE) {
                 this$Estimated$hessian <- tmp$hessian
                 vecSE <- vector("numeric")
                 try(vecSE <- sqrt(-diag(qr.solve(tmp$hessian))))

                 if(length(vecSE) > 0) {
                   this$Estimated$P1.se <- myUnVecl(vecSE[1:this@nr.covPars])
                   vecSE <- tail(vecSE,-this@nr.covPars)
                   this$Estimated$P2.se <- myUnVecl(vecSE[1:this@nr.covPars])
                   vecSE <- tail(vecSE,-this@nr.covPars)

                   this$Estimated$pars.se <- vecSE
                   if(this$shape != CORRshape$double) this$Estimated$pars.se <- c(this$Estimated$pars.se,NA)
                   names(this$Estimated$pars.se) <- names(this$pars)
                 }
               }
               this$Estimated$Pt <- .calc_Pt1(this)

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


## --- PRIVATE METHODS --- ####

# calcParamStats_STCC <- function(STCC,sig_level_percent=5) {
#
#   if(sig_level_percent < 1) {
#     warning("Please enter an integer % value, e.g. 5 for 5%")
#     return()
#   }
#
#   pars <- STCC$Estimated$parsVector
#   errs <- STCC$Estimated$stderr
#
#   STCC$Estimated$tStats <- vector(mode="numeric",length=length(pars))
#   STCC$Estimated$PValues <- vector(mode="numeric",length=length(pars))
#   STCC$Estimated$fStat <- NULL
#
#   sig_level_percent = 5
#   sig_value = (100 - sig_level_percent)/100 * 2
#
#   # Calculate t-stat = test if param is significantly different from zero
#   STCC$Estimated$tStats <- round(pars/errs,6)
#
#   # Calculate P-Values based on the above t-Stats, for the significance level provided
#   FUN1 <- function(x) {round(sig_value*pnorm(-abs(x)),6)}
#   STCC$Estimated$PValues <- vapply(STCC$Estimated$tStats, FUN1, vector("numeric",length = 1))
#
#   # Add an asterix to each P-Value (in the name), to help users understand the parameter significance
#   FUN2 <- function(x) {ifelse((abs(x) < sig_level_percent/100), "*", "") }
#   names(STCC$Estimated$PValues) <- paste0(names(STCC$Estimated$PValues),vapply(STCC$Estimated$PValues, FUN2, vector("character",length = 1)))
#
#   #Return
#   STCC
#
# }


## --- Override Methods --- ####

# ## -- plot() ####
# setMethod("plot",signature = c(x="stcc1_class",y="missing"),
#           function(x, y, ...){
#             this <- x
#             plot(this@h, type='l', ylab = "Cond.Variance", ...)
#           })

