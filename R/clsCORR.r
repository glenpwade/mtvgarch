## --- Correlation Structures --- ####
##
## -- We have a number of Correlation objects, with many common properties
## -- This class file maintains all these objects, for example:
## -- CCC, CEC, DCC, STCC, STEC

## --- When created:

# Slots (internal variables for use in methods - should only be set by pkg_code)

    # corr@st             -- numeric: conditional variances
    # corr@shape          -- numeric: similar to tv shape

# properties (external variables, visible to user)

    # corr$type                  -- CCC, CEC, DCC, STCC, STEC
    # corr$pars                  -- matrix: structure depends on type
    # corr$optimcontrol          -- list: wrapper for the 'control' parameter sent to optim()

## --- After estimation:

    # corr$Estimated$value        -- scalar: log-likelihood value
    # corr$Estimated$error        -- logical: T/F indication if an error occurred during estimation
    # corr$Estimated$pars 
    # corr$Estimated$hessian 
    # corr$Estimated$se 
    # corr$Estimated$optimoutput  -- list: the full returned value from optim().  Only set if the verbose param = True

STCCtype = list(none=0,CCC=1,STCC=2)
TVshape = list(none=0,delta0only=1,single=2,double=3,double1loc=4)
TRspeedopt = list(none=0,gamma=1,gamma_std=2,eta=3,lamda2_inv=4)
err_output <- -1e10

## --- Utiltiy Functions & Constants --- ####
#### ===============  Anna's Tricky Matrix Functions  ================== ####


myVecl <- function(M){
  ## Returns the lower triangle of a square matrix in vector format.
  ## Note: This operation can be reversed using myUnVecl()
  idx = 0
  N <- ncol(M)
  vM <- matrix(0,N*(N-1)/2,1)
  for (idxCol in seq(1,N-1)){
    for (idxRow in seq(idxCol+1,N)){
      idx <- idx+1
      vM[idx,1] <- M[idxRow,idxCol]
    }
  }
  #Return:
  vM
}

myUnVecl <- function(vM){
  ## Returns a square matrix constructed using a vector of it's lower triangle.
  ## Note: This operation can be reversed using myVecl()
  k <- length(vM)
  N <- (1+sqrt(1+8*k))/2
  M <- matrix(0,N,N)
  M[lower.tri(M,diag=FALSE)] <- vM
  M <- M + t(M) + diag(N)
  #Return:
  M
}


myQ.EC <- function(N)
{
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

#### ===============  Anna's Tricky Functions  ================== ###




## --- corr_class Definition --- ####
corr <- setClass(Class = "corr_class",
               slots = c(h="numeric",nr.pars="integer",type="integer"), 
               contains = c("namedList")
               )

## Initialise with no params
setMethod("initialize","corr_class",
          function(.Object,...){
            .Object <- callNextMethod(.Object,...)
            .Object@nr.pars <- as.integer(0)
            .Object@type <- as.integer(0)
            # Return:
            .Object
          })

setGeneric(name="corr",
           function(type="numeric",order="numeric")standardGeneric("corr"),
           valueClass = "corr_class",
           signature = c("type","order")
          )

## Constructor with 1 param - creates a Garch(1,1)
setMethod("corr",signature = c("numeric","missing"),
          function(type){
            # Create a CORR(CCC) model
            #this <- corr(type,c(1,1))
            return(this)
          })

## Constructor with 2 params
setMethod("corr",signature = c("numeric","numeric"),
          function(type,order){
           
            this <- new("corr_class")
            this$type <- type
            
            this$optimcontrol <- list(fnscale = -1, maxit = 1500, reltol = 1e-7)

            return(this)
          })



## --- STCC class Definition  ####
# Define the GARCH object class:
newSTCC <- setClass("stcc_class",contains = c("namedList"))

setMethod("initialize","stcc_class",
          function(.Object,...){
            # Initialise fields with defaults:
            .Object$type <- STCCtype$STCC
            .Object$Tobs <- 0
            .Object$st <- NA
            .Object$P1 <- diag(4)
            .Object$P2 <- matrix(1,4,4)
            .Object$TRpars <- c(2.0,0.5,0.7)  #speed = 2, location = 0.5
            .Object$shape <- TVshape$single
            .Object$speedoption <- TRspeedopt$eta    
            .Object$optimcontrol <- list(fnscale = -1)
            # Return:
            .Object
          })

setGeneric(name="newSTCC",def=function(Tobs,st,type,p1,p2,trpars,shape,speedopt,...){
  standardGeneric("newSTCC")
}, valueClass = "stcc_class", signature = c("Tobs","st","type","p1","p2","trpars","shape","speedopt"))


# Create default constructor:
sig1 <- c("missing","numeric","numeric","matrix","matrix","numeric","numeric","numeric")
setMethod(f="newSTCC",signature = sig1,
          definition = function(st,type,p1,p2,trpars,shape,speedopt, ...){
            # Validate params
            if (is.null(st)) stop("Error: missing 'st'\nAll parameters are required to construct an STCC object")
            if (is.null(type)) stop("Error: missing 'type'\nAll parameters are required to construct an STCC object")
            if (is.null(p1)) stop("Error: missing 'p1'\nAll parameters are required to construct an STCC object")
            if (is.null(p2)) stop("Error: missing 'p2'\nAll parameters are required to construct an STCC object")
            if (is.null(trpars)) stop("Error: missing 'trpars'\nAll parameters are required to construct an STCC object")
            if (is.null(shape)) stop("Error: missing 'shape'\nAll parameters are required to construct an STCC object")
            if (is.null(speedopt)) stop("Error: missing 'speedopt'\nAll parameters are required to construct an STCC object")
            #
            stcc <- new("stcc_class") 
            stcc$st <- st
            stcc$Tobs <- NROW(st)
            stcc$type <- type
            stcc$P1 <- p1
            stcc$P2 <- p2
            stcc$TRpars <- trpars
            stcc$shape <- shape
            stcc$speedoption <- speedopt
            # Return:
            stcc
          })

# Create convenience 'quick' constructor: Linear Transition from 0 to 1, using eta
sig2 <- c("numeric","missing","missing","missing","missing","missing","missing","missing")
setMethod(f="newSTCC",signature = sig2,
          definition = function(Tobs){
            # Validate params
            if (is.null(Tobs)) stop("Error: missing 'Tobs'\nPlease enter a valid number of observations to create an STCC object")
            #
            stcc <- new("stcc_class") 
            stcc$Tobs <- Tobs
            stcc$st <- 1:Tobs/Tobs
            # Return:
            stcc
          })



## --- Public Methods --- ####

EstimateSTCC <- function(z,stcc,calcHess=FALSE,verbose=FALSE) {
  
  ###
  ### ---  Call optim to calculate the estimate --- ###
  ###
  
  STCC <- stcc
  numCovPars <- NROW(myVecl(STCC$P1))
  optimpars <- c(myVecl(STCC$P1),myVecl(STCC$P2),STCC$TRpars)
  if (verbose) STCC$optimcontrol$trace <- 10
  tmp <- NULL
  try(tmp <- optim(optimpars,myLogLik.stcc,z,STCC,gr=NULL,method="BFGS",control=STCC$optimcontrol,hessian=calcHess))
  
  ### ---  Interpret the response from optim --- ###
  # An unhandled error could result in a NULL being returned by optim()
  if (is.null(tmp)) {
    STCC$Estimated$value <- err_output 
    STCC$Estimated$error <- TRUE
    return(STCC)
  }
  
  #Optim converged successfully => we expect tmp$par to have good estimates!
  if (tmp$convergence==0) {
    if(STCC$shape==TVshape$double) numTRpars <- 3 else numTRpars <- 2
    if (!is.null(STCC$P1)) STCC$Estimated$P1 <- myUnVecl(tmp$par[1:numCovPars])
    if (!is.null(STCC$P2)) STCC$Estimated$P2 <- myUnVecl(tmp$par[(numCovPars+1):(2*numCovPars)])
    if (!is.null(STCC$TRpars)) STCC$Estimated$TRpars <- tail(tmp$par,numTRpars)
    # TODO: fix this parsVector code - remove try()
    try(STCC$Estimated$parsVector <- c(myVecl(STCC$Estimated$P1),myVecl(STCC$Estimated$P2),STCC$Estimated$TRpars))
    
    STCC$Estimated$value <- tmp$value
    STCC$Estimated$condcorrs <- myLogLik.stcc(tmp$par,z,STCC,return_ll=FALSE)
    STCC$Estimated$error <- FALSE
    if (calcHess) STCC$Estimated$hessian <- tmp$hessian
    
  } else { 
    #Failed to converge
    STCC$Estimated$error <- TRUE
    STCC$Estimated$value <- err_output
    STCC$Estimated$optimoutput <- tmp
  }
  
  if (verbose) STCC$Estimated$optimoutput <- tmp
  #Return:
  STCC
  
}  #End: EstimateSTCC()



## --- PRIVATE METHODS --- ####

calcStderr_STCC <- function(e,stcc) {
  STCC <- stcc
  STCC$Estimated$stderr <- NULL
  
  if(is.null(STCC$Estimated$hessian)) {
    warning("This method of generating Standard Errors is unreliable. \nPlease re-estimate using the 'calcHess=TRUE' parameter")
    optimpars <- NULL
    optimpars <- c(myVecl(STCC$Estimated$P1),myVecl(STCC$Estimated$P2),STCC$Estimated$TRpars)
    STCC$Estimated$hessian <- optimHess(optimpars,myLogLik.tv_univar,gr=NULL,e,STCC,control=STCC$optimcontrol)
  }
  
  # Using LU decomp:
  STCC$Estimated$stderr <- sqrt(-diag(solve(STCC$Estimated$hessian)))
  # Using Choleski
  # try(TV$Estimated$stderr <- sqrt(-diag(chol2inv(chol(hess)))))  
  if (is.null(STCC$Estimated$stderr)) {
    msgWarning <- "Failed to calculate Standard Error"
    warning(msgWarning)
  } else STCC$Estimated$stderr <- round(STCC$Estimated$stderr,6)
  
  #Return
  STCC
}

calcParamStats_STCC <- function(STCC,sig_level_percent=5) {
  
  if(sig_level_percent < 1) {
    warning("Please enter an integer % value, e.g. 5 for 5%")
    return()
  }
  
  pars <- STCC$Estimated$parsVector
  errs <- STCC$Estimated$stderr
  
  STCC$Estimated$tStats <- vector(mode="numeric",length=length(pars))
  STCC$Estimated$PValues <- vector(mode="numeric",length=length(pars))
  STCC$Estimated$fStat <- NULL
  
  sig_level_percent = 5
  sig_value = (100 - sig_level_percent)/100 * 2
  
  # Calculate t-stat = test if param is significantly different from zero
  STCC$Estimated$tStats <- round(pars/errs,6)
  
  # Calculate P-Values based on the above t-Stats, for the significance level provided
  FUN1 <- function(x) {round(sig_value*pnorm(-abs(x)),6)}
  STCC$Estimated$PValues <- vapply(STCC$Estimated$tStats, FUN1, vector("numeric",length = 1))
  
  # Add an asterix to each P-Value (in the name), to help users understand the parameter significance
  FUN2 <- function(x) {ifelse((abs(x) < sig_level_percent/100), "*", "") }
  names(STCC$Estimated$PValues) <- paste0(names(STCC$Estimated$PValues),vapply(STCC$Estimated$PValues, FUN2, vector("character",length = 1)))
  
  #Return
  STCC
  
}


myLogLik.stcc <- function(optimpars,z,stcc,return_ll=TRUE){
  # model: STCC
  # input: pars        -- c(speed,loc1,[loc2])
  #        z           -- volatility standardised returns (matrix TxN)
  #        stcc        -- list containing all the other parameters
  
  STCC <- stcc
  speedoption <- STCC$speedoption
  shape <- STCC$shape  
  st <- STCC$st
  
  Tobs <- NROW(z)
  N <- NCOL(z)
  
  #### ======== constraint checks ======== ####
  
  # Check 1: Confirm we have a valid shape & extract speed & location:
  if (shape==TVshape$none) return(err_output)
  if(shape==TVshape$double) numTRpars <- 3 else numTRpars <- 2
  TRpars <- tail(optimpars,numTRpars)
  speed <- TRpars[1]
  loc1 <- TRpars[2]
  if(numTRpars==3) loc2 <- TRpars[3] else loc2 <- NA
  
  # Check 2: Check the boundary values for speed params:
  #speedoption -- 1=gamma, 2=gamma/std(st), 3=exp(eta), 4=1/lambda^2
  maxSpeed <- switch(speedoption,1000,(1000/sd(st)),7.0,0.30)
  if (speed > maxSpeed) return(err_output)
  if (speed < 0) return(err_output)
  
  # Check 3: Check the locations fall within min-max values of st
  # We must have at least 1 loc1 to be inside this shape..loop, so no need to check if loc1 contains a valid value:
  if (loc1 < min(st)) return(err_output)
  if (loc1 > max(st)) return(err_output)
  if (numTRpars==3) {
    if (loc2 < min(st)) return(err_output)
    if (loc2 > max(st)) return(err_output)  
  }
  
  
  #### ======== calculate loglikelihood ======== ####
  # - - - CCC - - -
  if (STCC$type==STCCtype$CCC){
    vP <- c(0,0,0) #pars
    mP <- myUnVecl(vP)
    eig <- eigen(mP,symmetric=TRUE,only.values = TRUE)
    if (min(eig$values) <= 0) return(err_output)
    
    # - - - P(t) and loglik-value
    llt <- rep(0,Tobs)
    mPinv <- solve(mP)
    const <- -0.5*log(det(mP))
    for(t in seq(1,Tobs)) llt[t] <- const - 0.5*(t(z[t,])%*%(mPinv)%*%z[t,])
  }
  
  # - - - STCC - - -
  if (STCC$type==STCCtype$STCC) {
    numCovPars <- NROW(myVecl(STCC$P1))
    vP1 <- optimpars[1:numCovPars]
    vP2 <- optimpars[(numCovPars+1):(2*numCovPars)]
    
    mP1 <- myUnVecl(vP1)
    eig1 <- eigen(mP1,symmetric=TRUE,only.values=TRUE)
    # Check for SPD - positive-definite check:
    if (min(eig1$values) <= 0) return(err_output)
    mP2 <- myUnVecl(vP2)
    eig2 <- eigen(mP2,symmetric=TRUE,only.values=TRUE)
    # Check for SPD - positive-definite check:
    if (min(eig2$values) <= 0) return(err_output)
    
    # - - - Calculate G(t) - - -
    # if(shape == TVshape$none){}
    # if(shape == TVshape$delta0only){}
    # if(shape == TVshape$single){st_c <- st-loc1}
    # if(shape == TVshape$double){st_c <- (st-loc1)*(st-loc2)}
    # if(shape == TVshape$double1loc){st_c <- (st-loc1)^2}
    st_c <- switch(shape,st,st_c <- st-loc1,(st-loc1)*(st-loc2),(st-loc1)^2)
    G_inv <- switch(speedoption, 1+exp(-speed*st_c), 1+exp(-speed/sd(st)*st_c), 1+exp(-exp(speed)*st_c))
    Gt <- matrix(1/G_inv,nrow = Tobs,ncol = 1)
    
    # - - - P(t) and loglik-value
    #Pt <- matrix(0,nrow=Tobs,ncol=(N*(N-1)/2))
    llt <- NULL
    #calcPt <- function(X,P1,P2) ((1-X)*P1 + X*P2)
    #Pt <- t(apply(Gt,MARGIN = 1,FUN = calcPt, P1=vP1, P2=vP2))
    Pt <- t(apply(Gt,MARGIN = 1,FUN = function(X,P1,P2) ((1-X)*P1 + X*P2), P1=vP1, P2=vP2))
    for(t in seq(1,Tobs))
    {
      #Pt[t,] <- (1-Gt[t])*vP1 + Gt[t]*vP2  #This line has been replaced by the apply() function above
      mPt <- myUnVecl(Pt[t,])
      mPtinv <- chol2inv(chol(mPt))
      llt[t] <- -0.5*log(det(mPt)) -0.5*( t(z[t,])%*%(mPtinv)%*%z[t,])
    } # End: for loop
    
  } # End: if(STCC$type==STCCtype$STCC)
  
  # Return: 
  if (return_ll) return(sum(llt)) else {
    if (STCC$type==STCCtype$CCC) Pt <- matrix(vP,nrow=Tobs,ncol=length(vP),byrow=TRUE)
    return(Pt)
  }
  
}  #End: myLogLik.stcc()

## --- Override Methods --- ####

## -- plot() ####
setMethod("plot",signature = c(x="corr_class",y="missing"),
          function(x, y, ...){
            this <- x
            plot(this@h, type='l', ylab = "Cond.Variance", ...)
          })

