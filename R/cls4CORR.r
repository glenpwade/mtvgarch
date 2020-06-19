## --- Correlation Structures --- ####
##
## -- We have a number of Correlation objects, with many common properties
## -- This class file maintains all these objects, for example:
## -- CCC, CEC, DCC, STCC, STEC
## --- Structures for STCC.1 (One transition), STCC.2 (Two transitions, same transition variable),


CORRtype = list(CCC=1,STCC.1=2,STCC.2=3,DSTCC=4)
CORRshape = list(single=1,double=2,double1loc=3)
CORRspeedopt = list(gamma=1,gamma_std=2,eta=3)
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


## --- corr_class Definition --- ####
corr <- setClass(Class = "corr_class",
               slots = c(N="integer",st="numeric",nr.covPars="integer",Tobs="integer"), 
               contains = c("namedList")
               )

## Initialise with no params
setMethod("initialize","corr_class",
          function(.Object,...){
            .Object <- callNextMethod(.Object,...)
            .Object@N <- as.integer(0)
            # Return:
            .Object
          })

setGeneric(name="corr",
           function(type="numeric",N="numeric",st="numeric")standardGeneric("corr"),
           valueClass = "corr_class",
           signature = c("type","N","st")
          )

## Constructor - creates a correlation object as specified
setMethod("corr",signature = c("numeric","numeric","numeric"),
          function(type,N,st){
            this <- new("corr_class")
            this$type <- type
            this@N <- as.integer(N)
            this@st <- st
            
            if(type==CORRtype$CCC){
              this$P <- diag(nrow = N)
              this$optimcontrol <- list(fnscale = -1, maxit = 1500, reltol = 1e-5)
              return(this)
            }
            if(type==CORRtype$STCC.1){
            
               # Do validation checks:
               
               # End validation
               this$speedopt <- CORRspeedopt$eta
               this$st <- st
               this$shape <- CORRshape$single
               this$pars <- c(2,0.5,NA)
               names(this$pars) <- c("speed","loc1","loc2")
               this$P1 <- matrix(0.2,N,N)
               diag(this$P1) <- 1
               this$P2 <- matrix(0.7,N,N)
               diag(this$P2) <- 1
               this$optimcontrol <- list(fnscale = -1, maxit = 1500, reltol = 1e-5)
               return(this)
            }
            if(type==CORRtype$STCC.2){
              
              # Do validation checks:
              
              # End validation
              this$speedopt <- CORRspeedopt$eta
              this$st <- st
              this$shape <- c(CORRshape$single,CORRshape$single)
              this$pars1 <- c(2,0.3,NA)
              names(this$pars1) <- c("speed","loc1","loc2")
              this$pars2 <- c(2,0.6,NA)
              names(this$pars2) <- c("speed","loc1","loc2")
              this$P1 <- matrix(0.2,N,N)
              diag(this$P1) <- 1
              this$P2 <- matrix(0.7,N,N)
              diag(this$P2) <- 1
              this$P3 <- matrix(0.5,N,N)
              diag(this$P3) <- 1
              this$optimcontrol <- list(fnscale = -1, maxit = 1500, reltol = 1e-5)
              return(this)
            }

            return(this)
          })


## -- .loglik.stcc.1() ####
setGeneric(name=".loglik.stcc.1",
           function(optimpars="numeric",z="matrix",stcc="corr_class")standardGeneric(".loglik.stcc.1"),
           valueClass = "numeric",
           signature = c("optimpars","z","stcc")
)

setMethod(".loglik.stcc.1",signature = c("numeric","matrix","corr_class"),
          function(optimpars,z,stcc){

              # input: optimpars   -- c(speed,loc1,[loc2])
              #        z           -- volatility standardised returns (matrix TxN)
              #        stcc        -- object containing all the other parameters
            if(stcc$type != CORRtype$STCC.1) warning("Wrong STCC type")
            
              err_output <- -1e10  
            
              STCC <- stcc
              st <- STCC$st
              shape <- STCC$shape
              speedoption <- STCC$speedopt
              Tobs <- STCC@Tobs
              

              #### ======== constraint checks ======== ####
              
              # # Check 1: Confirm we have a valid shape & extract speed & location:
              
              if(shape==CORRshape$double) numTRpars <- 3 else numTRpars <- 2
              TRpars <- tail(optimpars,numTRpars)
              speed <- TRpars[1]
              loc1 <- TRpars[2]
              if(numTRpars==3) loc2 <- TRpars[3] else loc2 <- NA
              
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

              vP1 <- tmp.par[1:STCC@nr.covPars]
              mP <- myUnVecl(vP1)
              tmp.par <- tail(tmp.par,-STCC@nr.covPars)
              eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
              # Check for SPD - positive-definite check:
              if (min(eig$values) <= 0) return(err_output)
              
              vP2 <- tmp.par[1:STCC@nr.covPars]
              mP <- myUnVecl(vP2)
              tmp.par <- tail(tmp.par,-STCC@nr.covPars)
              eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
              # Check for SPD - positive-definite check:
              if (min(eig$values) <= 0) return(err_output)
              
              st_c <- 0
              if(shape == CORRshape$single) { st_c <- st-loc1 }
              if(shape == CORRshape$double) { st_c <- (st-loc1)*(st-loc2) }
              if(shape == CORRshape$double1loc) { st_c <- (st-loc1)^2 }
              
              G <- 0
              if(speedoption == CORRspeedopt$gamma) { G <- 1/(1+exp(-speed*st_c)) }
              if(speedoption == CORRspeedopt$gamma_std) { G <- 1/(1+exp(-speed/sd(st)*st_c)) }
              if(speedoption == CORRspeedopt$eta) { G <- 1/(1+exp(-exp(speed)*st_c)) }

              Gt <- matrix(G,nrow = Tobs,ncol = 1)
              Pt <- t(apply(Gt,MARGIN = 1,FUN = function(X,P1,P2) ((1-X)*P1 + X*P2), P1=vP1, P2=vP2))
              
              llt <- vector("numeric")
              for(t in 1:Tobs) {
                mPt <- myUnVecl(Pt[t,])
                llt[t] <- -0.5*log(det(mPt)) -0.5*( t(z[t,])%*%(qr.solve(mPt))%*%z[t,])
              } 
              # Return: 
              return(sum(llt))  
               
          })


## -- .loglik.stcc.2() ####
setGeneric(name=".loglik.stcc.2",
           function(optimpars="numeric",z="matrix",stcc="corr_class")standardGeneric(".loglik.stcc.2"),
           valueClass = "numeric",
           signature = c("optimpars","z","stcc")
)

setMethod(".loglik.stcc.2",signature = c("numeric","matrix","corr_class"),
          function(optimpars,z,stcc){
            
            # input: optimpars   -- c(speed,loc1,[loc2])
            #        z           -- volatility standardised returns (matrix TxN)
            #        stcc        -- object containing all the other parameters
            if(stcc$type != CORRtype$STCC.2) warning("Wrong STCC type")
            
            err_output <- -1e10  
            
            STCC <- stcc
            st <- STCC$st
            shape <- STCC$shape
            speedoption <- STCC$speedopt
            Tobs <- STCC@Tobs
            
            
            #### ======== constraint checks ======== ####
            
            # # Check 1: Confirm we have a valid shape & extract speed & location:
            
            if(shape[1]==CORRshape$double) numTRpars1 <- 3 else numTRpars1 <- 2
            if(shape[2]==CORRshape$double) numTRpars2 <- 3 else numTRpars2 <- 2
            
            tmp.par <- optimpars
            TRpars <- tail(tmp.par,numTRpars2)
            speed2 <- TRpars[1]
            loc2.1 <- TRpars[2]
            if(numTRpars2 == 3) loc2.2 <- TRpars[3] else loc2.2 <- NA
            
            tmp.par <- tail(tmp.par,-numTRpars2)
            TRpars <- tail(tmp.par,numTRpars1)
            speed1 <- TRpars[1]
            loc1.1 <- TRpars[2]
            if(numTRpars1== 3) loc1.2 <- TRpars[3] else loc1.2 <- NA
            
            # # Check 2: Check the boundary values for speed params:
            # #speedoption -- 1=gamma, 2=gamma/std(st), 3=exp(eta)
            maxSpeed <- switch(speedoption,1000,(1000/sd(st)),7.0,0.30)
            if (speed1 >= maxSpeed) return(err_output)
            if (speed1 <= 0) return(err_output)
            if (speed2 >= maxSpeed) return(err_output)
            if (speed2 <= 0) return(err_output)
            # 
            # # Check 3: Check the locations fall within min-max values of st
            # # We must have at least 1 loc1 to be inside this shape..loop, so no need to check if loc1 contains a valid value:
            if (loc2.1 < loc1.1) return(err_output)
            if (loc1.1 <= 0) return(err_output)
            if (loc2.1 <= 0) return(err_output)
            
            #### ======== calculate loglikelihood ======== ####
            
            tmp.par <- optimpars
            
            vP1 <- tmp.par[1:STCC@nr.covPars]
            mP <- myUnVecl(vP1)
            tmp.par <- tail(tmp.par,-STCC@nr.covPars)
            eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
            # Check for SPD - positive-definite check:
            if (min(eig$values) <= 0) return(err_output)
            
            vP2 <- tmp.par[1:STCC@nr.covPars]
            mP <- myUnVecl(vP2)
            tmp.par <- tail(tmp.par,-STCC@nr.covPars)
            eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
            # Check for SPD - positive-definite check:
            if (min(eig$values) <= 0) return(err_output)
            
            vP3 <- tmp.par[1:STCC@nr.covPars]
            mP <- myUnVecl(vP3)
            tmp.par <- tail(tmp.par,-STCC@nr.covPars)
            eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
            # Check for SPD - positive-definite check:
            if (min(eig$values) <= 0) return(err_output)
            
            st_c <- 0
            if(shape[1] == CORRshape$single) { st_c <- st-loc1.1 }
            if(shape[1] == CORRshape$double) { st_c <- (st-loc1.1)*(st-loc1.2) }
            if(shape[1] == CORRshape$double1loc) { st_c <- (st-loc1.1)^2 }
            
            G1 <- 0
            if(speedoption == CORRspeedopt$gamma) { G1 <- 1/(1+exp(-speed1*st_c)) }
            if(speedoption == CORRspeedopt$gamma_std) { G1 <- 1/(1+exp(-speed1/sd(st)*st_c)) }
            if(speedoption == CORRspeedopt$eta) { G1 <- 1/(1+exp(-exp(speed1)*st_c)) }
            
            Gt1 <- matrix(G1,nrow = Tobs,ncol = 1)
            
            st_c <- 0
            if(shape[2] == CORRshape$single) { st_c <- st-loc2.1 }
            if(shape[2] == CORRshape$double) { st_c <- (st-loc2.1)*(st-loc2.2) }
            if(shape[2] == CORRshape$double1loc) { st_c <- (st-loc2.1)^2 }
            
            G2 <- 0
            if(speedoption == CORRspeedopt$gamma) { G2 <- 1/(1+exp(-speed2*st_c)) }
            if(speedoption == CORRspeedopt$gamma_std) { G2 <- 1/(1+exp(-speed2/sd(st)*st_c)) }
            if(speedoption == CORRspeedopt$eta) { G2 <- 1/(1+exp(-exp(speed2)*st_c)) }
            
            Gt2 <- matrix(G2,nrow = Tobs,ncol = 1)

            Pt1 <- t(apply(Gt1,MARGIN = 1,FUN = function(X,P1,P2) ((1-X)*P1 + X*P2), P1=vP1, P2=vP2))
            Pt2 <- (1-G2)*Pt1
            Pt3 <- t(apply(Gt2,MARGIN = 1,FUN = function(X,P3) (X*P3), P3=vP3) )
            
            Pt <- Pt2 + Pt3
            
            llt <- vector("numeric")
            for(t in 1:Tobs) {
              mPt <- myUnVecl(Pt[t,])
              llt[t] <- -0.5*log(det(mPt)) -0.5*( t(z[t,])%*%(qr.solve(mPt))%*%z[t,])
            } 
            # Return: 
            return(sum(llt))  
            
          })


## === .calculate_Pt.1() === ####            
setGeneric(name=".calculate_Pt.1",
           function(corr="corr_class")standardGeneric(".calculate_Pt.1"),
           valueClass = "matrix",
           signature = c("corr")
)

setMethod(".calculate_Pt.1",signature = c("corr_class"),
          function(corr){
            # -- Validation -- #
            if (corr$type==CORRtype$CCC) return(corr)
            
            st <- corr@st
            shape <- corr$shape
            speedoption <- corr$speedopt
            speed <- corr$Estimated$pars["speed"]
            loc1 <- corr$Estimated$pars["loc1"]
            loc2 <- corr$Estimated$pars["loc2"]
            
            vP1 <- myVecl(corr$Estimated$P1)
            vP2 <- myVecl(corr$Estimated$P2)
            
            st_c <- 0
            if(shape == CORRshape$single) { st_c <- st-loc1 }
            if(shape == CORRshape$double) { st_c <- (st-loc1)*(st-loc2) }
            if(shape == CORRshape$double1loc) { st_c <- (st-loc1)^2 }
            
            G <- 0
            if(speedoption == CORRspeedopt$gamma) { G <- 1/(1+exp(-speed*st_c)) }
            if(speedoption == CORRspeedopt$gamma_std) { G <- 1/(1+exp(-speed/sd(st)*st_c)) }
            if(speedoption == CORRspeedopt$eta) { G <- 1/(1+exp(-exp(speed)*st_c)) }
            
            Gt <- matrix(G,nrow = corr@Tobs,ncol = 1)
            Pt <- t(apply(Gt,MARGIN = 1,FUN = function(X,P1,P2) ((1-X)*P1 + X*P2), P1=vP1, P2=vP2))
            
            return(Pt)

          })

## === .calculate_Pt.2() === ####            
setGeneric(name=".calculate_Pt.2",
           function(corr="corr_class")standardGeneric(".calculate_Pt.2"),
           valueClass = "matrix",
           signature = c("corr")
)

setMethod(".calculate_Pt.2",signature = c("corr_class"),
          function(corr){
            # -- Validation -- #
            if (corr$type==CORRtype$CCC) return(corr)
            
            Tobs <- corr@Tobs
            st <- corr@st
            shape <- corr$shape
            speedoption <- corr$speedopt
            
            speed1 <- corr$Estimated$pars1["speed"]
            speed2 <- corr$Estimated$pars2["speed"]
            loc1.1 <- corr$Estimated$pars1["loc1"]
            loc1.2 <- corr$Estimated$pars1["loc2"]
            loc2.1 <- corr$Estimated$pars2["loc1"]
            loc2.2 <- corr$Estimated$pars2["loc2"]
            
            vP1 <- myVecl(corr$Estimated$P1)
            vP2 <- myVecl(corr$Estimated$P2)
            vP3 <- myVecl(corr$Estimated$P3)
            
            st_c <- 0
            if(shape[1] == CORRshape$single) { st_c <- st-loc1.1 }
            if(shape[1] == CORRshape$double) { st_c <- (st-loc1.1)*(st-loc1.2) }
            if(shape[1] == CORRshape$double1loc) { st_c <- (st-loc1.1)^2 }
            
            G1 <- 0
            if(speedoption == CORRspeedopt$gamma) { G1 <- 1/(1+exp(-speed1*st_c)) }
            if(speedoption == CORRspeedopt$gamma_std) { G1 <- 1/(1+exp(-speed1/sd(st)*st_c)) }
            if(speedoption == CORRspeedopt$eta) { G1 <- 1/(1+exp(-exp(speed1)*st_c)) }

            Gt1 <- matrix(G1,nrow = corr@Tobs,ncol = 1)
            
            st_c <- 0
            if(shape[2] == CORRshape$single) { st_c <- st-loc2.1 }
            if(shape[2] == CORRshape$double) { st_c <- (st-loc2.1)*(st-loc2.2) }
            if(shape[2] == CORRshape$double1loc) { st_c <- (st-loc2.1)^2 }
            
            G2 <- 0
            if(speedoption == CORRspeedopt$gamma) { G2 <- 1/(1+exp(-speed2*st_c)) }
            if(speedoption == CORRspeedopt$gamma_std) { G2 <- 1/(1+exp(-speed2/sd(st)*st_c)) }
            if(speedoption == CORRspeedopt$eta) { G2 <- 1/(1+exp(-exp(speed2)*st_c)) }
            
            Pt1 <- t(apply(Gt1,MARGIN = 1,FUN = function(X,P1,P2) ((1-X)*P1 + X*P2), P1=vP1, P2=vP2))
            Pt2 <- (1-G2)*Pt1
            Pt3 <- matrix(G2*vP3,nrow = Tobs, ncol = 6)
            Pt <- Pt2 + Pt3
            
            return(Pt)
            
          })

## --- Public Methods --- ####

EstimateSTCC.1 <- function(z,stccObj,calcSE=FALSE,verbose=FALSE) {
  
  STCC <- stccObj
  STCC$Estimated <- list()

  optimpars <- myVecl(STCC$P1) 
  numCovPars <- length(optimpars)
  optimpars <- c(optimpars, myVecl(STCC$P2) )
  optimpars <- na.omit( c(optimpars, STCC$pars["speed"], STCC$pars["loc1"], STCC$pars["loc2"] ) )
  
  # Set some variables that we will need later
  STCC@nr.covPars <- as.integer(numCovPars)
  STCC@Tobs <- as.integer(NROW(z))
  
  ### ---  Call optim to calculate the estimate --- ###
  if (verbose) STCC$optimcontrol$trace <- 10
  
  tmp <- NULL
  try(tmp <- optim(optimpars,.loglik.stcc.1,z,STCC,gr=NULL,method="BFGS",control=STCC$optimcontrol,hessian=calcSE))
  
  ### ---  Interpret the response from optim --- ###
  # An unhandled error could result in a NULL being returned by optim()
  if (is.null(tmp)) {
    STCC$Estimated$value <- err_output 
    STCC$Estimated$error <- TRUE
    return(STCC)
  }
  
  #Optim converged successfully => we expect tmp$par to have good estimates!
  if (tmp$convergence==0) {
    STCC$Estimated$error <- FALSE
    STCC$Estimated$value <- tmp$value
    
    tmp.par <- tmp$par
    STCC$Estimated$P1 <- myUnVecl(tmp.par[1:numCovPars])
    tmp.par <- tail(tmp.par,-numCovPars)
    STCC$Estimated$P2 <- myUnVecl(tmp.par[1:numCovPars])
    tmp.par <- tail(tmp.par,-numCovPars)
    STCC$Estimated$pars <- tmp.par
    if(STCC$shape != CORRshape$double) STCC$Estimated$pars <- c(STCC$Estimated$pars,NA)
    names(STCC$Estimated$pars) <- names(STCC$pars)

    if (calcSE) {
      STCC$Estimated$hessian <- tmp$hessian
      vecSE <- vector("numeric")
      try(vecSE <- sqrt(-diag(qr.solve(tmp$hessian))))
      
      if(length(vecSE) > 0) {
        STCC$Estimated$P1.se <- myUnVecl(vecSE[1:numCovPars])
        vecSE <- tail(vecSE,-numCovPars)
        STCC$Estimated$P2.se <- myUnVecl(vecSE[1:numCovPars])
        vecSE <- tail(vecSE,-numCovPars)
        
        STCC$Estimated$pars.se <- vecSE
        if(STCC$shape != CORRshape$double) STCC$Estimated$pars.se <- c(STCC$Estimated$pars.se,NA)
        names(STCC$Estimated$pars.se) <- names(STCC$pars)  
      }
    }
    STCC$Estimated$Pt <- .calculate_Pt.1(STCC)
    
  } else { 
    #Failed to converge
    STCC$Estimated$error <- TRUE
    STCC$Estimated$value <- err_output
    STCC$Estimated$optimoutput <- tmp
  }
  
  if (verbose) STCC$Estimated$optimoutput <- tmp
  #Return:
  STCC
  
}  #End: EstimateSTCC.1()


EstimateSTCC.2 <- function(z,stccObj,calcSE=TRUE,verbose=TRUE) {
  
  STCC <- stccObj
  STCC$Estimated <- list()
  
  optimpars <- myVecl(STCC$P1) 
  numCovPars <- length(optimpars)
  optimpars <- c(optimpars, myVecl(STCC$P2) )
  optimpars <- c(optimpars, myVecl(STCC$P3) )
  optimpars <- na.omit( c(optimpars, STCC$pars1["speed"], STCC$pars1["loc1"], STCC$pars1["loc2"] ) )
  optimpars <- na.omit( c(optimpars, STCC$pars2["speed"], STCC$pars2["loc1"], STCC$pars2["loc2"] ) )
  
  # Set some variables that we will need later
  STCC@nr.covPars <- as.integer(numCovPars)
  STCC@Tobs <- as.integer(NROW(z))
  
  ### ---  Call optim to calculate the estimate --- ###
  if (verbose) STCC$optimcontrol$trace <- 10
  
  tmp <- NULL
  try(tmp <- optim(optimpars,.loglik.stcc.2,z,STCC,gr=NULL,method="BFGS",control=STCC$optimcontrol,hessian=calcSE))
  
  ### ---  Interpret the response from optim --- ###
  # An unhandled error could result in a NULL being returned by optim()
  if (is.null(tmp)) {
    STCC$Estimated$value <- err_output 
    STCC$Estimated$error <- TRUE
    return(STCC)
  }
  
  #Optim converged successfully => we expect tmp$par to have good estimates!
  if (tmp$convergence==0) {
    STCC$Estimated$error <- FALSE
    STCC$Estimated$value <- tmp$value
    
    tmp.par <- tmp$par
    STCC$Estimated$P1 <- myUnVecl(tmp.par[1:numCovPars])
    tmp.par <- tail(tmp.par,-numCovPars)
    STCC$Estimated$P2 <- myUnVecl(tmp.par[1:numCovPars])
    tmp.par <- tail(tmp.par,-numCovPars)
    STCC$Estimated$P3 <- myUnVecl(tmp.par[1:numCovPars])
    tmp.par <- tail(tmp.par,-numCovPars)
    
    nrTRpars <- length(na.omit(STCC$pars2))
    STCC$Estimated$pars2 <- tail(tmp.par,nrTRpars)
    if(STCC$shape[2] != CORRshape$double) STCC$Estimated$pars2 <- c(STCC$Estimated$pars2,NA)
    tmp.par <- tail(tmp.par,-nrTRpars)

    STCC$Estimated$pars1 <- tmp.par
    if(STCC$shape[1] != CORRshape$double) STCC$Estimated$pars1 <- c(STCC$Estimated$pars1,NA)
    
    if (calcSE) {
      STCC$Estimated$hessian <- tmp$hessian
      vecSE <- vector("numeric")
      try(vecSE <- sqrt(-diag(qr.solve(tmp$hessian))))
      
      if(length(vecSE) > 0) {
        STCC$Estimated$P1.se <- myUnVecl(vecSE[1:numCovPars])
        vecSE <- tail(vecSE,-numCovPars)
        STCC$Estimated$P2.se <- myUnVecl(vecSE[1:numCovPars])
        vecSE <- tail(vecSE,-numCovPars)
        STCC$Estimated$P3.se <- myUnVecl(vecSE[1:numCovPars])
        vecSE <- tail(vecSE,-numCovPars)
        
        nrTRpars <- na.omit(length(STCC$pars2))
        STCC$Estimated$pars2.se <- tail(vecSE,nrTRpars)
        if(STCC$shape[2] != CORRshape$double) STCC$Estimated$pars2.se <- c(STCC$Estimated$pars2.se,NA)
        vecSE <- tail(vecSE,-nrTRpars)
        
        STCC$Estimated$pars1.se <- vecSE
        if(STCC$shape[1] != CORRshape$double) STCC$Estimated$pars1.se <- c(STCC$Estimated$pars1.se,NA)
      }
    }
    #STCC$Estimated$Pt <- .calculate_Pt.2(STCC)
    
  } else { 
    #Failed to converge
    STCC$Estimated$error <- TRUE
    STCC$Estimated$value <- err_output
    STCC$Estimated$optimoutput <- tmp
  }
  
  if (verbose) STCC$Estimated$optimoutput <- tmp
  #Return:
  STCC
  
}  #End: EstimateSTCC.2()


## --- PRIVATE METHODS --- ####
# myLogLik.stcc <- function(optimpars,z,stcc,return_ll=TRUE){
#   # model: STCC
#   # input: pars        -- c(speed,loc1,[loc2])
#   #        z           -- volatility standardised returns (matrix TxN)
#   #        stcc        -- list containing all the other parameters
#   
#   STCC <- stcc
#   speedoption <- STCC$speedoption
#   shape <- STCC$shape  
#   st <- STCC$st
#   
#   Tobs <- NROW(z)
#   N <- NCOL(z)
#   
#   #### ======== constraint checks ======== ####
#   
#   # Check 1: Confirm we have a valid shape & extract speed & location:

#   if(shape==CORRshape$double) numTRpars <- 3 else numTRpars <- 2
#   TRpars <- tail(optimpars,numTRpars)
#   speed <- TRpars[1]
#   loc1 <- TRpars[2]
#   if(numTRpars==3) loc2 <- TRpars[3] else loc2 <- NA
#   
#   # Check 2: Check the boundary values for speed params:
#   #speedoption -- 1=gamma, 2=gamma/std(st), 3=exp(eta), 4=1/lambda^2
#   maxSpeed <- switch(speedoption,1000,(1000/sd(st)),7.0,0.30)
#   if (speed > maxSpeed) return(err_output)
#   if (speed < 0) return(err_output)
#   
#   # Check 3: Check the locations fall within min-max values of st
#   # We must have at least 1 loc1 to be inside this shape..loop, so no need to check if loc1 contains a valid value:
#   if (loc1 < min(st)) return(err_output)
#   if (loc1 > max(st)) return(err_output)
#   if (numTRpars==3) {
#     if (loc2 < min(st)) return(err_output)
#     if (loc2 > max(st)) return(err_output)  
#   }
#   
#   
#   #### ======== calculate loglikelihood ======== ####
#   # - - - CCC - - -
#   if (STCC$type==CORRtype$CCC){
#     vP <- c(0,0,0) #pars
#     mP <- myUnVecl(vP)
#     eig <- eigen(mP,symmetric=TRUE,only.values = TRUE)
#     if (min(eig$values) <= 0) return(err_output)
#     
#     # - - - P(t) and loglik-value
#     llt <- rep(0,Tobs)
#     mPinv <- solve(mP)
#     const <- -0.5*log(det(mP))
#     for(t in seq(1,Tobs)) llt[t] <- const - 0.5*(t(z[t,])%*%(mPinv)%*%z[t,])
#   }
#   
#   # - - - STCC - - -
#   if (STCC$type==CORRtype$STCC) {
#     numCovPars <- NROW(myVecl(STCC$P1))
#     vP1 <- optimpars[1:numCovPars]
#     vP2 <- optimpars[(numCovPars+1):(2*numCovPars)]
#     
#     mP1 <- myUnVecl(vP1)
#     eig1 <- eigen(mP1,symmetric=TRUE,only.values=TRUE)
#     # Check for SPD - positive-definite check:
#     if (min(eig1$values) <= 0) return(err_output)
#     mP2 <- myUnVecl(vP2)
#     eig2 <- eigen(mP2,symmetric=TRUE,only.values=TRUE)
#     # Check for SPD - positive-definite check:
#     if (min(eig2$values) <= 0) return(err_output)
#     
#     # - - - Calculate G(t) - - -
#     # if(shape == CORRshape$single){st_c <- st-loc1}
#     # if(shape == CORRshape$double){st_c <- (st-loc1)*(st-loc2)}
#     # if(shape == CORRshape$double1loc){st_c <- (st-loc1)^2}
#     st_c <- switch(shape,st-loc1,(st-loc1)*(st-loc2),(st-loc1)^2)
#     G_inv <- switch(speedoption, 1+exp(-speed*st_c), 1+exp(-speed/sd(st)*st_c), 1+exp(-exp(speed)*st_c))
#     Gt <- matrix(1/G_inv,nrow = Tobs,ncol = 1)
#     
#     # - - - P(t) and loglik-value
#     #Pt <- matrix(0,nrow=Tobs,ncol=(N*(N-1)/2))
#     llt <- NULL
#     #calcPt <- function(X,P1,P2) ((1-X)*P1 + X*P2)
#     #Pt <- t(apply(Gt,MARGIN = 1,FUN = calcPt, P1=vP1, P2=vP2))
#     Pt <- t(apply(Gt,MARGIN = 1,FUN = function(X,P1,P2) ((1-X)*P1 + X*P2), P1=vP1, P2=vP2))
#     for(t in seq(1,Tobs))
#     {
#       #Pt[t,] <- (1-Gt[t])*vP1 + Gt[t]*vP2  #This line has been replaced by the apply() function above
#       mPt <- myUnVecl(Pt[t,])
#       mPtinv <- chol2inv(chol(mPt))
#       llt[t] <- -0.5*log(det(mPt)) -0.5*( t(z[t,])%*%(mPtinv)%*%z[t,])
#     } # End: for loop
#     
#   } # End: if(STCC$type==CORRtype$STCC)
#   
#   # Return: 
#   if (return_ll) return(sum(llt)) else {
#     if (STCC$type==CORRtype$CCC) Pt <- matrix(vP,nrow=Tobs,ncol=length(vP),byrow=TRUE)
#     return(Pt)
#   }
#   
# }  #End: myLogLik.stcc()

calcStderr_STCC <- function(e,stccObj) {
  STCC <- stccObj
  STCC$Estimated$stderr <- NULL
  
  if(is.null(STCC$Estimated$hessian)) {
    warning("This method of generating Standard Errors is unreliable. \nPlease re-estimate using the 'calcSE=TRUE' parameter")
    optimpars <- NULL
    optimpars <- c(myVecl(STCC$Estimated$P1),myVecl(STCC$Estimated$P2),STCC$Estimated$TRpars)
    STCC$Estimated$hessian <- optimHess(optimpars,myLogLik.stcc,gr=NULL,e,STCC,control=STCC$optimcontrol)
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




## --- Override Methods --- ####

## -- plot() ####
setMethod("plot",signature = c(x="corr_class",y="missing"),
          function(x, y, ...){
            this <- x
            plot(this@h, type='l', ylab = "Cond.Variance", ...)
          })

