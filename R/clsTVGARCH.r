# *************************************************************************************************************** #
#                   Class definitions & methods for the garch_class, tv_class and tvgarch_class                ####
# *************************************************************************************************************** #

# ---  TV --- SECTION  ####

{## Help doco: tv() ####
  #' @noRd
#' @title
  #' A time-varying
  #'
  #' @description
  #' `tv` is a
  #'
  #' @usage tv(st,shape)
  #'
  #' @param st A smooth
  #' @param shape Use enum, tvshape$...
  #'
  #' @details
  #' This object
  #'
  #' ```
  #'   myTV = tv(st,tvshape$single)
  #' ```
  #'
  #'
  #' @returns A tv_class object.
  #'
  #' @note
  #' I am a note
  #'
  #'
}

{## TV_CLASS Definition ----

tv <- setClass(Class = "tv_class",
               slots = c(Tobs="integer",st="numeric",g="numeric",delta0free="logical",nr.pars="integer", nr.transitions="integer"),
               contains = c("namedList")
)

setMethod("initialize","tv_class",
          function(.Object,...){
            .Object <- callNextMethod(.Object,...)
            # Slots
            .Object@st <- c(NaN)
            .Object@g <- 1
            .Object@delta0free <- TRUE
            .Object@nr.pars <- as.integer(1)
            .Object@nr.transitions <- as.integer(0)
            .Object@Tobs <- as.integer(0)
            # Properties
            .Object$shape <- tvshape$delta0only
            .Object$speedopt <- speedopt$none
            .Object$delta0 <- 1.0
            .Object$pars <- matrix(NA,4,1)
            .Object$optimcontrol <- list(fnscale = -1, reltol = 1e-8)  # Maximiser with default tolerance

            # Return:
            .Object
          })
}

.setInitialPars = function(tvObj){
  this <- tvObj

  nrLoc <- this@nr.transitions + length(this$shape[this$shape==tvshape$double])
  locNum <- 1
  locDen <- nrLoc + 1
  pars <- NULL
  parNames <- NULL
  for(n in 1:this@nr.transitions){

    loc1 <- round(locNum/locDen,4)

    if(this$shape[n] == tvshape$double) {
      loc2 <- round((locNum+1)/locDen,4)
      locNum <- locNum + 2
    } else {
      loc2 <- NA
      locNum <- locNum + 1
    }
    pars <- c(pars,1,3,loc1,loc2)

  }
  this$pars <- matrix(pars,nrow=4)
  return(this)
}

{## Constructor: tv() ----
  setGeneric(name="tv",
             valueClass = "tv_class",
             signature = c("st","shape"),
             def = function(st,shape){

               # Validate shape:
               if(length(shape) > 1){
                 if(any(shape == tvshape$delta0only)) stop("Invalid shape: delta0only / 0 is not a valid transition shape")
               }

               this <- new("tv_class")
               this$shape <- shape
               this@st <- st
               this@Tobs <- length(st)
               this@g <- rep(this$delta0,this@Tobs)

               if(shape[1] == tvshape$delta0only){
                 this@nr.transitions <- as.integer(0)
                 this$optimcontrol$ndeps <- c(1e-3)
                 this$optimcontrol$parscale <- c(1)
               }else {
                 this$speedopt <- speedopt$eta
                 this@nr.transitions <- as.integer(length(shape))
                 # Create the starting Pars matrix
                 this  <- .setInitialPars(this)
                 rownames(this$pars) <- c("deltaN","speedN","locN1","locN2")
                 this@nr.pars <- as.integer(length(this$pars[!is.na(this$pars)]) + 1)  # +1 for delta0
                 this$optimcontrol$ndeps <- rep(1e-3,this@nr.pars)


                 #TODO: Improve the parScale to better manage different multiple location Options
                 parScale_exDelta0 <- c(rep(c(4,3,0.5),this@nr.transitions) )
                 # Tricky bit of 'maths' below to produce NA's in the NA locations
                 # TODO!!! Bugs in code below - FIX for when loc2 is not NA
                 parScale_exDelta0 <- parScale_exDelta0 + (as.vector(head(this$pars,-1)) - as.vector(head(this$pars,-1)))
                 this$optimcontrol$parscale <- c(1,parScale_exDelta0[!is.na(parScale_exDelta0)])

               }
               return(this)
             }
  )
}

{## plot(tv) override ----
  setMethod("plot",signature = c(x="tv_class",y="missing"),
            function(x, y,...){
              plot.default(x=sqrt(x@g), type='l', ylab = "sqrt(g)", ...)
            })
}

{## summary(tv) override ----
  setMethod("summary",signature="tv_class",
            function(object,...){
              this <- object

              delta0Summary <- NULL
              parsSummary <- NULL

              if(!is.null(this$Estimated)){
                # Calculate significance indicator for delta0
                if(is.null(this$Estimated$delta0_se)) this$Estimated$delta0_se <- NaN
                se <- this$Estimated$delta0_se
                d0 <- this$Estimated$delta0
                d0Sig <- ""
                if(!is.nan(se)){
                  if(se*2.576 < abs(d0)) { d0Sig <- "***" }
                  else if(se*1.96 < abs(d0)) { d0Sig <- "** " }
                  else if(se*1.645 < abs(d0)) { d0Sig <- "*  " }
                }
                delta0Summary <- paste0(round(this$Estimated$delta0,6),"    se0 = ",round(this$Estimated$delta0_se,6),d0Sig)

                # Format tv parameters if present:
                parsSummary <- NULL
                if(this@nr.transitions > 0){
                  parsVec <-  round(as.vector(this$Estimated$pars),6)
                  if(!is.null(this$Estimated$se) ){
                    seVec <- round(as.vector(this$Estimated$se),6)
                    seVecSig <- vector("character", length(seVec))

                    for(n in seq_along(parsVec)){
                      if(is.nan(seVec[n])) {
                        seVecSig[n] <- "   "
                      } else {
                        # Calculate a significance indicator
                        if(seVec[n]*2.576 < abs(parsVec[n]) ) { (seVecSig[n] <- "***") }
                        else if(seVec[n]*1.96 < abs(parsVec[n]) ) { (seVecSig[n] <- "** ") }
                        else if(seVec[n]*1.645 < abs(parsVec[n]) ) { (seVecSig[n] <- "*  ") }
                      }
                    }
                  } else {
                    seVec <- rep(NaN,length(this$pars))
                    seVecSig <- rep("   ", length(seVec))
                  }

                  seMat <- matrix(seVec,nrow=4)
                  colnames(seMat) <- paste("se" ,1:this@nr.transitions,sep = "")
                  # Build parsSummary table and insert the significance indicators
                  parsSummary <- data.frame(NA,stringsAsFactors = FALSE)
                  for (n in 1:NCOL(this$Estimated$pars)){
                    sig <- matrix(seVecSig[1:4],nrow=4)
                    parsSummary <- cbind(parsSummary,round(this$Estimated$pars[,n,drop=F],6),seMat[,n,drop=F],sig)
                    seVecSig <- tail(seVecSig,-4)
                  }
                }
              }

              cat("\n\nTV OBJECT\n")
              cat("\nTransition Shapes: ", this$shape ,"\n")
              if(!is.null(this$Estimated)){
                cat("\nEstimation Results:\n")
                cat("\nDelta0 =",delta0Summary,"\n\n")
                if(this@nr.transitions > 0) print(parsSummary[,-1])
                cat("\nLog-likelihood value(TV): ",this$Estimated$value)
              }

            })
}

# ---  ESTIMATE TV --- SECTION  ####

{##  Help doco: estimateTV(e,tv,garch,ctrl) ####
  #' @noRd
#' @title
  #' Estimates a tv model
  #'
  #' @description
  #' `estimateTV` is a
  #'
  #' @usage estimateTV(e,tvObj,estimationControl)
  #'
  #' @param e An estimated TV model
  #' @param tvObj Use Enum, garchtype$...
  #' @param estimationControl A list
  #'
  #' @details
  #' This object
  #'
  #' ```
  #'   myTv = estimateTV(e,myTv,estimationControl)
  #' ```
  #'
  #'
  #' @returns A tv_class object.
  #'
  #' @note
  #' I am a note
  #'
  #'
}

{##  estimateTV(e,tv,garch,ctrl) ----
#  estimateTV <- function(e,tvObj,garchObj,estimationControl){0}
  .estimateTV <- function(e,tvObj,garchObj,estimationControl){

    this <- tvObj
    # debug:
    #this <- TVspec
    #garchObj <- GARCH
    #estimationControl <- estCtrl

    # If there is no existing this$Estimated$, then create one
    if(!("Estimated" %in% names(this))) { this$Estimated <- list() }

    # Set verbose tracing:
    if (isTRUE(estimationControl$verbose)) {
      this$optimcontrol$trace <- 10
      cat("\nEstimating TV object...\n")
    } else this$optimcontrol$trace <- 0

    # Check for the simple case of just delta0 provided, no TV$pars
    if(this@nr.transitions == 0){ return( .estimateTV_noPars(e,this) ) }

    # Get the Optimpars from the tvObj (and update the tvObj as needed)
    tvOptimpars <- .getOptimpars_FromTV(this,estimationControl)
    optimpars <- tvOptimpars$Optimpars
    this <- tvOptimpars$tvObj

    # TODO: Really want to handle noGarch -> Need to set h(t) on a noGarch model
    # Check we have a garchObj, if not create one (so optim() doesn't complain):
    if(is.null(garchObj)){ garchObj <- garch(garchtype$general)}
    # Check we have a valid h(t) vector (default = 1), if not create one (so optim() doesn't complain):
    if(length(garchObj@h) == 1) { garchObj@h <- rep(1,NROW(e)) }

    ## --- Now call optim(.loglik.tv.univar) --- ##
    tmp <- NULL
    try(tmp <- optim(optimpars,.loglik.tv.univar,gr=NULL,e,this,garchObj,method="BFGS",control=this$optimcontrol))

    ## --- Attach results of estimation to the object --- ##

    if(!(is.null(tmp))){
      # Add the optim output to the estimated model
      this$Estimated$optimoutput <- tmp

      # Add the final results from optim() back into the estimated model
      this <- .setEstimatedPars_TV(this,tmp)

      # Get the conditional variances
      this@g <- calculate_g(this)

      # Calculate the parameter standard errors, if requested
      if (isTRUE(estimationControl$calcSE)) { this <- .setStdErrors_TV(this,tmp,garchObj) }
    }else{
      # optim() returned NULL.  What to do?
      # Return the object passed in with the Error = TRUE
      this$Estimated$error <- TRUE
      warning("estimateTV() generated a NULL return from optim()")
    }

    return(this)
  }
}

{## loglik.tv.univar() ----
  .loglik.tv.univar = function(optimpars,e,tvObj,garchObj){
    this <- tvObj

    error <- -1e6

    # Copy the optimpars into a local tv_object
    if (isTRUE(this@delta0free)) {
      this$Estimated$delta0 <- optimpars[1]
      this$Estimated$pars <- .estimatedParsToMatrix(this,tail(optimpars,-1))
    } else{
      if(!("delta0" %in% names(this$Estimated))) this$Estimated$delta0 <- this$delta0
      this$Estimated$pars <- .estimatedParsToMatrix(this,optimpars)
    }

    # Do paramater boundary checks:
    # Check 1: Check that delta0 is positive
    if (this$Estimated$delta0 < 0) return(error)

    if (this@nr.transitions > 0){
      # We have some Tv$pars
      vecSpeed <- this$Estimated$pars["speedN",(1:this@nr.transitions)]
      vecLoc1 <- this$Estimated$pars["locN1",(1:this@nr.transitions)]
      vecLoc2 <- this$Estimated$pars["locN2",(1:this@nr.transitions)]

      # Check 2: Check the boundary values for speed params:
      #speedoptions: 1=gamma, 2=gamma/std(st), 3=exp(eta), 4=1/lambda^2
      maxSpeed <- switch(this$speedopt,1000,(1000/sd(this@st)),7.0,0.30)
      if (max(vecSpeed) > maxSpeed) return(error)
      if (min(vecSpeed) < 0) return(error)

      # Check 3: Check the loc1 locations fall within min-max values of st
      # We must have at least 1 loc1 to be inside this shape..loop, so no need to check if loc1 contains a valid value:
      if (min(vecLoc1) < min(this@st)) return(error)
      if (max(vecLoc1) > max(this@st)) return(error)

      # Check 4: Check that loc1.1 < loc1.2 .. locN.1 < locN.2 for all G(i)
      # Method: Subtract loc1_pos vector from loc2_pos vector and ensure it is positive:
      tmp <- vecLoc2 - vecLoc1
      # Note: tmp will contain NA wherever a loc2 element was NA - we can strip these out:
      if (sum(tmp < 0,na.rm = TRUE) > 0) return(error)

      # Check 5: Check the loc2 locations fall within min-max values of st
      # Confirm we have at least one valid numeric loc 2, before checking min & max:
      if (any(!is.na(vecLoc2))) {
        if (min(vecLoc2,na.rm = TRUE) < min(this@st)) return(error)
        if (max(vecLoc2,na.rm = TRUE) > max(this@st)) return(error)
      }

      # Check 6: Check that loc1.1 < loc2.1 where 2 locations exist... for all G(i)
      # We do need to have at least 2 locations for this error check
      if (NROW(vecLoc1) > 1) {
        v1 <- head(vecLoc1,-1)
        v2 <- tail(vecLoc1,-1)
        if (sum(v2-v1 < 0) > 0) return(error)
      }

    }# End: paramater boundary checks:


    ll <- tryCatch({
      g <- calculate_g(this)
      if (min(g,na.rm = TRUE) <= 0) return(error)
      h <- garchObj@h
      #Return the LogLiklihood value:
      loglik.tvgarch.univar(e,g,h)

    }, warning = function(w){
      message("loglik.tv.univar warning: ", w$message)
    },error = function(err){
      message("loglik.tv.univar() error: ", err$message)
      ll <- error
    },finally = {
      # Any cleanup reqd.?
    })

    return(ll)

  }
}

{## estimateTV( - overloads - ) ----

  #TODO: Use noGarch instead of general - After fixing noGarch issues, @h, etc.

  setGeneric("estimateTV", valueClass = "tv_class",
             function(e, tvObj, garchObj, estimationControl) {
               standardGeneric("estimateTV")
             }
  )

  # TV <- estimateTV(e,tvSpec,garchSpec,estCtrl)
  setMethod("estimateTV",
            signature = c(e="numeric",tvObj="tv_class",garchObj="garch_class",estimationControl="list"), #note: All lower case
            function(e,tvObj,garchObj,estimationControl){
              .estimateTV(e,tvObj,garchObj,estimationControl)
            }
  )


  # TV <- estimateTV(e,tvSpec)
  setMethod("estimateTV",
            signature = c(e="numeric",tvObj="tv_class",garchObj="missing",estimationControl="missing"),
            function(e,tvObj){
              garchObj <- garch(garchtype$general)
              estimationControl <- list(calcSE=FALSE, verbose=FALSE, maxIter=100, fixStartPars=FALSE, startparAdjust=10)
              .estimateTV(e,tvObj,garchObj,estimationControl)
            }
  )

  # TV <- estimateTV(e,tvSpec,estCtrl):  Note: estCtrl can be passed-in by position, utilising the garchObj par
  setMethod("estimateTV",
            signature = c(e="numeric",tvObj="tv_class",garchObj="list",estimationControl="missing"),  #note: We need to switch the garchObj to accept the estCtrl list()
            function(e,tvObj,garchObj){
              estimationControl <- garchObj
              garchMissing <- garch(garchtype$general)  # Was not provided by user, they put estControl in this position
              estimationControl <- list(calcSE=FALSE, verbose=FALSE, maxIter=100, fixStartPars=FALSE, startparAdjust=10)
              .estimateTV(e,tvObj,garchMissing,estimationControl)
            }
  )

  # TV <- estimateTV(e,tvSpec,garchSpec):
  setMethod("estimateTV",
            signature = c(e="numeric",tvObj="tv_class",garchObj="garch_class",estimationControl="missing"),  #note: We need to switch the garchObj to accept the estCtrl list()
            function(e,tvObj,garchObj){
              estimationControl <- list(calcSE=FALSE, verbose=FALSE, maxIter=100, fixStartPars=FALSE, startparAdjust=10)
              .estimateTV(e,tvObj,garchObj,estimationControl)
            }
  )


}

{## estimateTV() - helpers ----
.estimateTV_noPars <- function(e,tvObj){

  this <- tvObj

  if(this@delta0free){

    this@nr.pars <- as.integer(1)
    this$Estimated$delta0 <- var(e)
  } else {

    # Estimate TV (no transitions AND delta0free = FALSE)
    this@nr.pars <- as.integer(0)
    # If there is no existing estimated this$Estimated$delta0, then use the starting param
    if(!("Estimated" %in% names(this))) {
      this$Estimated <- list()
      this$Estimated$delta0 < this$delta0}
  }

  # Now set g(t) based on the correct delta0
  this@g <- rep(this$Estimated$delta0,this@Tobs)

  this$Estimated$value <- sum(-0.5*log(2*pi) - 0.5*log(this@g) - (0.5*e^2)/this@g)
  this$Estimated$error <- FALSE

  this$Estimated$pars <- c(NA,NA,NA,NA)
  if(isTRUE(estimationControl$calcSE)) { this$Estimated$delta0_se <- NaN }

  # RETURN:
  return(this)
}

.getOptimpars_FromTV <- function(tvObj,estimationControl){

  # Derive the optimpars vector, allowing for NA's and delta0free:
  # Return both the optimpars vector and the tvObj, with corrected internals, nr.pars, optimControl list.

  this <- tvObj
  estCtrl <- estimationControl

  # Step1. Check for iterative-estimator controls, Set defaults if not exists:
  if(!("fixStartPars" %in% names(estCtrl))) { estCtrl$fixStartPars <- FALSE }
  if(!("startParAdjust" %in% names(estCtrl))) { estCtrl$startParAdjust <- 50 }

  # Get the starting pars (parsVec):
  if(isTRUE(estCtrl$fixStartPars) || (!("pars" %in% names(this$Estimated))) ){
    # Get starting pars from model spec if it hasn't been estimated yet, OR fixStartPars is TRUE
    parsVec <- c(this$delta0, as.vector(this$pars))
  }else{
    # Get last Estimates:
    parsVec <- c(this$Estimated$delta0, as.vector(this$Estimated$pars))  #TODO: Confirm this is correct for TVGARCH estimator
  }

  # Step2. Remove any NA's (loc2 params) - Note: This is CRITICAL as length(parsVec is used below)
  parsVec <- parsVec[!is.na(parsVec)]

  # Step3. When delta0 is NOT a free param, remove it from the parsVec and optimcontrol's
  if(isTRUE(this@delta0free)){
    this@nr.pars <- as.integer(length(parsVec))
    # TODO: BUG here if user switches delta0free ON, then OFF, then ON again the optimcontrol gets whacked!
    if(length(this$optimcontrol$ndeps) < length(parsVec)) {this$optimcontrol$ndeps <- c(1e-3,this$optimcontrol$ndeps)}
    if(length(this$optimcontrol$parscale) < length(parsVec)) {this$optimcontrol$parscale <- c(1.0,this$optimcontrol$parscale)}

  }else{
    # Drop delta0 from the parsVec:
    parsVec <- tail(parsVec,-1)
    this@nr.pars <- as.integer(length(parsVec))
    # TODO: BUG here if user switches delta0free On, then OFF, then ON again the optimcontrol gets whacked!
    if(length(this$optimcontrol$ndeps) != length(parsVec)) this$optimcontrol$ndeps <- tail(this$optimcontrol$ndeps,this@nr.pars)
    if(length(this$optimcontrol$parscale) != length(parsVec)) this$optimcontrol$parscale <- tail(this$optimcontrol$parscale,this@nr.pars)
  }

  # Shake the starting pars (parsVec):
  if(isFALSE(estCtrl$fixStartPars) ){
    # "shake": Move a few steps away:
    boundSize <- estCtrl$startParAdjust *  this$optimcontrol$ndeps    # Should be in range 10 - 100
    boundLo <- parsVec - boundSize
    boundHi <- parsVec + boundSize

    # # Anna's method
    # boundSize <- estCtrl$startParAdjust  # Should be small, in the order of 1e-2
    # boundLo <- parsVec/this$optimcontrol$parscale - boundSize
    # boundHi <- parsVec/this$optimcontrol$parscale + boundSize
    # boundLo <- boundLo * this$optimcontrol$parscale
    # boundHi <- boundHi * this$optimcontrol$parscale

    parsVec <- runif(this@nr.pars,boundLo,boundHi)
  }


  #RETURN:
  rtnList <- list()
  rtnList$Optimpars <- parsVec
  rtnList$tvObj <- this
  return(rtnList)

}

.setEstimatedPars_TV <- function(tvObj,optimTmp){

  this <- tvObj
  tmp <- optimTmp

  # An unhandled error could result in a NULL being returned by optim()
  if (is.null(tmp)) {
    this$Estimated$value <- NA
    this$Estimated$error <- TRUE
    warning("estimateTV() - optim failed and returned NULL. Check the optim controls & starting params")
    return(this)
  }
  if (tmp$convergence != 0) {
    this$Estimated$value <- NA
    this$Estimated$error <- TRUE
    this$Estimated$optimoutput <- tmp
    warning("estimateTV() - failed to converge. Check the optim controls & starting params")
    return(this)
  }

  # No optim issues, so set output
  this$Estimated$value <- tmp$value
  this$Estimated$error <- FALSE

  #Update the TV object parameters using optimised pars:
  if (isTRUE(this@delta0free)){
    # delta0 is the first optimised param
    this$Estimated$delta0 <- as.numeric(tmp$par[1])
    # All other params are in a matrix, with a column per transition
    this$Estimated$pars <- .estimatedParsToMatrix(this,tail(tmp$par,-1))
  } else{
    # delta0 is not an estimated param, use first round estimate. TODO: Confirm this logic:
    # Should always be there - or something is badly wrong
    # if(!("delta0" %in% names(this$Estimated))) this$Estimated$delta0 <- this$delta0    #Use the starting param if no estimated value exists
    this$Estimated$pars <- .estimatedParsToMatrix(this,tmp$par)
  }
  colnames(this$Estimated$pars) <- paste("st" ,1:this@nr.transitions,sep = "")

  #RETURN:
  return(this)
}

.setStdErrors_TV <- function(tvObj,optimTmp,garchObj){

  this <- tvObj
  tmp <- optimTmp

  cat("\nCalculating TV standard errors...\n")
  this$Estimated$hessian <- NULL
  this$Estimated$se <- NULL
  stdErrors <- NULL

  # Get the hessian from the optimiser:
  try(this$Estimated$hessian <- optimHess(tmp$par,.loglik.tv.univar,gr=NULL,e,this,garchObj,control=this$optimcontrol))

  # Handle optimHess errors
  try(stdErrors <- sqrt(-diag(invertHess(this$Estimated$hessian))))
  if(!is.null(stdErrors)){
    parsVec <-  as.vector(this$pars)

    if (isTRUE(this@delta0free)){
      this$Estimated$delta0_se <- stdErrors[1]
      stdErrors <- tail(stdErrors,-1)
    } else {
      this$Estimated$delta0_se <- NaN
      stdErrors <- tail(stdErrors)
    }

    seIndex <- 1
    for(n in seq_along(parsVec)){
      if(!is.na(parsVec[n])) {
        this$Estimated$se[n] <- stdErrors[seIndex]
        seIndex <- seIndex + 1
      } else this$Estimated$se[n] <- NaN
    }
    this$Estimated$se <- matrix(this$Estimated$se,nrow = 4)
    colnames(this$Estimated$se) <- paste("se" ,1:this@nr.transitions,sep = "")
  }

  return(this)

}

.estimatedParsToMatrix <- function(tvObj,optimpars){

  this <- tvObj

  # TODO: Do we want a hard 'stop' here?  Or a warning and soft return, e.g. pars=NULL?
  # I don't think we can ever get here with deltaZero_Only - remove it?
  if(this@nr.transitions == 0) stop("There are no parameters on this tv object")

  # The optimpars passed-in had all NA's stripped out, so to make our TV Object whole again we need to
  # Add NA's for all missing loc2 pars:
  naPars <- NULL
  for (i in seq_along(this$shape)) {

    if (this$shape[i] == tvshape$double) {
      naPars <- c(naPars,optimpars[1:4])
      optimpars <- optimpars[-(1:4)]
    } else {
      naPars <- c(naPars,optimpars[1:3],NA)
      optimpars <- optimpars[-(1:3)]
    }
  }
  # RETURN:
  matrix(naPars,nrow=4,ncol=NROW(this$shape),dimnames=list(c("deltaN","speedN","locN1","locN2"),NULL))

}

calculate_g = function(tvObj){
  this <- tvObj

  # 1. Initialise g(t) to a constant variance = delta0
  if(!("delta0" %in% names(this$Estimated))){
    # Set defaults if the TV object has not been estimated yet
    g <- rep(this$delta0,this@Tobs)
    this$Estimated$pars <- this$pars
  }else {
    g <- rep(this$Estimated$delta0,this@Tobs)
  }

  # 2. Update based on any transition parameters in the model
  if (this@nr.transitions > 0){
    st_c <- 0
    Gi <- 0
    # calulate 'g'
    for (i in 1:this@nr.transitions) {
      if(this$shape[i] == tvshape$single) st_c <- this@st - this$Estimated$pars["locN1",i]
      if(this$shape[i] == tvshape$double) st_c <- (this@st - this$Estimated$pars["locN1",i]) * (this@st - this$Estimated$pars["locN2",i])
      if(this$shape[i] == tvshape$double1loc) st_c <- (this@st - this$Estimated$pars["locN1",i])^2

      if(this$speedopt == speedopt$gamma) Gi <- 1/(1+exp(-this$Estimated$pars["speedN",i] * st_c))
      if(this$speedopt == speedopt$gamma_std) Gi <- 1/(1+exp(-this$Estimated$pars["speedN",i] * st_c/sd(this@st)))
      if(this$speedopt == speedopt$eta) Gi <- 1/(1+exp(-exp(this$Estimated$pars["speedN",i]) * st_c))

      g <- g + this$Estimated$pars["deltaN",i]*Gi
    }
  }
  #Return:
  g
}

get_g = function(Obj){

  objType <- class(Obj)
  if(objType[1] == "tv_class"){
    rtn <- calculate_g(Obj)
    return(rtn)
  }
  #
  if(objType[1] != "stcc1_class"){
    this <- Obj
    this$Estimated$pars <- this$pars
    this$Estimated$P1 <- this$P1
    this$Estimated$P2 <- this$P2
    rtn <- calc.Gt(this)
    return(as.vector(rtn))
  }
  # Else:
  warning("Only tv & stcc1 objects are supported")
  return(NaN)

}

.loglik.tv.univar = function(optimpars,e,tvObj,garchObj){
  this <- tvObj
  error <- -Inf

  # Copy the optimpars into a local tv_object
  if (isTRUE(this@delta0free)) {
    this$Estimated$delta0 <- optimpars[1]
    this$Estimated$pars <- .estimatedParsToMatrix(this,tail(optimpars,-1))
  } else{
    if(!("delta0" %in% names(this$Estimated))) this$Estimated$delta0 <- this$delta0
    this$Estimated$pars <- .estimatedParsToMatrix(this,optimpars)
  }

  # Do paramater boundary checks:
  # Check 1: Check that delta0 is positive
  if (this$Estimated$delta0 < 0) return(error)

  if (this@nr.transitions > 0){
    # We have some Tv$pars
    vecSpeed <- this$Estimated$pars["speedN",(1:this@nr.transitions)]
    vecLoc1 <- this$Estimated$pars["locN1",(1:this@nr.transitions)]
    vecLoc2 <- this$Estimated$pars["locN2",(1:this@nr.transitions)]

    # Check 2: Check the boundary values for speed params:
    #speedoptions: 1=gamma, 2=gamma/std(st), 3=exp(eta), 4=1/lambda^2
    maxSpeed <- switch(this$speedopt,1000,(1000/sd(this@st)),7.0,0.30)
    if (max(vecSpeed) > maxSpeed) return(error)
    if (min(vecSpeed) < 0) return(error)

    # Check 3: Check the loc1 locations fall within min-max values of st
    # We must have at least 1 loc1 to be inside this shape..loop, so no need to check if loc1 contains a valid value:
    if (min(vecLoc1) < min(this@st)) return(error)
    if (max(vecLoc1) > max(this@st)) return(error)

    # Check 4: Check that loc1.1 < loc1.2 .. locN.1 < locN.2 for all G(i)
    # Method: Subtract loc1_pos vector from loc2_pos vector and ensure it is positive:
    tmp <- vecLoc2 - vecLoc1
    # Note: tmp will contain NA wherever a loc2 element was NA - we can strip these out:
    if (sum(tmp < 0,na.rm = TRUE) > 0) return(error)

    # Check 5: Check the loc2 locations fall within min-max values of st
    # Confirm we have at least one valid numeric loc 2, before checking min & max:
    if (any(!is.na(vecLoc2))) {
      if (min(vecLoc2,na.rm = TRUE) < min(this@st)) return(error)
      if (max(vecLoc2,na.rm = TRUE) > max(this@st)) return(error)
    }

    # Check 6: Check that loc1.1 < loc2.1 where 2 locations exist... for all G(i)
    # We do need to have at least 2 locations for this error check
    if (NROW(vecLoc1) > 1) {
      v1 <- head(vecLoc1,-1)
      v2 <- tail(vecLoc1,-1)
      if (sum(v2-v1 < 0) > 0) return(error)
    }

  }# End: paramater boundary checks:

  g <- calculate_g(this)
  if (min(g,na.rm = TRUE) <= 0) return(error)

  h <- garchObj@h

  #Return the LogLiklihood value:
  ll <- loglik.tvgarch.univar(e,g,h)
  return(ll)

}
}


# ---  GARCH --- SECTION  ####

{## Help doco: garch() ####
  #' @noRd
#' @title
  #' A garch object
  #'
  #' @description
  #' `garch` is a
  #'
  #' @usage garch(type,order)
  #'
  #' @param type Use Enum, garchtype$...
  #' @param order Optional. Integer, defaults to 1, e.g. Default is a Garch(1,1) model.
  #'
  #' @details
  #' This object
  #'
  #' ```
  #'   myGarch = garch(garchtype$standard)
  #' ```
  #'
  #' This package provide support for standard Garch and GJR-Garch specifications
  #'
  #' @returns A garch_class object.
  #'
  #' @note
  #' I am a note
  #'
  #'
  #'
}

{## GARCH_CLASS Definition ----
garch <- setClass(Class = "garch_class",
                  slots = c(h="numeric",nr.pars="integer",omegafree="logical",order="numeric"),
                  contains = c("namedList")
)

setMethod("initialize","garch_class",
          function(.Object,...){
            .Object <- callNextMethod(.Object,...)
            .Object@h <- 1
            .Object@nr.pars <- as.integer(0)
            .Object@order <- 0
            .Object$type <- garchtype$noGarch
            .Object@omegafree <- TRUE
            .Object$pars <- matrix(NA,4,1)
            .Object$optimcontrol <- list(fnscale = -1, reltol = 1e-8) # Maximiser with default tolerance

            # Return:
            .Object
          })
}

{## Constructor: garch() ----
  setGeneric(name="garch",
             valueClass = "garch_class",
             signature = c("type","order"),
             def = function(type,order){

               garchtype <- list(noGarch=0,general=1,gjr=2)

               this <- new("garch_class")
               this$type <- type
               if(type == garchtype$noGarch) return(this)

               this@order <- order
               this <- .setInitPars(this)
               this$optimcontrol$ndeps <- rep(1e-3,this@nr.pars)
               if(type == garchtype$general) this$optimcontrol$parscale <- c(0.05,0.05,0.90)
               if(type == garchtype$gjr) this$optimcontrol$parscale <- c(0.005, 0.095, 0.8, 0.1) # Maximiser with default tolerance
               return(this)
             }
  )

  ## Constructor with 1 param - creates a Garch(1,1)
  setMethod("garch",signature = c("numeric","missing"),
            function(type){
              # Create a GARCH(1,1) model
              garch(type,c(1,1))
            })

}

{## plot(garch) override ----
setMethod("plot",signature = c(x="garch_class",y="missing"),
          function(x, y, ...){
            plot.default(x=sqrt(x@h), type='l', ylab = "sqrt(h)", ...)
          }
)
}

{## summary(garch) override ----
setMethod("summary",signature="garch_class",
          function(object,...){
            this <- object

            parsSummary <- NULL
            TypeNames <- c("No Garch","General","GJR Garch")

            if(!is.null(this$Estimated)){

              parsVec <-  round(as.vector(this$Estimated$pars),6)
              parsRows <- NROW(this$Estimated$pars)

              if(!is.null(this$Estimated$se) ){
                seVec <- round(as.vector(this$Estimated$se),6)
                seVecSig <- vector("character", length(seVec))

                for(n in seq_along(parsVec)){
                  if(is.na(seVec[n])) {
                    seVecSig[n] <- "   "
                  } else {
                    # Calculate a significance indicator
                    if(seVec[n]*2.576 < abs(parsVec[n]) ) { (seVecSig[n] <- "***") }
                    else if(seVec[n]*1.96 < abs(parsVec[n]) ) { (seVecSig[n] <- "** ") }
                    else if(seVec[n]*1.645 < abs(parsVec[n]) ) { (seVecSig[n] <- "*  ") }
                  }
                }
              } else {
                seVec <- rep(NaN,length(this$pars))
                seVecSig <- rep("   ", length(seVec))
              }

              seMat <- matrix(seVec,nrow=parsRows)
              colnames(seMat) <- paste("se" ,1:max(this@order),sep = "")
              # Build parsSummary table and insert the significance indicators
              parsSummary <- data.frame(NA,stringsAsFactors = FALSE)
              for (n in 1:NCOL(this$Estimated$pars)){
                sig <- matrix(seVecSig[1:parsRows],nrow=parsRows)
                parsSummary <- cbind(parsSummary,round(this$Estimated$pars[,n,drop=F],6),seMat[,n,drop=F],sig)
                seVecSig <- tail(seVecSig,-parsRows)
              }
            }

            cat("\nGARCH OBJECT\n")
            cat("\nType: ",TypeNames[this$type+1])
            cat("\nOrder: (",this@order[1],",",this@order[2],")")
            if(!is.null(this$Estimated)){
              cat("\nEstimation Results:\n")
              cat("\nMethod: ",this$Estimated$method,"\n")
              print(parsSummary[,-1])
              cat("\nLog-likelihood value(GARCH): ",this$Estimated$value)
            }

          }
  )
}


# ---  ESTIMATE GARCH --- SECTION  ####


{## Help doco: estimateGARCH() ####
#' @noRd
#' @title
#' Estimates a garch model
#'
#' @description
#' `estimateGARCH` is a
#'
#' @usage estimateGARCH(e,garchObj,estimationControl)
#'
#' @param e An estimated TV model
#' @param garchObj Use Enum, garchtype$...
#' @param estimationControl A list
#'
#' @details
#' This object
#'
#' ```
#'   myGarch = estimateGARCH(e,myGarch,estimationControl)
#' ```
#'
#'
#' @returns A garch_class object.
#'
#' @note
#' I am a note
#'
#'
}


{## -- estimateGARCH(e,TV,Garch,estCtrl) ----
#  estimateGARCH <- function(e,tvObj,garchObj,estimationControl){0}
  .estimateGARCH <- function(e,tvObj,garchObj,estimationControl){

    this <- garchObj

    # debug:
    #this <- GARCHspec
    #tvObj <- TVspec
    #estimationControl <- estCtrl

    # If there is no existing this$Estimated$, then create one
    if(!("Estimated" %in% names(this))) { this$Estimated <- list() }

    # Set verbose tracing:
    if (isTRUE(estimationControl$verbose)) {
      this$optimcontrol$trace <- 10
      cat("\nEstimating GARCH object...\n")
    } else this$optimcontrol$trace <- 0

    if(this$type == garchtype$noGarch) {
      message("Cannot estimateGARCH for type: NoGarch")
      return(this)
    }

    # Get the Optimpars from the garchObj (and update the tvObj as needed)
    garchOptimpars <- .getOptimpars_fromGARCH(this,estimationControl)
    optimpars <- garchOptimpars$Optimpars
    this <- garchOptimpars$garchObj

    # Check we have a g(t) vector, if not create one (so optim() doesn't complain):
    if(length(tvObj@g) == 1){ tvObj@g <- rep(1,NROW(e)) }

    ## --- Now call optim(.loglik.garch.univar) --- ##
    tmp <- NULL
    try(tmp <- optim(optimpars,.loglik.garch.univar,gr=NULL,e,this,tvObj, method="BFGS",control=this$optimcontrol))

    ## --- Attach results of estimation to the object --- ##

    if(!(is.null(tmp))){
      # Add the optim output to the estimated model
      this$Estimated$optimoutput <- tmp

      # Add the final results from optim() back into the estimated model
      this <- .setEstimatedPars_GARCH(this,tmp)

      # Get the conditional garch
      this@h <- calculate_h(e/sqrt(tvObj@g),this)

      # Calculate the parameter standard errors, if requested
      if (isTRUE(estimationControl$calcSE)) { this <- .setStdErrors_GARCH(this,tmp,tvObj) }
    }else{
      # optim() returned NULL.  What to do?
      # Return the object passed in with the Error = TRUE
      this$Estimated$error <- TRUE
      warning("estimateGARCH() generated a NULL return from optim()")
    }

    return(this)

  }
}

{##  loglik.garch.univar ----
  .loglik.garch.univar =  function(optimpars,e,garchObj,tvObj){


    error <- -1e6

    this <- garchObj

    ## ======== constraint checks ======== ##
    # Check if any parameter is negative:
    if(min(optimpars,na.rm = TRUE) <= 0) return(error)


    ## ======== calculate loglikelihood ======== ##


    if (isTRUE(this$omegafree)){
      # All pars are estimated
      if (optimpars[2] + optimpars[3] >= 0.9999) return(error)    #  constraint check: alpha & beta in [2,3]
      this$Estimated$pars <- .parsVecToMatrix(this,optimpars)
      names(this$Estimated$pars) <- rownames(this$pars)
    }else{
      if (optimpars[1] + optimpars[2] >= 0.9999) return(error)    #  constraint check: alpha & beta in [1,2]
      # omega must be calculated - and inserted back into the pars matrix.
      omega <- (1 - optimpars[1] - optimpars[2])
      this$Estimated$pars <- .parsVecToMatrix(this,c(omega,optimpars))
      names(this$Estimated$pars) <- rownames(this$pars)
    }



    ll <- tryCatch({
      # Get the g(t) vector
      g <- tvObj@g
      h <- calculate_h(e/sqrt(g),this)
      if (min(h,na.rm = TRUE) <= 0) return(error)
      #Return the LogLiklihood value:
      loglik.tvgarch.univar(e,g,h)

    }, warning = function(w){
      message("loglik.garch.univar warning: ", w$message)
    },error = function(err){
      message("loglik.garch.univar() error: ", err$message)
      ll <- error
    },finally = {
      # Any cleanup reqd.?
    })


    return(ll)

  }
}

{## estimateGARCH( - overloads -)  ----

  setGeneric("estimateGARCH",
             valueClass = "garch_class",
             function(e, tvObj, garchObj, estimationControl) {
               standardGeneric("estimateGARCH")
             }
  )

  # GARCH <- estimateGARCH(e,TV,GARCH,estCtrl)
  setMethod("estimateGARCH",
            signature = c(e="numeric",tvObj="tv_class",garchObj="garch_class",estimationControl="list"),
            function(e,tvObj,garchObj,estimationControl){
              .estimateGARCH(e,tvObj,garchObj,estimationControl)
            }
  )
  # GARCH <- estimateGARCH(e,TV,GARCH)
  setMethod("estimateGARCH",
            signature = c(e="numeric",tvObj="tv_class",garchObj="garch_class",estimationControl="missing"),
            function(e,tvObj,garchObj){
              estimationControl <- list(calcSE=FALSE, verbose=FALSE, maxIter=100, fixStartPars=FALSE, startParAdjust=10)
              .estimateGARCH(e,tvObj,garchObj,estimationControl)
            }
  )

  # GARCH <- estimateGARCH(e,GARCH,estCtrl)
  setMethod("estimateGARCH",
            signature = c(e="numeric",tvObj="garch_class",garchObj="list",estimationControl="missing"),
            function(e,tvObj,garchObj){
              st = (1:NROW(e))/NROW(e)
              tvMissing <- tv(1,tvshape$single)
              garchObj <- tvObj
              estimationControl <- garchObj
              .estimateGARCH(e,tvMissing,garchObj,estimationControl)
            }
  )
  # GARCH <- estimateGARCH(e,GARCH)
  setMethod("estimateGARCH",
            signature = c(e="numeric",tvObj="garch_class",garchObj="missing",estimationControl="missing"),
            function(e,tvObj){
              st = (1:NROW(e))/NROW(e)
              tvMissing <- tv(1,tvshape$single)
              garchObj <- tvObj
              estimationControl <- list(calcSE=FALSE, verbose=FALSE, maxIter=100, fixStartPars=FALSE, startParAdjust=10)
              .estimateGARCH(e,tvMissing,garchObj,estimationControl)
            }
  )
}

{## estimateGARCH() - helpers ----
.getOptimpars_fromGARCH <- function(garchObj,estimationControl){

  # Derive the optimpars vector, allowing for NA's and omega0free:
  # Return both the optimpars vector and the garchObj, with corrected internals, e.g. nr.pars, optimControl list.

  this <- garchObj
  estCtrl <- estimationControl

  #debug:
  #this <- GARCH

  # Check for iterative-estimator controls, Set defaults if not exists:
  if(!("fixStartPars" %in% names(estCtrl))) { estCtrl$fixStartPars <- FALSE }
  if(!("startParAdjust" %in% names(estCtrl))) { estCtrl$startParAdjust <- 50 }

  # Step1. Get the starting pars (parsVec):
  if(isTRUE(estCtrl$fixStartPars)  || (!("pars" %in% names(this$Estimated))) ){
    # Do nothing
    parsVec <- as.vector(this$pars)
  }else{
    # Get last Estimates:
    parsVec <- as.vector(this$Estimated$pars)  #TODO: Confirm this is correct for TVGARCH estimator
  }

  # Step2. Remove any NA's (GJR.gamma params) - Note: This is CRITICAL as length(parsVec is used below)
  parsVec <- parsVec[!is.na(parsVec)]

  # Step3.  When omega is NOT a free parameter, remove it from parsVec and the optimcontrol's
  if(isTRUE(this@omegafree)){
    # parsVec is correct
    names(parsVec) <- rownames(this$pars)
    this@nr.pars <- as.integer(length((parsVec)))
    # TODO: BUG here if user switches omegafree ON, then OFF, then ON again the optimcontrol gets whacked!
    if(length(this$optimcontrol$ndeps) < length(parsVec)) { this$optimcontrol$ndeps <- c(1e-3,this$optimcontrol$ndeps) }
    if(length(this$optimcontrol$parscale) < length(parsVec)) { this$optimcontrol$parscale <- c(0.05,this$optimcontrol$parscale) }
  }else{
    # VarTargetting, calculate omega, estimate alpha & beta:
    parsVec <- tail(as.vector(this$pars), -1)
    names(parsVec) <- tail(rownames(this$pars),-1)
    this@nr.pars <- as.integer(length((parsVec)))
    # Drop omega from the OptimControl's
    if(length(this$optimcontrol$ndeps) != length(parsVec)) this$optimcontrol$ndeps <- tail(this$optimcontrol$ndeps,this@nr.pars)
    if(length(this$optimcontrol$parscale) != length(parsVec)) this$optimcontrol$parscale <- tail(this$optimcontrol$parscale,this@nr.pars)
  }

  # Shake the starting pars (parsVec):
  boundSize <- estCtrl$startParAdjust *  this$optimcontrol$ndeps    # Should be in range 10 - 100
  boundLo <- parsVec - boundSize
  boundHi <- parsVec + boundSize

  # # Anna's method
  # boundSize <- estCtrl$startParAdjust  # Should be small, in the order of 1e-2
  # boundLo <- parsVec/this$optimcontrol$parscale - boundSize
  # boundHi <- parsVec/this$optimcontrol$parscale + boundSize
  # boundLo <- boundLo * this$optimcontrol$parscale
  # boundHi <- boundHi * this$optimcontrol$parscale

  parsVec <- runif(this@nr.pars,boundLo,boundHi)

  #RETURN:
  rtnList <- list()
  rtnList$Optimpars <- parsVec
  rtnList$garchObj <- this
  return(rtnList)

}

.setEstimatedPars_GARCH <- function(garchObj,optimTmp){

  this <- garchObj
  tmp <- optimTmp

  # An unhandled error could result in a NULL being returned by optim()
  if (is.null(tmp)) {
    this$Estimated$value <- NA
    this$Estimated$error <- TRUE
    warning("estimateGARCH() - optim failed unexpectedly and returned NULL. Check the optim controls & starting params")
    return(this)
  }
  if (tmp$convergence!=0) {
    this$Estimated$value <- NA
    this$Estimated$error <- TRUE
    this$Estimated$optimoutput <- tmp
    warning("estimateGARCH() - failed to converge. Check the optim controls & starting params")
    return(this)
  }

  # No optim issues, so set output
  this$Estimated$value <- tmp$value
  this$Estimated$error <- FALSE

  #Update the GARCH object parameters using optimised pars:
  if (isTRUE(this@omegafree)){
    this$Estimated$pars <- .parsVecToMatrix(this,tmp$par)
    this@nr.pars <- as.integer(3)

  } else{
    omega <- (1 - tmp$par[1] - tmp$par[2])
    this$Estimated$pars <- .parsVecToMatrix(this,c(omega,tmp$par))
    this@nr.pars <- as.integer(2)
  }
  colnames(this$Estimated$pars) <- paste("Est")

  #RETURN:
  return(this)

}

.setStdErrors_GARCH <- function(garchObj,optimTmp,tvObj){

  this <- garchObj
  tmp <- optimTmp

  cat("\nCalculating GARCH standard errors...\n")
  this$Estimated$hessian <- NULL
  this$Estimated$se <- NULL
  StdErrors <- NULL

  # this$Estimated$se <- matrix(NA,nrow=3,ncol=1)
  # rownames(this$Estimated$se) <- rownames(this$pars)
  # colnames(this$Estimated$se) <- "se"

  # Get the hessian from the optimiser:
  try(this$Estimated$hessian <- optimHess(tmp$par,.loglik.garch.univar,gr=NULL,e,this,tvObj,control=this$optimcontrol))

  # Attempt to invert it
  try(StdErrors <- sqrt(-diag(invertHess(this$Estimated$hessian))))

  # Handle invertHess errors
  if(!is.null(StdErrors)) {
    parsVec <-  as.vector(this$pars)

    if(isTRUE(this@omegafree)){
      this$Estimated$se <- matrix(StdErrors,nrow=3,ncol=1)
    }else{
      # omega is calculated, so has no stderror
      this$Estimated$se <- matrix(c(NA,StdErrors),nrow=3,ncol=1)

      #this$Estimated$se[1,1] <- NA
      #this$Estimated$se[2,1] <- StdErrors[1]
      #this$Estimated$se[3,1] <- StdErrors[2]
    }
    colnames(this$Estimated$se) <- "se"
  }
  #RETURN:
  return(this)
}

.setInitPars = function(garchObj){
  this <- garchObj

  if(this$type == garchtype$noGarch) {
    message("Cannot create Garch$pars for type: NoGarch")
    return(this)
  }

  maxLag <- max(this@order)

  # Set the row names:
  GarchparsRownames <- c("omega","alpha","beta","gamma")

  if(this$type == garchtype$general) {
    this@nr.pars <- as.integer(3)
    pars <- matrix(nrow = this@nr.pars,ncol = maxLag)
    rownames(pars) <- GarchparsRownames[1:this@nr.pars]
    for(n in 1:maxLag){
      pars["omega",n] <- 0.05
      pars["alpha",n] <- 0.05
      pars["beta",n] <- 0.90
    }
    this$pars <- pars
  }
  if(this$type == garchtype$gjr) {
    this@nr.pars <- as.integer(4)
    pars <- matrix(nrow = this@nr.pars,ncol = maxLag)
    rownames(pars) <- GarchparsRownames[1:this@nr.pars]
    for(n in 1:maxLag){
      pars["omega",n] <- 0.005
      pars["alpha",n] <- 0.095
      pars["beta",n] <- 0.80
      pars["gamma",n] <- 0.10
    }
    this$pars <- pars
  }
  ## TODO: Implement in future release
  # if(this$type == garchtype$gjr_alpha0) {
  #   this@nr.pars <- as.integer(3)
  #   pars <- matrix(nrow = this@nr.pars,ncol = maxLag)
  #   rownames(pars) <- GarchparsRownames[1:(this@nr.pars + 1)]
  #   for(n in 1:maxLag){
  #     pars["omega",n] <- 0.05
  #     pars["alpha",n] <- 0.0
  #     pars["beta",n] <- 0.85
  #     pars["gamma",n] <- 0.05
  #   }
  #   this$pars <- pars
  # }

  return(this)
}

.parsVecToMatrix <-function(garchObj,pars){
  this <- garchObj

  # convert the passed-in 'pars' vector to the Garch$$=pars matrix:

  if(this$type == garchtype$noGarch) {
    message("Cannot create Garch Params for type: NoGarch")
    return(this)
  }
  #maxLag <- max(this@order)
  maxLag <- 1

  # Set the row names:
  # TODO: remove GJR capability from package
  #garchparsRownames <- c("omega","alpha","beta","gamma")
  garchparsRownames <- c("omega","alpha","beta")

  # pars contains all Garch paramaters - the estimateGARCH() does this
  # Return the formatted matrix
  #matrix(pars,nrow = this@nr.pars ,ncol = maxLag,dimnames = list(garchparsRownames[1:this@nr.pars],"Est"))
  matrix(pars, nrow=3 ,ncol=1, dimnames = list(garchparsRownames,"Est"))

}

calculate_h <- function(e,garchObj){
  this <- garchObj

  if(this$type == garchtype$noGarch) return(this@h)

  # First Run = No this$Estimated$pars exists:
  if(!("pars" %in% names(this$Estimated))) { this$Estimated$pars <- this$pars }

  Tobs <- NROW(e)
  h <- rep(0,Tobs)
  h[1] <- sum(e*e)/Tobs

  # TODO: Remove GJR Garch
  for(t in 2:Tobs) {
    if(isTRUE(this@omegafree)){
      h[t] <- this$Estimated$pars["omega",1] + this$Estimated$pars["alpha",1]*(e[t-1])^2 + this$Estimated$pars["beta",1]*h[t-1]
    }else{
      h[t] <- (1 - this$Estimated$pars["alpha",1] - this$Estimated$pars["beta",1]) + this$Estimated$pars["alpha",1]*(e[t-1])^2 + this$Estimated$pars["beta",1]*h[t-1]
    }
    #if(this$type == garchtype$gjr) h[t] <- h[t] + this$Estimated$pars["gamma",1]*(min(e[t-1],0))^2
  }

  return(h)
}

get_h = function(garchObj,e){
  this <- garchObj

  if(this$type == garchtype$noGarch) return(this@h)

  Tobs <- NROW(e)
  h <- rep(0,Tobs)
  h[1] <- sum(e*e)/Tobs

  # TODO: Extend the below to handle more lags (higher order Garch)
  for(t in 2:Tobs) {
    h[t] <- this$pars["omega",1] + this$pars["alpha",1]*(e[t-1])^2 + this$pars["beta",1]*h[t-1]
    if(this$type == garchtype$gjr) h[t] <- h[t] + this$pars["gamma",1]*(min(e[t-1],0))^2
  }

  return(h)
}

}


{##  Help doco: estimateGARCH_RollingWindow  ####
  #' @noRd
#' @title
  #' Approximately estimates a garch model
  #'
  #' @description
  #' `estimateGARCH_RollingWindow` is a
  #'
  #' @usage estimateGARCH_RollingWindow(e,garchObj,estimationControl)
  #'
  #' @param e An estimated TV model
  #' @param garchObj Use Enum, garchtype$...
  #' @param estimationControl A list
  #'
  #' @details
  #' This object
  #'
  #' ```
  #'   estimationControl$vartargetWindow = 400
  #'   myGarch = estimateGARCH_RollingWindow(e,myGarch,estimationControl)
  #' ```
  #'
  #'
  #' @returns A garch_class object.
  #'
  #' @note
  #' I am a note
  #'
  #'
}

{## estimateGARCH_RollingWindow(e,garchObj,estCtrl)  ----
  estimateGARCH_RollingWindow <- function(e,garchObj,estimationControl){0}
  .estimateGARCH_RollingWindow <- function(e,garchObj,estimationControl){
    this <- garchObj

    # == Validations == #
    if(this$type == garchtype$noGarch) {
      message("Cannot estimateGARCH for type: NoGarch")
      return(this)
    }

    if(!is.list(estimationControl)){
      warning("A valid estimationControl list is required - see Help")
    }

    if(is.null(estimationControl$vartargetWindow)) {
      warning("A valid estimationControl$vartargetWindow length is required.\nA default value of 500 observations will be used - see Help for details.")
      vartargetWindow <- 500
    } else if(estimationControl$vartargetWindow <= 0) {
      vartargetWindow <- 500
    } else vartargetWindow <- estimationControl$vartargetWindow
    # == End: Validations == #

    # Attach results of estimation to the object
    this$Estimated <- list()

    #
    if(!is.null(estimationControl$calcSE)) calcSE <- estimationControl$calcSE else calcSE <- TRUE
    if(!is.null(estimationControl$verbose)) verbose <- estimationControl$verbose else verbose <- TRUE

    if (isTRUE(verbose)) {
      this$optimcontrol$trace <- 10
      cat("\nEstimating GARCH object...\n")
    } else this$optimcontrol$trace <- 0

    # vartargetWindow is only used for initial estimation of the Garch
    this$Estimated$method <- paste0("MLE, variance-targetting a rolling Window of ",vartargetWindow, " observations")

    # Get Optimpars from garch$pars & remove omega
    optimpars <- as.vector(this$pars)
    names(optimpars) <- rownames(this$pars)
    # When var-targetting, omega is calculated, so...
    optimpars <- tail(optimpars,-1)
    this$optimcontrol$ndeps <- tail(this$optimcontrol$ndeps,-1)
    this$optimcontrol$parscale <- tail(this$optimcontrol$parscale,-1)

    # Now call optim:
    tmp <- NULL
    try(tmp <- optim(optimpars,.loglik.garch.rollingWin,gr=NULL,e,this,vartargetWindow, method="BFGS",control=this$optimcontrol))

    # An unhandled error could result in a NULL being returned by optim()
    if (is.null(tmp)) {
      this$Estimated$value <- -Inf
      this$Estimated$error <- TRUE
      warning("estimateGARCH() - optim failed and returned NULL. Check the optim controls & starting params")
      return(this)
    }
    if (tmp$convergence!=0) {
      this$Estimated$value <- -Inf
      this$Estimated$error <- TRUE
      this$Estimated$optimoutput <- tmp
      warning("estimateGARCH() - failed to converge. Check the optim controls & starting params")
      return(this)
    }

    this$Estimated$value <- tmp$value
    this$Estimated$error <- FALSE

    #Update the GARCH object paramters using optimised pars:
    omega <- 1 - tmp$par[1] - tmp$par[2]
    this$Estimated$pars <- .parsVecToMatrix(this,c(omega,tmp$par))

    # Get conditional variance - 'h'
    Tobs <- NROW(e)
    h <- rep(0,Tobs)
    h[1] <- sum(e*e)/Tobs
    halfWindow <- round(vartargetWindow/2)

    for(t in 2:Tobs){

      if( t <= halfWindow || t > (Tobs-halfWindow) ) {
        h[t] <- var(e)
      } else {
        start <- t - halfWindow
        end <- min(t + halfWindow,Tobs)
        localVar <- var(e[start:end])
        this$Estimated$pars["omega",1] <- (1 - this$Estimated$pars["alpha",1] - this$Estimated$pars["beta",1]) * localVar
        h[t] <- this$Estimated$pars["omega",1] + this$Estimated$pars["alpha",1]*(e[t-1])^2 + this$Estimated$pars["beta",1]*h[t-1]
        if(this$type == garchtype$gjr) h[t] <- h[t] + this$Estimated$pars["gamma",1]*(min(e[t-1],0))^2
      }
    }
    # End: Get conditional variance - 'h'


    # Calc Std Errors
    if (isTRUE(calcSE)) {
      this$Estimated$hessian <- NULL
      try(this$Estimated$hessian <- optimHess(tmp$par,.loglik.garch.rollingWin,gr=NULL,e,this,vartargetWindow, control=this$optimcontrol))
      # Handle optimHess errors
      StdErrors <- NULL
      try(StdErrors <- sqrt(-diag(invertHess(this$Estimated$hessian))))
      if(is.null(StdErrors)) {
        this$Estimated$se <- matrix(NA,nrow=(this@nr.pars))
      }else {
        StdErrors <- c(NA,StdErrors)
        this$Estimated$se <- matrix(StdErrors,nrow=(this@nr.pars))
      }
      rownames(this$Estimated$se) <- rownames(this$pars)
      colnames(this$Estimated$se) <- "se"
    }
    if (isTRUE(verbose)) this$Estimated$optimoutput <- tmp

    return(this)
  }
}

{## estimateGARCH_RollingWindow( - overloads -)  ----
  setGeneric("estimateGARCH_RollingWindow",valueClass = "garch_class")

  setMethod("estimateGARCH_RollingWindow",
            signature = c(e="numeric",garchObj="garch_class",estimationControl="list"),
            function(e,garchObj,estimationControl){
              .estimateGARCH_RollingWindow(e,garchObj,estimationControl)
            }
  )
  setMethod("estimateGARCH_RollingWindow",
            signature = c(e="numeric",garchObj="garch_class",estimationControl="missing"),
            function(e,garchObj){
              estimationControl <- list(calcSE=TRUE,verbose=TRUE,vartargetWindow=500)
              .estimateGARCH_RollingWindow(e,garchObj,estimationControl)
            }
  )
}







.loglik.garch.rollingWin =  function(optimpars,e,garchObj,vartargetWindow){

  error <- -Inf
  this <- garchObj

  ## ======== constraint checks ======== ##
  # Check if any parameter is negative:
  if(min(optimpars,na.rm = TRUE) <= 0) return(error)
  if (optimpars["alpha"] + optimpars["beta"] >= 1) return(error)

  ## ======== calculate loglikelihood ======== ##

  estPars <- optimpars
  omega <- 1 - optimpars[1] - optimpars[2]
  estPars <- c(omega,optimpars)
  this$Estimated$pars <- .parsVecToMatrix(this,estPars)

  # Get conditional variance - 'h'
  Tobs <- NROW(e)
  h <- rep(0,Tobs)
  h[1] <- sum(e*e)/Tobs
  halfWindow <- round(vartargetWindow/2)

  for(t in 2:Tobs){

    if( t <= halfWindow || t > (Tobs-halfWindow) ) {
      h[t] <- var(e)
    } else {
      start <- t - halfWindow
      end <- min(t + halfWindow,Tobs)
      localVar <- var(e[start:end])
      this$Estimated$pars["omega",1] <- (1 - this$Estimated$pars["alpha",1] - this$Estimated$pars["beta",1]) * localVar
      h[t] <- this$Estimated$pars["omega",1] + this$Estimated$pars["alpha",1]*(e[t-1])^2 + this$Estimated$pars["beta",1]*h[t-1]
      if(this$type == garchtype$gjr) h[t] <- h[t] + this$Estimated$pars["gamma",1]*(min(e[t-1],0))^2
    }
  }
  # End: Get conditional variance - 'h'

  #Return the LogLiklihood value:
  g <- rep(1,NROW(e))
  ll <- loglik.tvgarch.univar(e,g,h)
  return(ll)

}





# ---  TVGARCH --- SECTION  ####

{## Help doco: tvgarch()  ####
  #' @noRd
#' @title
  #' A multiplicative tvgarch object
  #'
  #' @description
  #' `tvgarch` is a
  #'
  #' @usage tvgarch(tvObj,garchType)
  #'
  #' @param tvObj An estimated TV model
  #' @param garchType Use Enum, garchtype$...
  #'
  #' @details
  #' This object
  #'
  #' ```
  #'   myTvGarch = tvgarch(myTV,garchtype$standard)
  #' ```
  #'
  #'
  #' @returns A tvgarch_class object.
  #'
  #' @note
  #' I am a note
  #'
  #'
}

{## TVGARCH_CLASS Definition ----
  tvgarch <- setClass(Class = "tvgarch_class",
                      slots = c(Tobs="integer",tvObj="tv_class",garchObj="garch_class",e="numeric",iterations="integer"),
                      contains = c("namedList")
  )

  setMethod("initialize","tvgarch_class",
            function(.Object){
              .Object@Tobs <- as.integer(0)
              .Object@e <- vector("numeric")
              # Number of iterations required for estimation to converge
              .Object@iterations <- as.integer(0)    # internal counter for iterative estimator

              #
              .Object$varTarget <- TRUE

              # TODO:  Should we drop all these copies & use the tvObj,garchObj instead?  Or keep these as the Estimated values??
              # Default TV properties
              .Object$tvshape <- tvshape$delta0only
              .Object$tvspeedopt <- speedopt$none
              .Object$tvdelta0 <- 1.0
              .Object$tvpars <- matrix(NA,4,1)
              .Object$tvOptimcontrol <- list(fnscale = -1, reltol = 1e-8)

              # Default GARCH properties
              .Object$garchtype <- garchtype$noGarch
              .Object$garchpars <- matrix(NA,4,1)
              .Object$garchOptimcontrol <- list(fnscale = -1, reltol = 1e-8)
              # user sets iteration-convergence threshold (loglik value)
              .Object$iterationReltol <- 1e-8
              # Name of Data Series - for plotting & historical reference
              .Object$data_desc <- "mtvgarch plot"

              # Return:
              .Object
            })
}

{## Constructor: tvgarch()  ----
  setGeneric(name="tvgarch",
             valueClass = "tvgarch_class",
             signature = c("tvObj","garchObj"),
             def = function(tvObj,garchObj){

               # Validate: Spit dummy if TV or GARCH are estimated
               if(("Estimated" %in% names(tvObj)) || ("Estimated" %in% names(tvObj)) ) {
                 message("tvgarch-class objects require the tv & garch parameters to be specified, but NOT estimated.")
                 return(this)
               }

               this <- new("tvgarch_class")
               this@tvObj <- tvObj
               this@garchObj <- garchObj

               this@Tobs <- tvObj@Tobs
               # TV
               this$tvshape <- tvObj$shape
               this$tvspeedopt <- tvObj$speedopt
               this$tvdelta0 <- tvObj$delta0
               this$tvpars <- tvObj$pars
               this$tvOptimcontrol <- this@tvObj$optimcontrol

               # GARCH
               this$garchtype <- this@garchObj$type
               this$garchpars <- this@garchObj$pars
               this$garchOptimcontrol <- this@garchObj$optimcontrol

               cat("\ntvgarch object created successfully!\n")

               return(this)
             }
  )
}

{## plot(tvgarch override ----
  setMethod("plot",signature = c(x="tvgarch_class",y="missing"),
            function(x, y,...){
              #TODO: Allow override of main/title by user
              if (is.na(x$e_desc)) title = "" else title = x$e_desc
              g <- x$Estimated$g
              h <- x$Estimated$h
              plot.default(x=sqrt(g), type='l', ylab = "sqrt(g)", main=title, ...)
              plot.default(x=sqrt(h), type='l', ylab = "sqrt(h)", main=title, ...)
            })
}

{## summary(tvgarch) override ----
  setMethod("summary",signature="tvgarch_class",
            function(object,...){
              this <- object
              if(is.null(this$Estimated)){
                #cat("\n\nPlease estimate the TVGARCH Model first")
                return("Please estimate the TVGARCH Model first")
              }

              cat("\n -- TVGARCH Model Specification --\n")
              cat("\nMultiplicative Model Log-Likelihood Value: ", this$Estimated$value)
              cat("\n\nTVGARCH Model Parameters:")
              ## TODO:  Fix the Summary - it stopped working  :(   Changing to 'call' didn't help
              call("summary",this@garchObj)
              call("summary",this@tvObj)
              cat("\n\n -- End of TVGARCH Model Specification --")

            }
  )
}

# ---  ESTIMATE TVGARCH --- SECTION  ####

{## Help doco:  estimateTVGARCH_Iterate ####
#' @noRd
#' @title
#' Estimates a tv model
#'
#' @description
#' `estimateTVGARCH_Iterate` is a
#'
#' @usage estimateTVGARCH_Iterate(e,tvgarchObj,estimationControl)
#'
#' @param e An
#' @param tvgarchObj Use
#' @param estimationControl A list
#'
#' @details
#' This object
#'
#' ```
#'   myTvGarch = estimateTVGARCH_Iterate(e,myTvGarch,estimationControl)
#' ```
#'
#'
#' @returns A tvgarch_class object.
#'
#' @note
#' I am a note
#'
#'
}


{# estimateTVGARCH ----
  .estimateTVGARCH <- function(e,tvgarchObj,estimationControl){

    this <- tvgarchObj
    # debug:
    #this <- mod
    #estimationControl <- estCtrl

    # Get the TV & GARCH specifications/models from the tvgarchObj param
    TV <- this@tvObj
    GARCH <- this@garchObj

    # Validate: Spit dummy if TV or GARCH are estimated
    if(("Estimated" %in% names(TV)) || ("Estimated" %in% names(GARCH)) ) {
      message(".estimateTVGARCH() requires the tv & garch parameters to be specified, but NOT estimated.")
      return(this)
    }else{
      # create the Estimated list now
      this$Estimated <- list()
    }

    # Cache the data
    this@e <- e

    # Get maxIterations from the control list
    # Note if the user didn't set estimation controls the following defaults are provided by an overload method:
    #estimationControl <- (calcSE=FALSE, verbose=FALSE, maxIter=100, fixStartPars=FALSE, startParAdjust=10)

    maxIterations <- estimationControl$maxIter
    if(!(is.numeric(this$iterationReltol) )) {this$iterationReltol <- 1e-5}

    # NOte: the fixStartPars and startParAdjust parameters will be passed down to the underlying Estimators

    # Configure variance targetting
    if(isTRUE(this$varTarget)){
      TV@delta0free <- TRUE
      GARCH@omegafree <- FALSE
    }else{
      TV@delta0free <- FALSE
      GARCH@omegafree <- TRUE
    }

    # Start Iterative TVGARCH Estimation ####
    cat("\nStarting Iterative TVGARCH Estimation...\n")

    # While..Loop => Ensures we always run at least once.
    # Update this@iterations after Estimates (checking for Convergence or MaxIterations)
    # If an estimation FAILS inside the loop: Set the Result = 'Failed' and Exit immediately using return(this)

    # Set starting values to ensure we do at least 1 loop:
    last.tvg.loglik <- -1e6
    tvg.loglik <- 1
    this@iterations <- as.integer(0)



    while(isTRUE(this@iterations <= maxIterations) && abs(tvg.loglik - last.tvg.loglik) > this$iterationReltol ){

      # 1. Copy current to last:
      last.tvg.loglik <- tvg.loglik

      # 2. Increment the iteration counter
      this@iterations <- this@iterations + as.integer(1)

      # 3. Run estimations and calculate a new tvg.loglik
      TV <- tryCatch({
        estimateTV(e,TV,GARCH,estimationControl)
      }, warning = function(w){
        message("Caught a warning: ", w$message)
        # TODO: Find a better way to bubble up the errors from inside the estimators!
        if(grepl("NULL", w$message)){
          # This is in-fact an error for this process
          this$Estimated$error <- TRUE
          this$Estimated$converged <- FALSE
          if(isTRUE(estimationControl$verbose)){ cat("\nLast TV Estimate generated a NULL optim() result - Stopping iterator now\n") }
          break
        }
      },error = function(err){
        message("estimateTVGARCH() error: ", err$message)
        this$Estimated$error <- TRUE
        this$Estimated$converged <- FALSE
        if(isTRUE(estimationControl$verbose)){ cat("\nLast TV Estimate generated an Error - Stopping iterator now\n") }
        break
      },finally = {
        # Any cleanup reqd.?
      })

      if(isTRUE(estimationControl$verbose)){cat(".")}  # Single line comment, TV done

      # The following GARCH estimation is done even if the prior TV did not improve the Loglik value
      GARCH <- tryCatch({
        estimateGARCH(e,TV,GARCH,estimationControl)
      }, warning = function(w){
        message("estimateTVGARCH() warning: ", w$message)
        if(grepl("NULL", w$message)){
          # This is in-fact an error for this process
          this$Estimated$error <- TRUE
          this$Estimated$converged <- FALSE
          if(isTRUE(estimationControl$verbose)){ cat("\nLast GARCH Estimate generated a NULL optim() result - Stopping iterator now\n") }
          break
        }
      },error = function(err){
        message("estimateTVGARCH() error: ", err$message)
        this$Estimated$error <- TRUE
        this$Estimated$converged <- FALSE
        if(isTRUE(estimationControl$verbose)){ cat("\nLast GARCH Estimate generated an Error\nStopping iterator now\n") }
        break
      },finally = {

        # 4. Calc the loglik value ONLY when TV AND GARCH have successfully estimated
        if(class(TV)=="tv_class" && class(GARCH)=="garch_class"){
          # Confirm the classes are valid
          tvg.loglik <- loglik.tvgarch.univar(e,g=TV@g,h=GARCH@h)
        }else{
          # What should we do if one of the estimator failed - but we got this far??
        }

      })

      # Calculate the initial/starting Loglik.Value (Use to determine convergence/divergence)
      if(this@iterations == as.integer(1)){ startLoglik <- loglik.tvgarch.univar(e,g=TV@g,h=GARCH@h) }
      #
      if(isTRUE(estimationControl$verbose)){cat(".")}  # Single line comment, Garch done

    }  # End: While() Loop

    # TODO: Tidy up the wrap-up logic below
    if(isTRUE(this$Estimated$error)) return(this)


    # If we have exceded the maxIterations, usually that means we diverged.
    # The exception to this rule is: User wants to do a 2-Step Estimation (allow user to specify 1-Step or 2-Step)
    if((this@iterations >= maxIterations) && (maxIterations > 2)){
      # Model diverged - return estimation failure
      this$Estimated$value <- tvg.loglik  # May be slow convergence, the value may still be good
      this$Estimated$error <- FALSE
      if(tvg.loglik > startLoglik) {
        # TODO: Check this logic with Anna
        # Slow convergence
        this$Estimated$converged <- TRUE
      } else {
        # Diverging
        this$Estimated$converged <- FALSE
      }
    }else{
      # Model Converged - return estimation Success
      this$Estimated$value <- tvg.loglik
      this$Estimated$error <- FALSE
      this$Estimated$converged <- TRUE
    }

    # Always Update the internal objects with the Estimated objects:
    this <- .updateEstimatedDetails(this,TV,GARCH)

    return(this)
  }
}


{# estimateTVGARCH() overrides ----

  setGeneric("estimateTVGARCH",
             valueClass = "tvgarch_class",
             function(e,tvgarchObj,estimationControl) {
               standardGeneric("estimateTVGARCH")
             }
  )
  # TVG <- estimateTVGARCH(e,TVG,estCtrl)
  setMethod("estimateTVGARCH",
            signature = c(e="numeric",tvgarchObj="tvgarch_class",estimationControl="list"),
            function(e,tvgarchObj,estimationControl){
              .estimateTVGARCH(e,tvgarchObj,estimationControl)
            }
  )
  # TVG <- estimateTVGARCH(e,TVG,estCtrl)
  setMethod("estimateTVGARCH",
            signature = c(e="numeric",tvgarchObj="tvgarch_class",estimationControl="missing"),
            function(e,tvgarchObj){
              estimationControl <- list(calcSE=FALSE, verbose=FALSE, maxIter=100, fixStartPars=FALSE, startParAdjust=10)
              .estimateTVGARCH(e,tvgarchObj,estimationControl)
            }
  )
}


loglik.tvgarch.univar = function(e,g,h){

  # Calculate the tvgarch logliklihood value.  Depends on g(t) & h(t)

  ll <- tryCatch({
    sum( -0.5*log(2*pi) - 0.5*log(g) - 0.5*log(h) - 0.5*(e^2/(g*h) ) )
  }, warning = function(w){
    message("loglik.tvgarch.univar() warning: ", w$message)
  },error = function(err){
    message("loglik.tvgarch.univar() error: ", err$message)
    ll <- -1e6  #TODO: Is this a good idea???  (Trying to allow the Iterative estimator to continue)
  },finally = {
    #
  })


  #ll <- sum( -0.5*log(2*pi) - 0.5*log(g) - 0.5*log(h) - 0.5*(e^2/(g*h) ) )
  #names(ll) <- "Loglik.Value"
  return(ll)
}

.updateEstimatedDetails <- function(tvgarchObj,TV,GARCH){

  this <- tvgarchObj

  if(!(is.null(this))){
    # Put this valid model into the Estimated list:
    this$Estimated$tv <- TV$Estimated
    this$Estimated$tv$g <- TV@g

    this$Estimated$garch <- GARCH$Estimated
    this$Estimated$garch$h <- GARCH@h

    # Populate the convenience attributes:
    this$Estimated$g <- TV@g
    this$Estimated$h <- GARCH@h
  }else{
    warning("estimateTVGARCH() completed with a NULL object")
  }

  return(this)

}



# ========================================== #
# TODO:  Add unit tests!!!


# ---  TESTS --- SECTION  ####

{# Help doco: test.LM.TR2(e,tv,ord) ####
  #' @noRd
#' @title
  #' LM TR-Squared Test
  #'
  #' @description
  #' `test.LM.TR2` is a
  #'
  #' @usage test.LM.TR2(e,tvObj,testOrder)
  #'
  #' @param e An estimated TV model
  #' @param tvObj Use Enum, garchtype$...
  #' @param testOrder An integer
  #'
  #' @details
  #' This object
  #'
  #' ```
  #'  testOrder = 1
  #'  testResult = test.LM.TR2(e,myTv,testOrder)
  #' ```
  #'
  #'
  #' @returns A numeric test statistic
  #'
  #' @note
  #' I am a note
  #'
  #'
}

test.LM.TR2=function(e,tvObj,testOrder){
  this <- tvObj

  if(testOrder <= 0){
    message("Cannot execute test with no alternate hypothesis. Please set a valid testOrder")
    return(NaN)
  }

  # Test Method: Regress psi2_1 on 1/gt*(dgdt and dgdt2)
  # 1. Calc derivatives of params dgdt = Tx1 or Tx4 or Tx7 or...
  #    NCOL(dgdt) increases with the order of TV function.
  dgdt <- .dg_dt(this)

  # 2. Calc derivatives of taylor pars (linearised component) under the null
  dgdt2 <- .dg_dt2(this@st,testOrder)

  g <- calculate_g(this)
  X <- cbind(dgdt,dgdt2)/g

  # 3. Invert crossprod(X) to calculate SSR1
  Xinv <- NULL
  try(Xinv <- solve(crossprod(X)))
  if (is.null(Xinv)){
    rm(g,dgdt,dgdt2)
    return(NaN)
  }

  # 4. Calculate psi2_1 to calculate SSR0
  psi2_1 <- matrix(data=(e^2/g-1),nrow = this@Tobs,ncol = 1)

  # 5. Calc the TestStat:
  SSR0 <- sum(psi2_1*psi2_1)    # Scalar
  SSR1 <- sum((psi2_1-X%*%Xinv%*%t(X)%*%psi2_1)^2)

  Result <- this@Tobs*(SSR0-SSR1)/SSR0

  # Tidy up & release memory before returning:
  rm(this,psi2_1,g,X,Xinv,dgdt,dgdt2)

  # Return:
  Result

}


{# Help doco: test.LM.Robust(e,tv,ord) ####
  #' @noRd
#' @title
  #' LM TR-Squared Test - Robust Version
  #'
  #' @description
  #' `test.LM.Robust` is a
  #'
  #' @usage test.LM.Robust(e,tvObj,testOrder)
  #'
  #' @param e An estimated TV model
  #' @param tvObj Use Enum, garchtype$...
  #' @param testOrder An integer
  #'
  #' @details
  #' This object
  #'
  #' ```
  #'  testOrder = 1
  #'  testResult = test.LM.Robust(e,myTv,testOrder)
  #' ```
  #'
  #'
  #' @returns A numeric test statistic
  #'
  #' @note
  #' I am a note
  #'
  #'
}

test.LM.Robust = function(e,tvObj,testOrder){
  this <- tvObj

  if(testOrder <= 0){
    message("Cannot execute test with no alternate hypothesis. Please set a valid testOrder")
    return(NaN)
  }
  # 1. Calc derivatives of params dgdt = Tx1 or Tx4 or Tx7 or...
  #    NCOL(dgdt) increases with the order of TV function.
  dgdt <- .dg_dt(this)

  # 2. Calc derivatives of taylor pars (linearised component) under the null
  dgdt2 <- .dg_dt2(this@st,testOrder)

  g <- calculate_g(this)
  X <- dgdt/g

  # 3. Invert crossprod(X) to calculate SSR1
  Xinv <- NULL
  try(Xinv <- solve(crossprod(X)))
  if (is.null(Xinv)){
    message("error")
    rm(g,X,dgdt,dgdt2)
    return(NaN)
  }

  XXXX <- X%*%Xinv%*%t(X)
  Y <- as.matrix(dgdt2/g)
  W <- as.matrix(Y-XXXX%*%Y)

  #4. Regress 1 on (psi2-1)*w, and compute SSR
  psi2_1 <- as.vector(e^2/g - 1)
  X <- psi2_1*W  #psi2_1 must be a vector for this!!

  #5. Compute test statistic:
  Xinv <- NULL
  try(Xinv <- solve(crossprod(X)))
  if(is.null(Xinv)) {
    message("error")
    rm(psi2_1,g,W,Y,X,XXXXdgdt,dgdt2)
    return(NaN)
  }

  Result <- this@Tobs-sum(diag(this@Tobs)-(X%*%Xinv%*%t(X)))

  # Tidy up & release memory before returning:
  rm(this,psi2_1,g,W,Y,X,XXXX,Xinv,dgdt,dgdt2)

  # Return:
  Result

}

.test.misSpec.Robust =  function(z,r1,r2){

  rtn <- list()
  z2_1 <- z
  Tobs <- NROW(z)

  SSR0 <- t(z2_1) %*% z2_1

  # 2: regress z2_1 on r1~r2, get SSR
  X <- cbind(r1,r2) # T x ncol(r1)+ncol(r2)
  Y <- z2_1  # T x 1
  #tXX <- solve(t(X),X,tol=1e-10)
  #b <- tXX %*% t(X) %*% Y  # vector len=ncol(X)
  b <- solve(t(X) %*% X) %*% t(X) %*% Y  # vector len=ncol(X)
  resid <- Y - (X %*% b)   # T x 1
  SSR1 <- t(resid) %*% resid # 1x1
  # LM stat
  rtn$LM <- as.numeric(Tobs * (SSR0-SSR1)/SSR0)
  rtn$pVal <- as.numeric(pchisq(rtn$LM,df=NCOL(r2),lower.tail=FALSE))

  # ROBUST
  # 1 as above
  # 2 regress r2 on r1, get residual vectors
  resid <- matrix(0,Tobs,NCOL(r2))
  X <- r1
  for (i in 1:NCOL(r2)){
    Y <- r2[,i,drop=FALSE]
    #tXX <- solve(t(X),X,tol=1e-10)
    #b <- tXX %*% t(X) %*% Y  # vector len=ncol(X)
    b <- solve(t(X) %*% X) %*% t(X) %*% Y  # vector len=ncol(X)
    resid[,i] <- Y-(X %*% b)   # T x 1
  }

  # regress 1 on (z2_1)resid, get SSR
  Y <- matrix(1,Tobs,1)
  X <- as.vector(z2_1) * resid
  #tXX <- solve(t(X),X,tol=1e-10)
  #b <- tXX %*% t(X) %*% Y  # vector len=ncol(X)
  b <- solve(t(X) %*% X) %*% t(X) %*% Y  # vector len=ncol(X)
  resid <- Y - (X %*% b)   # T x 1
  SSR <- t(resid) %*% resid # 1x1
  # LM Robust
  rtn$LMrob <- as.numeric(Tobs-SSR)
  rtn$pValrob <- as.numeric(pchisq(rtn$LMrob,df=NCOL(r2),lower.tail=FALSE))

  return(rtn)
}

getTestStats = function(e,tvObj,testOrder){
  this <- list()
  this$TR2 <- test.LM.TR2(e,tvObj,testOrder)
  this$Robust <- test.LM.Robust(e,tvObj,testOrder)
  return(this)
}



{# Help Doco: testStatDist() ####
  #' @noRd
#' @title
  #' Generates a test statistic distribution...
  #'
  #' @description
  #' `testStatDist` is a
  #'
  #' @usage testStatDist(refdata,tvObj,reftests,simcontrol)
  #'
  #' @param refdata A matrix
  #' @param tvObj Use Enum, garchtype$...
  #' @param reftests An integer
  #' @param simcontrol An integer
  #'
  #' @details
  #' This object
  #'
  #' ```
  #'  simcontrol$maxTestorder = 2
  #'  testResult = testStatDist(refdata,tvObj,reftests,simcontrol)
  #' ```
  #'
  #'
  #' @returns A list containing
  #'
  #' @note
  #' I am a note
  #'
  #'
}

testStatDist = function(refdata,tvObj,reftests,simcontrol){
  this <- tvObj

  if(is.null(simcontrol)){
    simControl <- list()
    simControl$saveAs <- paste("TestStatDist-",strftime(Sys.time(),format="%Y%m%d-%H%M%S",usetz = FALSE))
    simControl$numLoops <- 1100
    simControl$numCores <- parallel::detectCores() - 1
  }

  # Validate SimControl:
  if(is.null(simcontrol$maxTestorder)){
    warning("Maximum Test Order must be a valid number between 1 - 4")
    return(list())
  }else if((simcontrol$maxTestorder < 1) || (simcontrol$maxTestorder > 4) ){
    warning("Maximum Test Order must be a valid number between 1 - 4")
    return(list())
  }else if(length(reftests) != simcontrol$maxTestorder) {
    warning("Maximum Test Order is mis-matched with the number of Reference Tests provided")
    return(list())
  }

  # 1. Setup the default params


  if(!is.null(simcontrol$numLoops)) numLoops <- simcontrol$numLoops else numLoops <- 1100
  if(!is.null(simcontrol$numCores)) numCores <- simcontrol$numCores else numCores <- detectCores() - 1

  # 2. Load the generated data with Garch and add the 'g' from our TV object
  refdata <- refdata[1:this@Tobs,]*sqrt(this@g)

  # 3. Setup the parallel backend
  Sys.setenv("MC_CORES" = numCores)
  cl <- makeCluster(numCores)
  registerDoParallel(cl, cores = numCores)

  # 4. Set the estimation controls to suppress SE & console output
  estCtrl <- list(calcSE = FALSE, verbose = FALSE)

  # 5. Setup the timer to provide duration feedback
  tmr <- proc.time()
  timestamp(prefix = "Starting to build Test Stat Distribution - ",suffix = "\nPlease be patient as this may take a while...\n")

  # 6. Calculate Results for all test orders required (optimised for parallel processing)

  # First estimate all TV objects in parallel

  noGarch <- garch(garchtype$noGarch,1)  # Prevent the estimateTV method from creating a Garch object every time
  listTV <- foreach(a = 1:numLoops, .inorder=FALSE, .verbose = FALSE) %dopar%{
    estimateTV(refdata[,a],this,estCtrl,noGarch)    # Note: The tv params don't change, only the data changes
  }

  testStats <- foreach(b = 1:numLoops, .inorder=FALSE, .combine=rbind, .verbose = FALSE) %dopar% {
    if (isFALSE(listTV[[b]]$Estimated$error)) {

      foreach(testOrder = 1:simcontrol$maxTestorder, .inorder=FALSE, .combine=c, .verbose = FALSE) %do% {
        if(is.nan(reftests[[testOrder]]$TR2)) simTEST1 <- NA else simTEST1 <- test.LM.TR2(refdata[,b],listTV[[b]],testOrder)
        if(is.nan(reftests[[testOrder]]$Robust)) simTEST2 <- NA else simTEST2 <- test.LM.Robust(refdata[,b],listTV[[b]],testOrder)
        c(b,reftests[[testOrder]]$TR2,simTEST1,as.integer(simTEST1 > reftests[[testOrder]]$TR2),reftests[[testOrder]]$Robust,simTEST2,as.integer(simTEST2 > reftests[[testOrder]]$Robust),listTV[[b]]$Estimated$value)
        #runSimrow <- c(runSimrow,b,reftests[[testOrder]]$TR2,simTEST1,as.integer(simTEST1 > reftests[[testOrder]]$TR2),reftests[[testOrder]]$Robust,simTEST2,as.integer(simTEST2 > reftests[[testOrder]]$Robust),listTV[[b]]$Estimated$value)
      }

    }
  }


  # Debug: Internal Parallel Result:
  #  runSimrow
  # End: testStats <- foreach(b = 1:numloops,...

  # Extract Test P_Values from Results & express as decimal
  colnamesResults <- vector("character")
  if(simcontrol$maxTestorder >= 1) colnamesResults <- c(colnamesResults,"TestOrd_1","Ref$LMTR2","Stat_TR2.1","Pval_TR2.1","Ref$LMRobust","Stat_Robust.1","Pval_Robust.1","Estimated_LL")
  if(simcontrol$maxTestorder >= 2) colnamesResults <- c(colnamesResults,"TestOrd_2","Ref$LMTR2","Stat_TR2.2","Pval_TR2.2","Ref$LMRobust","Stat_Robust.2","Pval_Robust.2","Estimated_LL")
  if(simcontrol$maxTestorder >= 3) colnamesResults <- c(colnamesResults,"TestOrd_3","Ref$LMTR2","Stat_TR2.3","Pval_TR2.3","Ref$LMRobust","Stat_Robust.3","Pval_Robust.3","Estimated_LL")
  if(simcontrol$maxTestorder >= 4) colnamesResults <- c(colnamesResults,"TestOrd_4","Ref$LMTR2","Stat_TR2.4","Pval_TR2.4","Ref$LMRobust","Stat_Robust.4","Pval_Robust.4","Estimated_LL")
  #
  colnames(testStats) <- colnamesResults

  Results <- list()
  if(simcontrol$maxTestorder >= 1){
    Results$Order1 <- list()
    Results$Order1$pVal_TR2 <- mean(testStats[,"Pval_TR2.1"],na.rm = TRUE)
    Results$Order1$pVal_ROB <- mean(testStats[,"Pval_Robust.1"],na.rm = TRUE)
  }
  #
  if(simcontrol$maxTestorder >= 2){
    Results$Order2 <- list()
    Results$Order2$pVal_TR2 <- mean(testStats[,"Pval_TR2.2"],na.rm = TRUE)
    Results$Order2$pVal_ROB <- mean(testStats[,"Pval_Robust.2"],na.rm = TRUE)
  }
  #
  if(simcontrol$maxTestorder >= 3){
    Results$Order3 <- list()
    Results$Order3$pVal_TR2 <- mean(testStats[,"Pval_TR2.3"],na.rm = TRUE)
    Results$Order3$pVal_ROB <- mean(testStats[,"Pval_Robust.3"],na.rm = TRUE)
  }
  #
  if(simcontrol$maxTestorder >= 4){
    Results$Order4 <- list()
    Results$Order4$pVal_TR2 <- mean(testStats[,"Pval_TR2.4"],na.rm = TRUE)
    Results$Order4$pVal_ROB <- mean(testStats[,"Pval_Robust.4"],na.rm = TRUE)
  }
  Results$TestStatDist <- testStats

  # 7. Save the distribution
  if(!is.null(simcontrol$saveAs)) {
    # Create SimDist folder (if not there) & set Save filename
    if (!dir.exists(file.path(getwd(),"SimDist"))) dir.create(file.path(getwd(),"SimDist"))
    saveAs <- file.path("SimDist",simcontrol$saveAs)
    try(saveRDS(Results,saveAs))
  }

  # 8. Stop the parallel cluster
  stopCluster(cl)
  rm(cl)

  # 9. Print the time taken to the console:
  cat("\nTest Stat Distribution Completed \nRuntime:",(proc.time()-tmr)[3],"seconds\n")

  # 10. Attempt to release memory:
  rm(refdata,listTV,testStats)

  # Return:
  return(Results)

}


# TODO: Rename to setInitTVPars and change param to list()


.dg_dt =  function(tvObj){
  this <- tvObj

  rtn <- matrix(nrow=this@Tobs,ncol=this@nr.pars)
  col_idx <- 0

  if(isTRUE(this@delta0free)){
    col_idx <- col_idx + 1
    rtn[,col_idx] <- 1  # derivative of delta0
  }

  if (this@nr.transitions > 0) {
    # initialise some variables
    stdev_st <- sd(this@st)
    st_c <- speed_transf <- Gi <- 0

    for (i in 1:this@nr.transitions) {

      if(this$shape[i] == tvshape$single) st_c <- this@st - this$Estimated$pars["locN1",i]
      if(this$shape[i] == tvshape$double) st_c <- (this@st - this$Estimated$pars["locN1",i]) * (this@st - this$Estimated$pars["locN2",i])
      if(this$shape[i] == tvshape$double1loc) st_c <- (this@st - this$Estimated$pars["locN1",i])^2

      if(this$speedopt == speedopt$gamma) {
        speed_transf <- this$Estimated$pars["speedN",i]
        Gi <- 1/(1+exp(-this$Estimated$pars["speedN",i] * st_c))
      }
      if(this$speedopt == speedopt$gamma_std) {
        speed_transf <- this$Estimated$pars["speedN",i]/stdev_st
        Gi <- 1/(1+exp(-this$Estimated$pars["speedN",i] * st_c/stdev_st))
      }
      if(this$speedopt == speedopt$eta) {
        speed_transf <- exp(this$Estimated$pars["speedN",i])
        Gi <- 1/(1+exp(-exp(this$Estimated$pars["speedN",i]) * st_c))
      }

      deriv_const <- this$Estimated$pars["deltaN",i]*speed_transf*Gi*(1-Gi)

      col_idx <- col_idx + 1
      rtn[,col_idx] <- Gi    # derivative of delta1..n
      col_idx <- col_idx + 1
      rtn[,col_idx] <- deriv_const*st_c    # derivative of speed1..n

      if(this$shape[i] == tvshape$single){
        col_idx <- col_idx + 1
        rtn[,col_idx] <- -deriv_const    # derivative of loc1..n (shape=TVshape$single)
      }
      if(this$shape[i] == tvshape$double){
        col_idx <- col_idx + 1
        rtn[,col_idx] <- -deriv_const*(this@st-this$Estimated$pars["locN1",i])  # derivative of loc1..n (shape=TVshape$double)
        col_idx <- col_idx + 1
        rtn[,col_idx] <- -deriv_const*(this@st-this$Estimated$pars["locN2",i])  # derivative of loc2..n (shape=TVshape$double)
      }
      if(this$shape[i] == tvshape$double1loc){
        col_idx <- col_idx + 1
        rtn[,col_idx] <- -deriv_const*2*(this@st-this$Estimated$pars["locN1",i])    # derivative of loc1..n (shape=TVshape$double1loc)
      }

    } # End: for loop

  } # End: if (this@nr.transitions > 0)

  return(rtn)

}

.dg_dt2 =  function(st,testOrder){

  rtn <- matrix(nrow=NROW(st),ncol=testOrder)
  for(n in 1:testOrder){
    rtn[,n] <- st^n
  }
  # Return:
  return(rtn)
}










.dh_dg =  function(e,tvgarchObj){
  this <- tvgarchObj
  Tobs <- this@Tobs
  w <- e/sqrt(this$Estimated$tv$g)

  dhdg <- matrix(NA,nrow = Tobs,ncol = this@tvObj@nr.pars)  # Initialise return matrix = T x TV@nr.pars, value = 1

  if (this$garch$type == garchtype$noGarch){
    return(dhdg)
  } else if (this$garch$type == garchtype$general){
    v_garch <- cbind(c(0,rep(1,(Tobs-1))),c(0,w[1:(Tobs-1)]^2),c(0,this$Estimated$garch$h[1:(Tobs-1)])) # T x Num_garch_pars, each row = "1~w(i,t-1)^2~h(i,t-1)", i=1,...,N
  } else if (this$garch$type == garchtype$gjr){
    v_garch <- cbind(c(0,rep(1,(Tobs-1))),c(0,w[1:(Tobs-1)]^2),c(0,this$Estimated$garch$h[1:(Tobs-1)]),c(0,(min(w[1:(Tobs-1)],0))^2)) # T x Num_garch_pars, each row = "1~w(i,t-1)^2~h(i,t-1)", i=1,...,N
  }

  beta <- rep(this$Estimated$garch$Estimated$pars["beta",1],NCOL(v_garch))
  dhdg <- .ar1.Filter(v_garch,beta) # T x Num_garch_pars, each row = dh(t).dgarchpars
  return(dhdg)

}


.dh_dt =  function(e,tvgarchObj){
  this <- tvgarchObj
  Tobs <- this@tvObj@Tobs
  w <- e/sqrt(this$Estimated$tv$g)

  dhdt <- matrix(NA,nrow = Tobs,ncol = this@tvObj@nr.pars)  # Initialise return matrix = T x TV@nr.pars, value = 1

  dgdt <- .dg_dt(this$Estimated$tv)  # T x TV@nr.pars (includes d0's derivative if it is a free param)
  dgdt_l1 <- rbind(rep(0,NCOL(dgdt)),dgdt[(1:Tobs-1),])
  g <- this$Estimated$tv$g

  if (this$garch$type != garchtype$noGarch){
    v_tv <- (-this$Estimated$garch$Estimated$pars["alpha",1] * c(0,1/g[1:(Tobs-1)])*c(0,w[1:(Tobs-1)]^2)) * dgdt_l1
    beta <- rep(this$Estimated$garch$Estimated$pars["beta",1],NCOL(v_tv))
    dhdt <- .ar1.Filter(v_tv,beta) # T x Num_tv_pars, each row = dh(t).dtvpar
  }
  return(dhdt)
}

.df_dli =  function(tvObj,testOrder){
  this <- tvObj

  st <- this@st

  if(isTRUE(this@delta0free)){
    # We already have a derivative of the constant
    nrDerivs <- testOrder
    ret <- matrix(nrow=NROW(st),ncol=nrDerivs)
    for(n in 1:nrDerivs){ ret[,n] <- st^(n) }

  }else{

    nrDerivs <- testOrder + 1
    ret <- matrix(nrow=NROW(st),ncol=nrDerivs)
    for(n in 1:nrDerivs){ ret[,n] <- st^(n-1) }
  }
  return(ret)
}




# -- test.misSpec1(e,tvgarch,testOrd) ####
#' @noRd
#' @title
#' Tests for mis-specification of...
#'
#' @description
#' `test.misSpec1` is a
#'
#' @usage test.misSpec1(e,tvgarch,testOrd)
#'
#' @param e A numeric vector
#' @param tvgarch An estimated TvGarch model
#' @param testOrd Use Enum, garchtype$...
#'
#' @details
#' This object
#'
#'
#'
#' @returns A list...
#'
#' @note
#' I am a note
#'
#'
## Tests for mis-specification of the tv component of the model, e.g. is there another transition?
test.misSpec1 =  function(e,tvgarchObj,testOrder){
  this <- tvgarchObj
  rtn <- list()

  Tobs <- this@Tobs
  g <- this$Estimated$tv$g
  h <- this$Estimated$garch$h

  # derivatives:
  dg_dtv <- .dg_dt(this$Estimated$tv)        # T x nr.tv.pars
  dh_dtv <- .dh_dt(e,this)                   # T x nr.tv.pars
  dh_dga <- .dh_dg(e,this)                   # T x nr.garch.pars
  df_dli <- .df_dli(this@tvObj,testOrder)    # T x (testorder+1)

  z <- e/sqrt(g*h)
  z2 <- z^2
  z2_1 <- z^2 - 1

  # matrices
  u <- g^(-1) * dg_dtv     # (1/g)*(dg/dtvpars); T x nr.tv.pars
  x1 <- (h^(-1)) * dh_dtv  # (1/ht)*(dh/dtvpars); T x nr.tv.pars
  x2 <- (h^(-1)) * dh_dga  # (1/ht)*(dh/dgarchpars); T x nr.garch.pars
  #
  if(this$Estimated$garch$type == garchtype$noGarch){
    r1 <- u     # T x nr.tv.pars
  }else{
    r1 <- cbind(u+x1,x2)     # T x (tot.tv+garch.pars)
  }
  r2 <- (g^(-1)) * df_dli  # (1/gt)*(df/dlinpars); T x nr.lin.pars (=testorder+1)
  rtn <- .test.misSpec.Robust(z2_1,r1,r2)
  return(rtn)
}


# -- test.misSpec2(e,tvgarch,type) ####
#' @noRd
#' @title
#' Tests for mis-specification of...
#'
#' @description
#' `test.misSpec2` is a
#'
#' @usage test.misSpec2(e,tvgarch,type)
#'
#' @param e A numeric vector
#' @param tvgarch An estimated TvGarch model
#' @param type Use Enum, garchtype$...
#'
#' @details
#' This object
#'
#'
#'
#' @returns A list...
#'
#' @note
#' I am a note
#'
#'
## Tests for mis-specification of the Garch-order of the model

test.misSpec2 =  function(e,tvgarchObj,type){
  # type 1: GARCH(1,1) vs GARCH(1,2)
  # type 2: GARCH(1,1) vs GARCH(2,1)
  # type 3: GARCH vs STGARCH

  this <- tvgarchObj
  rtn <- list()

  Tobs <- this@Tobs
  g <- this$Estimated$tv$g
  h <- this$Estimated$garch$h
  w <- e/sqrt(g)
  z <- e/sqrt(g*h)
  z2 <- z^2
  z2_1 <- z^2 - 1

  # derivatives:
  dg_dtv <- .dg_dt(this$Estimated$tv)     # T x nr.tv.pars
  dh_dtv <- .dh_dt(e,this)                # T x nr.tv.pars
  dh_dga <- .dh_dg(e,this)                # T x nr.garch.pars

  # matrices
  u <- g^(-1) * dg_dtv     # (1/g)*(dg/dtvpars); T x nr.tv.pars
  x1 <- (h^(-1)) * dh_dtv  # (1/ht)*(dh/dtvpars); T x nr.tv.pars
  x2 <- (h^(-1)) * dh_dga  # (1/ht)*(dh/dgarchpars); T x nr.garch.pars

  r1 <- cbind(u+x1,x2)     # T x (tot.tv+garch.pars)
  if (type==1){ r2 <- (h^(-1)) * lag0(w^2,2)  }                 # T x 1
  if (type==2){ r2 <- (h^(-1)) * lag0(h,2)    }                 # T x 1
  if (type==3){ r2 <- (h^(-1)) * cbind(lag0(w,1),lag0(w^3,1)) } # T x 2

  rtn <- .test.misSpec.Robust(z2_1,r1,r2)
  return(rtn)
}


# -- test.misSpec3(e,tvgarch,maxLag) ####
#' @noRd
#' @title
#' Tests for mis-specification of...
#'
#' @description
#' `test.misSpec3` is a
#'
#' @usage test.misSpec3(e,tvgarch,maxLag)
#'
#' @param e A numeric vector
#' @param tvgarch An estimated TvGarch model
#' @param maxLag Use Enum, garchtype$...
#'
#' @details
#' This object
#'
#'
#'
#' @returns A list...
#'
#' @note
#' I am a note
#'
#'
## Tests for mis-specification of the remaining-ARCH in the model

test.misSpec3 =  function(e,tvgarchObj,maxLag){
  this <- tvgarchObj
  rtn <- list()

  Tobs <- this@tvObj@Tobs
  g <- this$Estimated$tv$g
  h <- this$Estimated$garch$h
  z <- e/sqrt(g*h)
  z2 <- z^2
  z2_1 <- z^2 - 1

  # derivatives:
  dg_dtv <- .dg_dt(this$Estimated$tv)     # T x nr.tv.pars
  dh_dtv <- .dh_dt(e,this)                # T x nr.tv.pars
  dh_dga <- .dh_dg(e,this)                # T x nr.garch.pars

  # matrices
  u <- g^(-1) * dg_dtv         # (1/g)*(dg/dtvpars); T x nr.tv.pars
  x1 <- (h^(-1)) * dh_dtv      # (1/ht)*(dh/dtvpars); T x nr.tv.pars
  x2 <- (h^(-1)) * dh_dga      # (1/ht)*(dh/dgarchpars); T x nr.garch.pars
  #
  if(this$Estimated$garch$type == garchtype$noGarch){
    r1 <- u     # T x nr.tv.pars
  }else{
    r1 <- cbind(u+x1,x2)     # T x (tot.tv+garch.pars)
  }
  #
  r2 <- lag0(z2,(1:maxLag))    # T x maxlag

  rtn <- .test.misSpec.Robust(z2_1,r1,r2)
  return(rtn)
}


# -- test.misSpec4(e,tvgarch,exogVar) ####
#' @noRd
#' @title
#' Tests for mis-specification of...
#'
#' @description
#' `test.misSpec4` is a
#'
#' @usage test.misSpec4(e,tvgarch,exogVar)
#'
#' @param e A numeric vector
#' @param tvgarch An estimated TvGarch model
#' @param exogVar Use Enum, garchtype$...
#'
#' @details
#' This object
#'
#'
#'
#' @returns A list...
#'
#' @note
#' I am a note
#'
#'
## Tests for mis-specification of the additive exogenous variable(s) in the model

test.misSpec4 =  function(e,tvgarchObj,exogVar){
  this <- tvgarchObj
  rtn <- list()

  Tobs <- this@Tobs
  g <- this$Estimated$tv$g
  h <- this$Estimated$garch$h

  # derivatives:
  dg_dtv <- .dg_dt(this$Estimated$tv)        # T x nr.tv.pars
  dh_dtv <- .dh_dt(e,this)                   # T x nr.tv.pars
  dh_dga <- .dh_dg(e,this)                   # T x nr.garch.pars

  z <- e/sqrt(g*h)
  z2 <- z^2
  z2_1 <- z^2 - 1

  # matrices
  u <- g^(-1) * dg_dtv     # (1/g)*(dg/dtvpars); T x nr.tv.pars
  x1 <- (h^(-1)) * dh_dtv  # (1/ht)*(dh/dtvpars); T x nr.tv.pars
  x2 <- (h^(-1)) * dh_dga  # (1/ht)*(dh/dgarchpars); T x nr.garch.pars
  #
  if(this$Estimated$garch$type == garchtype$noGarch){
    r1 <- u     # T x nr.tv.pars
  }else{
    r1 <- cbind(u+x1,x2)     # T x (tot.tv+garch.pars)
  }
  #
  r2 <- (g^(-1)) * exogVar  # (1/gt)*(exogVar); T x nr. Exogenous Vars

  rtn <- .test.misSpec.Robust(z2_1,r1,r2)
  return(rtn)
}








