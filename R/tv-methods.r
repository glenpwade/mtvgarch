# *************************************************************************************************************** #
#                   Class methods for the tv_class                                                             ####
# *************************************************************************************************************** #

{## .setInitTVPars  ----
  .setInitTVPars = function(tvObj){
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

              # Tabulate the model parameters, with 3 rows: (Spec|Estimated|StdError)
              #  The first col is always delta0 (even if it is not an estimated par)
              tabSummary <- matrix(NaN, nrow=3, ncol=1+(4*this@nr.transitions) )

              rownames(tabSummary) <- c("StartPar","Est.Par","StdErr")
              colNames <- "delta0"
              if(this@nr.transitions > 0){
                for( n in 1:this@nr.transitions){
                  nextTrans <- c(paste0("d",n), paste0("spd",n), paste0("loc1.",n), paste0("loc2.",n))
                  colNames <- c(colNames, nextTrans )
                }
                colnames(tabSummary) <- colNames
              }

              # Check if this object has been estimated - if not, return the specification & start pars:
              if(!("Estimated" %in% names(this)) ){
                loglikValue <- NA

                # Populate the table with specification only:
                tabSummary[1,1] <- this$delta0
                # Transpose the pars matrix to display in row-format
                if(this@nr.transitions > 0){
                  for(n in 1:this@nr.transitions){
                    var <- 4*n
                    tabSummary[1,(var-2):(var+1)] <- this$pars[,n]
                  }
                }

              }else{
                # Object has been estimated:
                loglikValue <- this$Estimated$value

                # Populate the table:
                tabSummary[1,1] <- this$delta0
                tabSummary[2,1] <- this$Estimated$delta0
                # Check is StdErr was calculated:
                if("se" %in% names(this$Estimated) ) tabSummary[3,1] <- this$Estimated$delta0_se

                # Transpose the pars matrix to display in row-format
                if(this@nr.transitions > 0){
                  for(n in 1:this@nr.transitions){
                    var <- 4*n
                    tabSummary[1,(var-2):(var+1)] <- this$pars[,n]
                    tabSummary[2,(var-2):(var+1)] <- this$Estimated$pars[,n]
                    # Check is StdErr was calculated:
                    if("se" %in% names(this$Estimated) ) tabSummary[3,(var-2):(var+1)] <- this$Estimated$se[,n]
                  }
                }

                if("se" %in% names(this$Estimated) ) {
                  # StdErr was calculated:

                  # Add a significance indicator:
                  sigStars <- rep("",NCOL(tabSummary))  #Initialise the sigStars vector

                  for(n in 1:NCOL(tabSummary)){
                    # Calculate significance from estimated value & std.Err
                    if(  isTRUE((tabSummary[3,n]*2.576) < abs(tabSummary[2,n])) ) {
                      sigStars[n] <- "***"
                    }else{
                      if( isTRUE((tabSummary[3,n]*1.960) < abs(tabSummary[2,n])) ) {
                        sigStars[n] <- "**"
                      }else{
                        if( isTRUE((tabSummary[3,n]*1.645) < abs(tabSummary[2,n])) ) {sigStars[n] <- "*"}
                      }
                    }
                  }

                  # Round to default dp:
                  tabSummary <- round(tabSummary,getOption("digits"))

                  # Add significance Stars & braces to table:
                  for(n in 1:NCOL(tabSummary)){
                    tabSummary[2,n] <- paste0(tabSummary[2,n],sigStars[n])
                    tabSummary[3,n] <- paste0("(",tabSummary[3,n],")")
                  }


                }  # End: StdErr was calculated:


              }  # End: # Object has been estimated:

              # FINALLY: Print the Est Results as a Table
              cat("\n| TV Estimation |  Loglik:",loglikValue)
              print(kable(tabSummary))
              cat("\nDelta0Free: ",this@delta0free,"\n")

            })
}


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
    if(!("startParAdjust" %in% names(estCtrl))) { estCtrl$startParAdjust <- 10 }

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
      # boundSize <- estCtrl$startParAdjust * this$optimcontrol$ndeps    # startParAdjust Should be in range 1 - 10
      # boundLo <- parsVec - boundSize
      # boundHi <- parsVec + boundSize

      # # Anna's method - apply 'shake' to par-scaled param, the re-scale:    # startParAdjust Should be in range 1 - 100, default=50
      scaled_parsVec <- parsVec/this$optimcontrol$parscale
      bound <- estCtrl$startParAdjust * this$optimcontrol$ndeps
      scaled_parsVec <- runif(this@nr.pars,(scaled_parsVec - bound),(scaled_parsVec + bound))
      parsVec <- scaled_parsVec * this$optimcontrol$parscale

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


{##  estimateTV(e,tv,garch,ctrl) ----
  .estimateTV <- function(e,tvObj,garchObj,estimationControl){

    this <- tvObj
    # debug:
    # this <- TVspec
    # garchObj <- GARCHspec
    # estimationControl <- estCtrl

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

{## estimateTV( - overloads - ) ----

  #TODO: Use noGarch instead of general - After fixing noGarch issues, @h, etc.

  # TV <- estimateTV(e,tvSpec,garchSpec,estCtrl)
  setMethod("estimateTV",
            signature = c(e="numeric",tvObj="tv_class",garchObj="garch_class",estimationControl="list"), #note: All lower case
            function(e,tvObj,garchObj,estimationControl){
              # Ensure estimationControl has all elements:
              if(!("calcSE" %in% names(estimationControl)) ) estimationControl$calcSE <- FALSE
              if(!("verbose" %in% names(estimationControl)) ) estimationControl$verbose <- FALSE
              if(!("maxIter" %in% names(estimationControl)) ) estimationControl$verbose <- 1
              if(!("fixStartPars" %in% names(estimationControl)) ) estimationControl$fixStartPars <- FALSE
              if(!("startParAdjust" %in% names(estimationControl)) ) estimationControl$startParAdjust <- 10
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
