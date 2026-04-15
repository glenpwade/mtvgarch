# *************************************************************************************************************** #
#                   Class methods for the garch_class                                                          ####
# *************************************************************************************************************** #


{## .setInitGARCHPars ----

  .setInitGARCHPars = function(garchObj){
    this <- garchObj

    if(this$type == garchtype$noGarch) {
      message("Cannot create Garch$pars for type: NoGarch")
      return(this)
    }

    maxLag <- max(this@order)

    # Set the row names:
    GarchparsRownames <- c("omega","alpha","beta","gamma")  #TODO: Remove or fix GJR (gamma)

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

    ## TODO: Implement in future release
    # if(this$type == garchtype$gjr) {
    #   this@nr.pars <- as.integer(4)
    #   pars <- matrix(nrow = this@nr.pars,ncol = maxLag)
    #   rownames(pars) <- GarchparsRownames[1:this@nr.pars]
    #   for(n in 1:maxLag){
    #     pars["omega",n] <- 0.005
    #     pars["alpha",n] <- 0.095
    #     pars["beta",n] <- 0.80
    #     pars["gamma",n] <- 0.10
    #   }
    #   this$pars <- pars
    # }
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
}

{## Constructor with 1 param - creates a Garch(1,1)  ----
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
          })
}

{## summary(garch) override ----
setMethod("summary",signature="garch_class",
      function(object,...){
        this <- object

        # Tabulate the model parameters, with 3 rows: (Spec|Estimated|StdError)
        tabSummary <- matrix(NaN, nrow=3, ncol=3)

        rownames(tabSummary) <- c("StartPar","Est.Par","StdErr")
        colnames(tabSummary) <- c("omega","alpha","beta")


        # Check if this object has been estimated - if not, return the specification & start pars:
        if(!("Estimated" %in% names(this)) ){
          loglikValue <- NA
          # Populate the table with specification only:
          # Transpose the pars matrix to display in row-format
          tabSummary[1,] <- this$pars[,1]

        }else{
          # Object has been estimated:
          loglikValue <- this$Estimated$value

          # Populate the table:

          # Transpose the pars matrix to display in row-format
          tabSummary[1,] <- this$pars[,1]
          tabSummary[2,] <- this$Estimated$pars[,1]

          # Check is StdErr was calculated:
          if("se" %in% names(this$Estimated) ) {
            # StdErr was calculated:

            tabSummary[3,] <- this$Estimated$se[,1]

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
        cat("\n| GARCH Estimation |  Loglik:",loglikValue)
        print(kable(tabSummary))
        cat("\nOmegaFree: ",this@omegafree,"\n")

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
    if(!("startParAdjust" %in% names(estCtrl))) { estCtrl$startParAdjust <- 10 }

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
    # boundSize <- estCtrl$startParAdjust *  this$optimcontrol$ndeps    # Should be in range 10 - 100
    # boundLo <- parsVec - boundSize
    # boundHi <- parsVec + boundSize

    # # Anna's method - apply 'shake' to par-scaled param, the re-scale:    # startParAdjust Should be in range 1 - 100, default=50
    scaled_parsVec <- parsVec/this$optimcontrol$parscale
    bound <- estCtrl$startParAdjust * this$optimcontrol$ndeps
    scaled_parsVec <- runif(this@nr.pars,(scaled_parsVec - bound),(scaled_parsVec + bound))
    parsVec <- scaled_parsVec * this$optimcontrol$parscale

    #RETURN:
    rtnList <- list()
    rtnList$Optimpars <- parsVec
    rtnList$garchObj <- this
    return(rtnList)

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

  .setEstimatedPars_GARCH <- function(garchObj,optimTmp){

    this <- garchObj
    tmp <- optimTmp

    if(tmp$convergence!=0) {
      this$Estimated$value <- NA
      this$Estimated$error <- TRUE
      this$Estimated$optimoutput <- tmp
      warning("estimateGARCH() - optim failed to converge. Check the optim controls & starting params")
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
    # TODO: Implement a better try-catch

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

{## -- estimateGARCH(e,TV,Garch,estCtrl) ----
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
      this$Estimated$value <- NA
      warning("estimateGARCH() generated a NULL return from optim()")
    }

    return(this)

  }
}


{## estimateGARCH( - overloads -)  ----

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


