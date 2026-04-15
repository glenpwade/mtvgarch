# *************************************************************************************************************** #
#                   Class methods for the tvgarch_class                                                        ####
# *************************************************************************************************************** #

# ---  TVGARCH --- SECTION  ####


{## setVarianceTargetting ----
  setVarianceTargetting <- function(tvgarchObj,On_Off){

    this <- tvgarchObj

    if(On_Off==1){
      # Variance Targetting ON:
      this$varTarget <- TRUE
      this@tvObj@delta0free <- TRUE
      this@garchObj@omegafree <- FALSE

    }else{
      # OFF:
      this$varTarget <- FALSE
      this@tvObj@delta0free <- FALSE
      this@garchObj@omegafree <- TRUE

    }
    return(this)
  }

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

{## Help doco:  estimateTVGARCH ####
#' @noRd
#' @title
#' Estimates a tv model
#'
#' @description
#' `estimateTVGARCH` is a
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
#'   myTvGarch = estimateTVGARCH(e,myTvGarch,estimationControl)
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
    if("maxIter" %in% names(estimationControl)) maxIterations <- estimationControl$maxIter else maxIterations <- 100
    if(!(is.numeric(this$iterationReltol) )) {this$iterationReltol <- 1e-5}

    # NOte: the fixStartPars and startParAdjust parameters will be passed down to the underlying Estimators

    # Configure variance targetting:  TODO: Test estimators with var Target OFF!
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

      # 2. Run estimations and calculate a new tvg.loglik
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

        # 3. Calc the loglik value ONLY when TV AND GARCH have successfully estimated
        if(class(TV)=="tv_class" && class(GARCH)=="garch_class"){
          # Confirm the classes are valid
          tvg.loglik <- loglik.tvgarch.univar(e,g=TV@g,h=GARCH@h)
        }else{
          # What should we do if one of the estimator failed - but we got this far??
        }
        #
        if(isTRUE(estimationControl$verbose)){cat(".")}  # Single line comment, Garch done

      })

      # Calculate the initial/starting Loglik.Value (Use to determine convergence/divergence)
      if(this@iterations == as.integer(1)){ startLoglik <- loglik.tvgarch.univar(e,g=TV@g,h=GARCH@h) }

      # 4. Increment the iteration counter
      this@iterations <- this@iterations + as.integer(1)

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

