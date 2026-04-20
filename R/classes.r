# *************************************************************************************************************** #
#                   Class definitions for the tv_class, garch_class and tvgarch_class                          ####
# *************************************************************************************************************** #

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

# *************************************************************************************************************** #
#                   Correlation Model classes                        ####
# *************************************************************************************************************** #

{## Multi-variate (N) model class ----
  ## --- Multivariate N-TV-GARCH List Class --- ##

  ## --- ntvgarch_class Definition --- ##
  ntvgarch <- setClass(Class = "ntvgarch_class",
                       contains = c("namedList")
  )

  ## === Initialise  ===####
  setMethod("initialize","ntvgarch_class",
            function(.Object,...){
              .Object <- callNextMethod(.Object,...)

              # Return:
              .Object
            })

}


{## CCC class ----
  ## -- The MTVGARCH package supports a number of Correlation objects
  ## -- This class file maintains the structure for CCC (Constant Conditional Correlation)

  ## Note:  The CTC (Constant Touplitz-Correlation)
  ##        Model can also be implemented using this class.

  ccc <- setClass(Class = "ccc_class",
                  slots = c(ntvg="ntvgarch_class",z="matrix"),
                  contains = c("namedList")
  )

  ## --- Initialise --- ##
  setMethod("initialize","ccc_class",
            function(.Object,...){
              .Object <- callNextMethod(.Object,...)

              # Default initial values
              .Object$N <- 0
              .Object$e <- matrix("numeric")
              .Object$Tobs <- 0
              .Object$nr.covPars <- 0
              # CCC$e will hold the data to be used by the model.
              .Object$P <- matrix()
              # Return:
              .Object
            })

}

{## STCC class ----
  ## -- The MTVGARCH package supports a number of Correlation objects
  ## -- This class file maintains the structure for STCC1 (STCC with One Transition) & STCC2 (STCC with Two Transitions)
  ## -- Both classes only support single transitions at this time


  # --- stcc1_class Definition --- #
  stcc1 <- setClass(Class = "stcc1_class",
                    slots = c(ntvg="ntvgarch_class",nr.corPars="integer",z="matrix"),
                    contains = c("namedList")
  )

  ## --- Initialise --- ##
  setMethod("initialize","stcc1_class",
            function(.Object,...){
              .Object <- callNextMethod(.Object,...)

              # Default initial values
              .Object$N <- 0
              .Object$e <- matrix("numeric")
              .Object$Tobs <- 0

              .Object$shape <- corrshape$single
              .Object$speedopt <- corrspeedopt$eta
              .Object$optimcontrol <- list(fnscale = -1, reltol = 1e-5, trace = 10)

              .Object$P1 <- matrix("numeric")
              .Object$P2 <- matrix("numeric")
              .Object$pars <- c(2.5,0.5)
              names(.Object$pars) <- c("speed","loc")


              # Return:
              .Object
            })

}

{## CDC class ----
  ## -- The MTVGARCH package supports a number of Correlation objects
  ## -- This class file maintains the structure for CDC (Constant Distance Correlation)

  ## Note:  The CDC model uses distance as a proxy for correlation

  cdc <- setClass(Class = "cdc_class",
                  slots = c(distData="data.frame",droppedObs="numeric",g="matrix",h="matrix",z="matrix"),
                  contains = c("namedList")
  )

  ## --- Initialise --- ##
  setMethod("initialize","cdc_class",
            function(.Object,...){
              .Object <- callNextMethod(.Object,...)

              # Default initial values
              .Object@distData <- data.frame()
              .Object$optimcontrol <- list(fnscale = -1, reltol = 1e-9)
              .Object$pars <- c(1)

              # Return:
              .Object
            }
  )


}

{## CEC class ----
  ## -- The MTVGARCH package supports a number of Correlation objects
  ## -- This class file maintains the structure for CEC (Constant Equi-Correlation)

  cec <- setClass(Class = "cec_class",
                  slots = c(N="integer",Tobs="integer",nr.covPars="integer"),
                  contains = c("namedList")
  )

  ## --- Initialise --- ##
  setMethod("initialize","cec_class",
            function(.Object,...){
              .Object <- callNextMethod(.Object,...)

              # Default initial values
              .Object@N <- as.integer(0)
              .Object@Tobs <- as.integer(0)
              .Object@nr.covPars <- as.integer(0)
              .Object$optimcontrol <- list(fnscale = -1, reltol = 1e-5)
              .Object$P <- matrix()

              .Object$pars <- c(1) #rho

              # Return:
              .Object
            }
  )
}

{## DCC class ----
  ## -- The MTVGARCH package supports a number of Correlation objects
  ## -- This class file maintains the structure for DCC (Dynamic Conditional Correlation)

  dcctype <- list(General=1,Dynamic=2)

  ## --- dcc_class Definition --- ##
  dcc <- setClass(Class = "dcc_class",
                  slots = c(st="numeric",nr.corPars="integer",nr.trPars="integer",Tobs="integer",N="integer",e="matrix"),
                  contains = c("namedList")
  )

  ## --- Initialise --- ##
  setMethod("initialize","dcc_class",
            function(.Object,...){
              .Object <- callNextMethod(.Object,...)

              # Default initial values
              .Object@N <- as.integer(0)
              .Object@Tobs <- as.integer(0)
              .Object@e <- matrix("numeric")
              .Object$ntvgarch <- list()
              .Object$speedopt <- corrspeedopt$eta
              .Object$optimcontrol <- list(fnscale = -1, reltol = 1e-5)
              .Object$pars <- c(0.05,0.80)
              names(.Object$pars) <- c("alpha","beta")

              # Return:
              .Object
            })
}
