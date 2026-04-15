# *************************************************************************************************************** #
#                   Class definitions for the garch_class, tv_class and tvgarch_class                          ####
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



