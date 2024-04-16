## --- Multivariate N-TV-GARCH List Class --- ####
##


## --- ntvgarch_class Definition --- ####
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

## === Constructor: ntvgarch  ===####
#' @title
#' A list of estimated mtvgarch models
#'
#' @description
#' `ntvgarch` is a convenience object designed to hold a list of estimated mtvgarch models, to be used for multivariate correlation modelling
#'
#' @usage ntvgarch("tvgarch_list","series.names")
#'
#' @param tvgarch_list A standard R list() containing 2 or more univariate estimated mtvgarch models
#' @param series.names Optional. A character vector containing the names of the series in order.
#'
#' @details
#' This object was created to avoid issues with simple named lists.  By forcing the user to create this object, we can
#' ensure that all relevant validations are run and that the multivariate correlation models that use this will function as expected
#'
#' ```
#'   myNTVG = ntvgarch("tvgarch_list","series.names")
#' ```
#'
#'
#' @returns An ntvgarch_class object. Basically a named list with some useful attributes added.
#'
#' @note
#' I am a note
#'
#'
setGeneric(name="ntvgarch",
           valueClass = "ntvgarch_class",
           signature = c("tvgarch_list","series.names"),
           def = function(tvgarch_list,series.names){
             this <- new("ntvgarch_class")

             # Do basic validation checks:
             if(length(tvgarch_list) < 2) {
               warning("multivariate ntvgarch objects require a minimum of 2 data series.")
               return(this)
             }
             if(length(tvgarch_list) != length(series.names) ) {
               warning("The number of series.names provided must match the number of data series in the list.")
               return( new("ntvgarch_class") )
             }
             # End validation

             ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####
             ##
             ##    It is critical that the index of this object == index of tvgarch_list    ##
             ##
             ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####   ####

             for(n in 1:length(tvgarch_list)){
               # Validate the class_type - exit on error
               objType <- class(tvgarch_list[[n]])
               if(objType[1] != "tvgarch_class") {
                 warning("multivariate ntvgarch objects require the tvgarch_list param to be a list of tvgarch_class objects")
                 return(this)
               }
              # Add the tv_list items
              this[[series.names[n]]] <- tvgarch_list[[n]]
             }

             # Now add all the other Object attributes
             this$N <- length(tvgarch_list)
             this$Tobs <- tvgarch_list[[1]]@Tobs    # TODO: Validate all series are same length!!

             this$z <- this$w <- this$e <- matrix(NA,nrow = this$Tobs, ncol = this$N)
             this$h <- this$g <- matrix(1,nrow = this$Tobs, ncol = this$N)
             this$beta <- matrix(1,nrow = 1, ncol = this$N)
             #this$nr.tv.pars <- matrix(NA,nrow = 1, ncol = this$N)
             #this$nr.garch.pars <- matrix(NA,nrow = 1, ncol = this$N)


             for(n in 1:this$N){
               # Get the data from all the univar objects
               this$e[,n] <- tvgarch_list[[n]]@e
               # Filter the data & save it
               this$w[,n] <- this$e[,n]/sqrt(tvgarch_list[[n]]$Estimated$tv$g)
               this$z[,n] <- this$w[,n]/sqrt(tvgarch_list[[n]]$Estimated$garch$h)

               # Get the g & h data from the univar objects and compile into matrices (defaults for noGarch set above == 1)
               this$g[,n] <- tvgarch_list[[n]]$Estimated$tv$g
               if (tvgarch_list[[n]]$garchtype!=garchtype$noGarch) this$h[,n] <- tvgarch_list[[n]]$Estimated$garch$h
               if (tvgarch_list[[n]]$garchtype!=garchtype$noGarch) this$beta[1,n] <- tvgarch_list[[n]]$Estimated$garch$pars["beta",1]
               #this$nr.tv.pars[,n] <- tvgarch_list[[n]]@tvObj@nr.pars
               #this$nr.garch.pars[,n] <- tvgarch_list[[n]]@garchObj@nr.pars

             }

             return(this)
           }
)

setMethod("ntvgarch", signature = c("list","missing"),
          function(tvgarch_list){
            series.names = as.character(seq(1:length(tvgarch_list)))
            return(ntvgarch(tvgarch_list,series.names))
          }
)


