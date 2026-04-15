# *************************************************************************************************************** #
#                   Class generics for the garch_class, tv_class and tvgarch_class                ####
# *************************************************************************************************************** #

# --- GENERAL ----
setGeneric("summary")
setGeneric("plot")


# ---  TV --- SECTION  ####
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
                 this  <- .setInitTVPars(this)
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

{## estimateTV( - overloads - ) ----

  #TODO: Use noGarch instead of general - After fixing noGarch issues, @h, etc.

  setGeneric("estimateTV", valueClass = "tv_class",
             function(e, tvObj, garchObj, estimationControl) {
               standardGeneric("estimateTV")
             }
  )

}


# ---  GARCH --- SECTION  ####
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
               this <- .setInitGARCHPars(this)
               this$optimcontrol$ndeps <- rep(1e-3,this@nr.pars)
               if(type == garchtype$general) this$optimcontrol$parscale <- c(0.05,0.05,0.90)
               #if(type == garchtype$gjr) this$optimcontrol$parscale <- c(0.005, 0.095, 0.8, 0.1)
               return(this)
             }
  )
}

{## estimateGARCH( - overloads -)  ----

  setGeneric("estimateGARCH",
             valueClass = "garch_class",
             function(e, tvObj, garchObj, estimationControl) {
               standardGeneric("estimateGARCH")
             }
  )

}

{## estimateGARCH_RollingWindow( - overloads -)  ----
  setGeneric("estimateGARCH_RollingWindow",
             valueClass = "garch_class",
             function(e, garchObj, estimationControl) {
               standardGeneric("estimateGARCH_RollingWindow")
             }
  )

}


# ---  TVGARCH --- SECTION  ####

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

               # VarTargetting
               this <- setVarianceTargetting(this,1)  # Default = ON

               cat("\ntvgarch object created successfully!\n")

               return(this)
             }
  )
}

{## estimateTVGARCH() overrides ----

  setGeneric("estimateTVGARCH",
             valueClass = "tvgarch_class",
             function(e,tvgarchObj,estimationControl) {
               standardGeneric("estimateTVGARCH")
             }
  )

}

