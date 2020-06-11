## --- mtvgarch_class Structure --- ####

## --- When created:

# Slots (internal variables for use in methods - should only be set by pkg_code)

# mtvgarch@tvObj                     -- initial tv object used to create the multiplicative object
# mtvgarch@garchObj                  -- initial garch object used to create the multiplicative object
# mtvgarch@value                  -- "numeric" - starting log-liklihood value
# mtvgarch@e                      -- "numeric" - starting data set

# properties (external variables, visible to user)

# mtvgarch$Estimated$value        -- scalar: log-liklihood value
# mtvgarch$Estimated$tv
# mtvgarch$Estimated$garch
# mtvgarch$results[[1..n]]        -- "list" - contains a tv,garch & ll_value for each estimation iteration

## --- tv & garch class definitions are needed to create mtvgarch --- ####
tv <- setClass(Class = "tv_class",
               slots = c(st="numeric",g="numeric",delta0free="logical",nr.pars="integer", nr.transitions="integer",Tobs="integer",taylor.order="integer"),
               contains = c("namedList")
)
garch <- setClass(Class = "garch_class",
                  slots = c(h="numeric",nr.pars="integer",order="numeric"),
                  contains = c("namedList")
)
## --- MTVGARCH_CLASS Definition --- ####
mtvgarch <- setClass(Class = "mtvgarch_class",
               slots = c(tvObj="tv_class",garchObj="garch_class",value="numeric",e="numeric"),
               contains = c("namedList")
)

## Initialise with no params
setMethod("initialize","mtvgarch_class",
          function(.Object,...){
            .Object <- callNextMethod(.Object,...)
            .Object$results <- list()
            .Object$Estimated <- list()
            # Return:
            .Object
          })

setGeneric(name="mtvgarch",
           valueClass = "mtvgarch_class",
           signature = c("e","tvObj","garchObj"),
           def = function(e,tvObj,garchObj){

             this <- new("mtvgarch_class")

             # Validate: Spit dummy if TV is not estimated
             if(is.null(tvObj$Estimated) || is.null(garchObj$Estimated)) {
               message("mtvgarch-class objects require the tv & garch components to be estimated before initialising.")
               this$tv <- this$garch <- NULL
               return(this)
             }

             this@e <- e
             this@tvObj <- tvObj
             this@garchObj <- garchObj
             this@value <- loglik.mtvgarch.univar(e,garchObj@h,tvObj@g)
             this$initial_value <-  this@value

             # Configure the tv object, based on Garch type
             if(this@garchObj$type != garchtype$noGarch){
               this@tvObj@delta0free <- FALSE
               this@tvObj@nr.pars <- this@tvObj@nr.pars - as.integer(1)
               this@tvObj$optimcontrol$ndeps <- tail(this@tvObj$optimcontrol$ndeps,-1)
               this@tvObj$optimcontrol$parscale <- tail(this@tvObj$optimcontrol$parscale,-1)
             }
             cat("\nmtvgarch object created successfully!\n")
             cat("\nNext Steps:\n")
             cat("\n1. Copy the TV component from this object into a local TV variable, TV <- mtvgarch@tvObj")
             cat("\n2. Copy the GARCH component from this object into a local GARCH variable, GARCH <- mtvgarch@garchObj")
             cat("\n3. Filter the original data by dividing it by 'h',  e <- e/sqrt(GARCH@h) ")
             cat("\n4. Estimate the local TV using the filtered data, TV <- estimateTV(e,TV) ")
             cat("\n5. Filter the original data by dividing it by 'g' obtained in step 4.  e <- e/sqrt(TV@g) ")
             cat("\n6. Estimate the local GARCH using this filtered data, GARCH <- estimateGARCH(e,GARCH) ")
             cat("\n7. Make sure the estimated results have valid std errors.")
             cat("\n8. Calculate the Log-liklihood of the updated model: mtvgarch <- calcLogLik(mtvgarch,TV,GARCH) ")
             cat("\n9. Keep estimating TV & GARCH in pairs until the LogLiklihood is maximized ")

             return(this)
           }
)


## -- loglik.mtvgarch.univar() ####
setGeneric(name="loglik.mtvgarch.univar",
           valueClass = "numeric",
           signature = c("e","g","h"),
           def = function(e,g,h){
             sum( -0.5*log(2*pi) - 0.5*log(g) - 0.5*log(h) - 0.5*(e^2/(h*g) ) )
           }
)

## -- getTargetValue(e,TV,GARCH) ####
setGeneric(name="getTargetValue",
           valueClass = "numeric",
           signature = c("e","tvObj","garchObj"),
           def = function(e,tvObj,garchObj){
             ll <- NaN
             if(is.null(garchObj$Estimated)){
               g <- .calculate_g(tvObj)
               ll <- sum( -0.5*log(2*pi) - 0.5*log(g) - 0.5*(e^2)/g )
             } else {
               G1 <- calculate_h(garchObj,e)
               ll <- sum( -0.5*log(2*pi) - 0.5*log(G1@h) - 0.5*(e^2)/G1@h )
             }
             return(ll)
           }
)

setMethod("getTargetValue",signature = c("numeric","tv_class","missing"),
          function(e,tvObj){
            g1 <- garch(garchtype$noGarch)
            getTargetValue(e,tvObj,g1)
          })
setMethod("getTargetValue",signature = c("numeric","missing","garch_class"),
          function(e,garchObj){
            t1 <- new("tv_class")
            getTargetValue(e,t1,garchObj)
          })

## -- calcLoglik() ####
setGeneric(name="calcLoglik",
           valueClass = "mtvgarch_class",
           signature = c("mtvgarchObj","tvObj","garchObj"),
           def = function(mtvgarchObj,tvObj,garchObj){
             this <- mtvgarchObj

             NextResult <- length(this$results) + 1
             this$results[[NextResult]] <- list()
             this$results[[NextResult]]$tv <- tvObj
             this$results[[NextResult]]$garch <- garchObj

             ll_val <- loglik.mtvgarch.univar(this@e,tvObj@g,garchObj@h)
             this$results[[NextResult]]$value <- ll_val

             cat("\nThe LogLik value is:",ll_val)
             cat("\n\nIf this is better than any previous estimate then do another iteration,")
             cat("\notherwise, we are done.")

             return(this)
           }
)



