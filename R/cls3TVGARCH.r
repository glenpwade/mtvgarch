## --- tvgarch_class Structure --- ####

## --- When created:

# Slots (internal variables for use in methods - should only be set by pkg_code)

# tvgarch@tvObj                     -- initial tv object used to create the multiplicative object
# tvgarch@garchObj                  -- initial garch object used to create the multiplicative object
# tvgarch@value                  -- "numeric" - starting log-liklihood value
# tvgarch@e                      -- "numeric" - starting data set

# properties (external variables, visible to user)

# tvgarch$Estimated$value        -- scalar: log-liklihood value
# tvgarch$Estimated$tv
# tvgarch$Estimated$garch
# tvgarch$results[[1..n]]        -- "list" - contains a tv,garch & ll_value for each estimation iteration


## --- tvgarch_CLASS Definition --- ####
tvgarch <- setClass(Class = "tvgarch_class",
               slots = c(tvObj="tv_class",garchObj="garch_class",value="numeric",e="numeric"),
               contains = c("namedList")
)

## Initialise with no params
setMethod("initialize","tvgarch_class",
          function(.Object,...){
            .Object <- callNextMethod(.Object,...)
            .Object$results <- list()
            .Object$Estimated <- list()
            # Return:
            .Object
          })

setGeneric(name="tvgarch",
           valueClass = "tvgarch_class",
           signature = c("e","tvObj","garchObj"),
           def = function(e,tvObj,garchObj){

             this <- new("tvgarch_class")

             # Validate: Spit dummy if TV is not estimated
             if(is.null(tvObj$Estimated) || is.null(garchObj$Estimated)) {
               message("tvgarch-class objects require the tv & garch components to be estimated before initialising.")
               this$tv <- this$garch <- NULL
               return(this)
             }

             this@e <- e
             this@tvObj <- tvObj
             this@garchObj <- garchObj
             this@value <- loglik.tvgarch.univar(e,garchObj@h,tvObj@g)
             this$initial_value <-  this@value

             # Configure the tv object, based on Garch type
             if(this@garchObj$type != garchtype$noGarch){
               this@tvObj@delta0free <- FALSE
               this@tvObj@nr.pars <- this@tvObj@nr.pars - as.integer(1)
               this@tvObj$optimcontrol$ndeps <- tail(this@tvObj$optimcontrol$ndeps,-1)
               this@tvObj$optimcontrol$parscale <- tail(this@tvObj$optimcontrol$parscale,-1)
             }
             cat("\ntvgarch object created successfully!\n")
             cat("\nNext Steps:\n")
             cat("\n1. Copy the TV component from this object into a local TV variable, TV <- tvgarch@tvObj")
             cat("\n2. Copy the GARCH component from this object into a local GARCH variable, GARCH <- tvgarch@garchObj")
             cat("\n3. Filter the original data by dividing it by 'h',  e <- e/sqrt(GARCH@h) ")
             cat("\n4. Estimate the local TV using the filtered data, TV <- estimateTV(e,TV) ")
             cat("\n5. Filter the original data by dividing it by 'g' obtained in step 4.  e <- e/sqrt(TV@g) ")
             cat("\n6. Estimate the local GARCH using this filtered data, GARCH <- estimateGARCH(e,GARCH) ")
             cat("\n7. Make sure the estimated results have valid std errors.")
             cat("\n8. Calculate the Log-liklihood of the updated model: tvgarch <- calcLogLik(tvgarch,TV,GARCH) ")
             cat("\n9. Keep estimating TV & GARCH in pairs until the LogLiklihood is maximized ")

             return(this)
           }
)


## -- loglik.tvgarch.univar() ####
setGeneric(name="loglik.tvgarch.univar",
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
           valueClass = "tvgarch_class",
           signature = c("tvgarchObj","tvObj","garchObj"),
           def = function(tvgarchObj,tvObj,garchObj){
             this <- tvgarchObj

             nextResult <- length(this$results) + 1
             this$results[[nextResult]] <- list()
             this$results[[nextResult]]$tv <- tvObj
             this$results[[nextResult]]$garch <- garchObj

             ll_val <- loglik.tvgarch.univar(this@e,tvObj@g,garchObj@h)
             this$results[[nextResult]]$value <- ll_val

             cat("\nThe LogLik value is:",ll_val)
             cat("\n\nIf this is better than any previous estimate then do another iteration,")
             cat("\notherwise, we are done.")

             return(this)
           }
)



