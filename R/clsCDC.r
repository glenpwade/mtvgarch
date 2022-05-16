## -- The MTVGARCH package supports a number of Correlation objects
## -- This class file maintains the structure for CDC (Constant Distance Correlation)

## Note:  The CDC model uses distance as a proxy for correlation

cdc <- setClass(Class = "cdc_class",
                  slots = c(distData="data.frame",droppedObs="numeric",g="matrix",h="matrix",z="matrix"),
                  contains = c("namedList")
)

## --- Initialise --- ####
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

## -- .calc.distCorr -- ####
setGeneric(name=".calc.distCorr",
           valueClass = "matrix",
           signature = c("pars","distData"),
           def =   function(pars,distData){

             # # Test:
             # pars = CDC$pars
             # distData = CDC@distData

             # Validate the number of beta-pars matches the distance vectors

             # Remove alpha from the pars vector, for convenience
             alpha = pars[1]
             betaPars <- pars[-1]  # drop alpha to leave just beta pars
             expArg = 0

             # Now calculate the exponent argument using all provided data
             for(n in seq_along(betaPars)) expArg <- expArg - betaPars[n]*distData[,n]
             #
             veclCorr <- alpha*exp(expArg)
             distCorr <- unVecL(veclCorr)

           return(distCorr)
           }
)

## -- Constructor:cdc -- ####
setGeneric(name="cdc",
           valueClass = "cdc_class",
           signature = c("ntvgarchObj","distanceData","pars"),
           def = function(ntvgarchObj,distanceData,pars){
             this <- new("cdc_class")

             if(!is.null(ntvgarchObj)){
               ## -- Do validation checks -- ##
               objType <- class(ntvgarchObj)
               if(objType[1] != "ntvgarch_class"){
                 warning("a valid instance of the ntvgarch_class is required to create a cdc model")
                 return(this)
               }
             }
             if(isFALSE(all.equal(NCOL(distanceData),(length(pars)-1))) ){
               warning("Size mis-match: distanceData, pars")
               return(this)
             }
             # Need to check that the distanceData represents N series, i.e. matches ntvgarch
             # End validation

             # Set Default Values:
             this$ntvgarch <- ntvgarchObj
             this@distData <- as.data.frame(distanceData)
             this$pars <- pars

             # Now build the scaled Distance-Correlation matrices

             # Step 1: Get the minimum common data length
             # Capture the series lengths into a vector
             vecTobs = vapply(this$ntvgarch,function(X) attributes(X)[["Tobs"]], vector("numeric",1),USE.NAMES = FALSE)
             # Find the minimum value
             commonLen <- min(vecTobs)
             # Calculate how many obs will be dropped
             this@droppedObs <- max(vecTobs) - min(vecTobs)

             # Step 2: Extract e, g, h from ntvgarch object & calculate z
             this@z <- this@g <- this@h <- e <- matrix(NA,commonLen,length(this$ntvgarch))
             for(n in seq_along(this$ntvgarch)){
               #TODO: Add support for garch@h (may need to update the tvgarch class)
               this@g[,n] <- this$ntvgarch[[n]]$Estimated$tv$g[1:commonLen]
               this@h[,n] <- rep(1,commonLen)
               e[,n] <- this$ntvgarch[[n]]@e[1:commonLen]
             }
             this@z <- e / sqrt(this@g)

             # Step 4: Create the pseudo-Corr matrix (starting values)
             this$P <- .calc.distCorr(this$pars, this@distData)

             # check PosDef:
             posDef <- FALSE
             try( { eig <- eigen(this$P,only.values = TRUE)
             posDef = isTRUE(min(eig$values) >= 0) }, silent = TRUE
             )
             if (isFALSE(posDef)) {
               warning("The parameters supplied did not generate a pos-def starting matrix. \nPlease try changing the gamma pars and/or check the distanceData. ")
             }

             return(this)
           }
)

## -- .loglik.cdc.ml() --####
setGeneric(name=".loglik.cdc.ml",
           valueClass = "numeric",
           signature = c("optimpars","cdcObj"),
           def = function(optimpars,cdcObj){

             err_output <- -1e20
             this <- cdcObj
             N = NCOL(this$P)

             #### ======== constraint checks ======== ####

             # # Check : Check the boundary values for gamma params:
             if(isTRUE( min(optimpars) <= 0) ) return(err_output)
             if(isTRUE( max(optimpars[1]) > 1) ) return(err_output)

             mP <- .calc.distCorr(optimpars,this@distData)

             eig <- NULL
             try( eig <- eigen(mP,symmetric=TRUE,only.values=TRUE) )
             if(is.null(eig)) return(err_output)

             # Check for SPD - positive-definite check:
             #if (min(eig$values) <= 0) return(err_output)
             if (isTRUE(min(eig$values) <= 0) ) {
               res <- nearPD(mP,corr = TRUE, maxit = 250)
               mP <- matrix(res$mat@x,N,N)
             }

             #### ======== calculate loglikelihood using Maximum Liklihood ======== ####
             Tobs <- NROW(this@z)
             llt <- vector("numeric",Tobs)
             mPinv <- solve(mP)
             # llConst:
             logDet_mP_2Pi <- log(det(mP)) + log(2*pi)
             for(t in seq.int(1,Tobs)) llt[t] <-  0.5 * ( -logDet_mP_2Pi -(t(this@z[t,]) %*% mPinv %*% this@z[t,]) -sum(log(this@g[t,])) -sum(log(this@h[t,])) )

             # Return:
             return(sum(llt,na.rm = TRUE))

           }
)


## -- .loglik.cdc.ls() -- ####
setGeneric(name=".loglik.cdc.ls",
           valueClass = "numeric",
           signature = c("optimpars","cdcObj"),
           def = function(optimpars,cdcObj){

             err_output <- 1e10
             this <- cdcObj

             #### ======== constraint checks ======== ####

             # # Check : Check the boundary values for gamma params:
             if(isTRUE( min(optimpars) <= 0) ) return(err_output)
             if(isTRUE( max(optimpars[1]) > 1) ) return(err_output)

             #### ======== calculate loglikelihood using Least-Squares ======== ####

             sampleCor <- cor(this@z)
             sampleVecL <- vecL(sampleCor)

             estCorr <- .calc.distCorr(optimpars,this@distData)
             estCorrVecL <- vecL(estCorr)

             ls <- (sampleVecL - estCorrVecL)^2

             # Return:
             return(sum(ls))

           }
)


## --- estimateCDC --- ####
setGeneric(name="estimateCDC",
           valueClass = "cdc_class",
           signature = c("cdcObj","estimationCtrl","estMethod"),
           def = function(cdcObj,estimationCtrl,estMethod){
             this <- cdcObj

             calcSE <- estimationCtrl$calcSE
             verbose <- estimationCtrl$verbose

             this$Estimated <- list()
             optimpars <- this$pars
             optimpars <- this$pars[!is.na(this$pars)]

             tmp <- NULL
             if(toupper(estMethod)=="ML"){
               # Maximum-Likelihood
               this$optimcontrol$fnscale = -1
               try(tmp <- optim(optimpars,.loglik.cdc.ml,this,gr=NULL,method="BFGS",control=this$optimcontrol,hessian=calcSE))
             }
             if(toupper(estMethod)=="LS"){
               # Minimise Least-Squared-Error
               this$optimcontrol$fnscale = 1
               try(tmp <- optim(optimpars,.loglik.cdc.ls,this,gr=NULL,method="BFGS",control=this$optimcontrol,hessian=calcSE))
             }

             if(!is.null(tmp)){
               # Check convergence
               if(tmp$convergence < 2){
                 this$Estimated$pars <- tmp$par
                 if(isTRUE(calcSE)) {
                   try(this$Estimated$se <- sqrt(diag(solve(-tmp$hessian))), silent = TRUE)
                 }
                 this$Estimated$P <- .calc.distCorr(tmp$par,this@distData)
                 if(toupper(estMethod)=="ML") this$Estimated$value <- .loglik.cdc.ml(tmp$par,this)
                 if(toupper(estMethod)=="LS") this$Estimated$value <- .loglik.cdc.ls(tmp$par,this)
                 if(isTRUE(verbose)) this$Estimated$tmp <- tmp
                 if(this@droppedObs > 0){
                   warnMsg <- paste0("Not all data series are the same length. ",this@droppedObs," observations have been dropped.")
                   message(warnMsg)
                 }
               } else this$Estimated$tmp <- tmp

             }

             return(this)
           }
)



