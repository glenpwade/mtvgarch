## -- The MTVGARCH package supports a number of Correlation objects
## -- This class file maintains the structure for CEC (Constant Equi-Correlation)

cec <- setClass(Class = "cec_class",
                slots = c(N="integer",Tobs="integer",nr.covPars="integer"),
                contains = c("namedList")
)

## --- Initialise --- ####
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

## -- Constructor:cec -- ####
setGeneric(name="cec",
           valueClass = "cec_class",
           signature = c("nr.series","ntvgarchObj"),
           def = function(nr.series,ntvgarchObj){
             this <- new("cec_class")
             if(!is.null(ntvgarchObj)){
               ## -- Do validation checks -- ##
               objType <- class(ntvgarchObj)
               if(objType[1] != "ntvgarch_class"){
                 warning("a valid instance of the ntvgarch_class is required to create a cec model")
                 return(this)
               }
               # End validation

               # Set Default Values:
               this@N <- ntvgarchObj@N
               this@Tobs <- ntvgarchObj@Tobs
               N <- this@N

               # Add the Estimated components from the ntvgarch
               this$ntvgarch <- list()
               for(n in 1:N){
                 this$ntvgarch[[n]] <- list()
                 this$ntvgarch[[n]]$tv <- ntvgarchObj[[n]]@tvObj
                 this$ntvgarch[[n]]$garch <- ntvgarchObj[[n]]@garchObj
               }
               names(this$ntvgarch) <- names(ntvgarchObj)
             }

             if(!is.null(nr.series)){
               N <- this@N <- as.integer(nr.series)
             }
             #TODO: Set default starting value for rho:

             this@nr.covPars <- as.integer((N^2-N)/2)
             this$P <- matrix(0.5,N,N)
             diag(this$P) <- 1
             message("cec object created. Default correlation is 0.5  \nDon't forget to set the starting-parameter $P matrix.")

             return(this)
           }
)


## --- estimateCEC  --- ####
setGeneric(name="estimateCEC",
           valueClass = "cec_class",
           signature = c("e","cecObj","estimationCtrl"),
           def = function(e,cecObj,estimationCtrl){
             this <- cecObj
             this$Estimated <- list()

             if(is.null(this$ntvgarch)){
               # Estimate if no ntvgarch provided..?
               # TODO: For CCC we just called cor() - what do we do for CEC?
               warning("For CCC we just called cor() - what do we do for CEC?")
             }else{
               # Estimate the CEC model, using the provided ntvgarch
               # 1. Validate model & data sizes
               if(isFALSE( all.equal(NCOL(e),length(this$ntvgarch)) ) ) {
                 warning("size mismatch: e, cecObj$ntvgarch")
                 return(this)
               }

               # 2. Standardise data
               #TODO: Replace for(loop) with apply() - check FilterData() function in clsNTVGarch
               z <- w <- e
               for(n in 1:this@N){
                 w[,n] <- e[,n]/sqrt(this$ntvgarch[[n]]$tv@g)
                 z[,n] <- w[,n]/sqrt(this$ntvgarch[[n]]$garch@h)
               }

               # 3. Call optim()
               optimpars <- this$pars
               try(tmp <- optim(optimpars,.loglik.cec,z,this,gr=NULL,method="BFGS",control=this$optimcontrol,hessian=estimationCtrl$calcSE))

               ## --- Attach results of estimation to the object --- ##

               # An unhandled error could result in a NULL being returned by optim()
               if (is.null(tmp)) {
                 this$Estimated$value <- NA
                 this$Estimated$error <- TRUE
                 warning("estimateCEC() - optim failed and returned NULL. Check the optim controls & starting params")
                 return(this)
               }
               if (tmp$convergence != 0) {
                 this$Estimated$value <- NA
                 this$Estimated$error <- TRUE
                 this$Estimated$optimoutput <- tmp
                 warning("estimateCEC() - failed to converge. Check the optim controls & starting params")
                 return(this)
               }

               this$Estimated$value <- tmp$value
               this$Estimated$error <- FALSE

               #Update the CEC object parameters using optimised pars:
               this$Estimated$pars <- tmp$par
               #colnames(this$Estimated$pars) <- "rho"

               # Calc the std errors
               if (isTRUE(estimationCtrl$calcSE)){

                 cat("\nCalculating CEC standard errors...\n")
                 this$Estimated$se <- NULL
                 stdErrors <- NULL
                 this$Estimated$hessian <- tmp$hessian
                 #colnames(this$Estimated$se) <- "rho_se"

                 try(stdErrors <- sqrt(-diag(qr.solve(tmp$hessian))))
                 if(!is.null(stdErrors)){
                     this$Estimated$se <- stdErrors
                   }else this$Estimated$se <- NaN

               if (isTRUE(estimationCtrl$verbose)) this$Estimated$optimoutput <- tmp
             }
             return(this)
             }
           }
)

.loglik.cec <- function(pars,z,cec){
  # input: pars        -- c(rho)
  #        z           -- volatility standardised returns (matrix TxN
  # model: CEC

  Tobs <- NROW(z)
  N <- NCOL(z)

  #Initialise return value to the Error:
  ll <- err_output <- -1e10

  # - - - CEC - - -
  rho <- pars
  Q <- .eigVec.EC(N) ## eigenvectors (deterministic)
  L <- .eigVal.EC(N,rho) ## eigenvalues for CEC
  P <- Q%*%diag(L)%*%t(Q)  ##
  if (min(L) <= 0) return(err_output)

  # - - - P and loglik-value
  llt <- rep(0,Tobs)
  Pinv <- Q%*%diag(1/L)%*%t(Q)
  detP <- prod(L)
  for(t in seq(1,Tobs))
  {
    llt[t] <- -0.5*log(detP)-0.5*( t(z[t,])%*%(Pinv)%*%z[t,])
  }

  ll <- sum(llt)
  return(ll)

} #End: myLogLik.cec




