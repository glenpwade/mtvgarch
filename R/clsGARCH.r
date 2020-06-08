## --- garch_class Structure --- ####

## --- When created:

# Slots (internal variables for use in methods - should only be set by pkg_code)

    # garch@h             -- numeric: conditional variances
    # garch@nr.pars       -- integer: number of parameters in model 
    # garch@order         -- numeric: Garch(1,1), Only first order (p=q=1) is implemented.

# properties (external variables, visible to user)

    # garch$type                  -- General & GJR supported
    # garch$pars                  -- matrix: 3(omega,alpha,beta) or 4(+gamma) x max lags in model, e.g. Garch(2,2) => 2 Columns
    # garch$optimcontrol          -- list: wrapper for the 'control' parameter sent to optim()

## --- After estimation:

    # garch$Estimated$value        -- scalar: log-likelihood value
    # garch$Estimated$error        -- logical: T/F indication if an error occurred during estimation
    # garch$Estimated$pars 
    # garch$Estimated$hessian 
    # garch$Estimated$se 
    # garch$Estimated$optimoutput  -- list: the full returned value from optim().  Only set if the verbose param = True


garchtype = list(noGarch=0,general=1,gjr=2)

## --- GARCH_CLASS Definition --- ####
garch <- setClass(Class = "garch_class",
               slots = c(h="numeric",nr.pars="integer",order="numeric"), 
               contains = c("namedList")
               )

## Initialise with no params
setMethod("initialize","garch_class",
          function(.Object,...){
            .Object <- callNextMethod(.Object,...)
            .Object@nr.pars <- as.integer(0)
            .Object@order <- as.integer(0)
            # Return:
            .Object
          })

setGeneric(name="garch",
           function(type="numeric",order="numeric")standardGeneric("garch"),
           valueClass = "garch_class",
           signature = c("type","order")
          )

## Constructor with 1 param - creates a Garch(1,1)
setMethod("garch",signature = c("numeric","missing"),
          function(type){
            # Create a GARCH(1,1) model
            this <- garch(type,c(1,1))
            return(this)
          })

## Constructor with 2 params
setMethod("garch",signature = c("numeric","numeric"),
          function(type,order){
           
            this <- new("garch_class")
            this$type <- type
            if(type == garchtype$noGarch) return(this)
            
            this@order <- order
            this <- .setInitPars(this)
            this$optimcontrol <- list(fnscale = -1, maxit = 1500, reltol = 1e-7)

            return(this)
          })


## --- Public Methods --- ####


## -- estimateGARCH() ####
setGeneric(name=".estimateGARCH",
           function(e="numeric",garch="garch_class",estimationControl="list")standardGeneric(".estimateGARCH"),
           valueClass = "garch_class",
           signature = c("e","garch","estimationControl")
)
setMethod(".estimateGARCH",signature = c("numeric","garch_class","list"),
          function(e,garch,estimationControl){
            this <- garch
            
            if(this$type == garchtype$noGarch) {
              message("Cannot estimateGARCH for type: NoGarch")
              return(this)
            }
              
            if(!is.null(estimationControl$calcSE)) calcSE <- estimationControl$calcSE else calcSE <- FALSE
            if(!is.null(estimationControl$verbose)) verbose <- estimationControl$verbose else verbose <- FALSE

            # Start the Estimation process:
            if (verbose) {
              this$optimcontrol$trace <- 10 
              cat("\nEstimating GARCH object...\n")
            } else this$optimcontrol$trace <- 0
            
            
            # Get Optimpars from garch$pars
            optimpars <- as.vector(this$pars)
            names(optimpars) <- rownames(this$pars)
            
            # Now call optim:
            tmp <- NULL
            try(tmp <- optim(optimpars,.loglik.garch.univar,gr=NULL,e,this,method="BFGS",control=this$optimcontrol,hessian=calcSE))      
          
            ## --- Attach results of estimation to the object --- ##
            this$Estimated <- list()
            
            # An unhandled error could result in a NULL being returned by optim()
            if (is.null(tmp)) {
              this$Estimated$value <- -1e10
              this$Estimated$error <- TRUE
              warning("estimateGARCH() - optim failed and returned NULL. Check the optim controls & starting params")
              return(this)
            }
            if (tmp$convergence!=0) { 
              this$Estimated$value <- -1e10
              this$Estimated$error <- TRUE
              this$Estimated$optimoutput <- tmp
              warning("estimateGARCH() - failed to converge. Check the optim controls & starting params")
              return(this)
            }
            
            this$Estimated$value <- tmp$value
            this$Estimated$error <- FALSE
            
            #Update the GARCH object paramters using optimised pars:
            this$Estimated$pars <- .parsVecToMatrix(this,tmp$par)
            # Get conditional variance
            this <- .calculate_h(this,e)
            
            # Calc Std Errors
            if (calcSE) {
              this$Estimated$hessian <- round(tmp$hessian,5)
              StdErrors <- NULL
              try(StdErrors <- sqrt(-diag(qr.solve(tmp$hessian))))
              if(is.null(StdErrors)) {
                this$Estimated$se <- matrix(NA,nrow=this@nr.pars)
              }else {
                this$Estimated$se <- matrix(StdErrors,nrow=this@nr.pars)
              }
              rownames(this$Estimated$se) <- rownames(this$pars)
              colnames(this$Estimated$se) <- "se"
            }
            if (verbose) this$Estimated$optimoutput <- tmp

            return(this)
          })

setGeneric(name="estimateGARCH",
           function(e="numeric",garch="garch_class",estimationControl="list")standardGeneric("estimateGARCH"),
           valueClass = "garch_class",
           signature = c("e","garch","estimationControl")
)
setMethod("estimateGARCH",signature = c("numeric","garch_class","list"),
          function(e,garch,estimationControl){
            return(.estimateGARCH(e,garch,estimationControl))
          })
setMethod("estimateGARCH",signature = c("numeric","garch_class","missing"),
          function(e,garch){
            estimationControl <- list()
            estimationControl$calcSE <- FALSE
            estimationControl$verbose <- FALSE
            return(.estimateGARCH(e,garch,estimationControl))
          })



## --- PRIVATE METHODS --- ####

## -- .setInitPars() -- ####

setGeneric(name=".setInitPars",
           function(garch="garch_class")standardGeneric(".setInitPars"),
           valueClass = "garch_class",
           signature = c("garch")
)
setMethod(".setInitPars",signature = c("garch_class"),
          function(garch){
            this <- garch
            
            if(this$type == garchtype$noGarch) {
              message("Cannot create Garch$pars for type: NoGarch")
              return(this)
            }
              
            maxLag <- max(this@order)
            
            # Set the row names:
            GarchparsRownames <- c("omega","alpha","beta","gamma")
            
            if(this$type == garchtype$general) {
              this@nr.pars <- as.integer(3)
              pars <- matrix(nrow = this@nr.pars,ncol = maxLag)
              rownames(pars) <- GarchparsRownames[1:this@nr.pars]
              for(n in 1:maxLag){
                pars["omega",n] <- 0.05
                pars["alpha",n] <- 0.05
                pars["beta",n] <- 0.85
              }
              this$pars <- pars
            }
            if(this$type == garchtype$gjr) {
              this@nr.pars <- as.integer(4)
              pars <- matrix(nrow = this@nr.pars,ncol = maxLag)
              rownames(pars) <- GarchparsRownames[1:this@nr.pars]
              for(n in 1:maxLag){
                pars["omega",n] <- 0.05
                pars["alpha",n] <- 0.05
                pars["beta",n] <- 0.85
                pars["gamma",n] <- 0.05
              }
              this$pars <- pars
            }
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
          })

## -- .parsVecToMatrix() ####
setGeneric(name=".parsVecToMatrix",
           function(garch="garch_class",pars="numeric")standardGeneric(".parsVecToMatrix"),
           valueClass = "matrix",
           signature = c("garch","pars")
)
setMethod(".parsVecToMatrix",signature = c("garch_class","numeric"),
          function(garch,pars){
            this <- garch
            
            if(this$type == garchtype$noGarch) {
              message("Cannot create Garch Params for type: NoGarch")  
              return(this)
            }
            
            maxLag <- max(this@order)
            
            # Set the row names:
            GarchparsRownames <- c("omega","alpha","beta","gamma")
            # Return the formatted matrix
            matrix(pars,nrow = this@nr.pars ,ncol = maxLag,dimnames = list(GarchparsRownames[1:this@nr.pars],"Est"))
            
          })


## -- .calculate_h() ####
setGeneric(name=".calculate_h",
           function(garchobj="garch_class",e="numeric")standardGeneric(".calculate_h"),
           valueClass = "garch_class",
           signature = c("garchobj","e")
)
setMethod(".calculate_h",signature = c("garch_class","numeric"),
          function(garchobj,e){
            
            this <- garchobj
            
            if(this$type == garchtype$noGarch) {
              message("Cannot calculate h(t) for type: NoGarch")
              return(this)
            }

            Tobs <- NROW(e)

            h <- rep(0,Tobs)
            h[1] <- sum(e*e)/Tobs
            # TODO: Extend the below to handle more lags
            for(t in 2:Tobs) {
              h[t] <- this$Estimated$pars["omega",1] + this$Estimated$pars["alpha",1]*(e[t-1])^2 + this$Estimated$pars["beta",1]*h[t-1] 
              if(this$type == garchtype$gjr) h[t] <- h[t] + this$Estimated$pars["gamma",1]*(min(e[t-1],0))^2
            }
            
            this@h <- h
            return(this)
          })

setGeneric(name="calculate_h",
           function(garchobj="garch_class",e="numeric")standardGeneric("calculate_h"),
           valueClass = "garch_class",
           signature = c("garchobj","e")
)

setMethod("calculate_h",signature = c("garch_class","numeric"),
          function(garchobj,e){
           .calculate_h(garchobj,e) 
          })

## -- .loglik.garch.univar() ####
setGeneric(name=".loglik.garch.univar",
           function(optimpars="numeric",e="numeric",garchobj="garch_class")standardGeneric(".loglik.garch.univar"),
           valueClass = "numeric",
           signature = c("optimpars","e","garchobj")
)

setMethod(".loglik.garch.univar",signature = c("numeric","numeric","garch_class"),
          function(optimpars,e,garchobj){
          
            error <- -1e10
            this <- garchobj

            ## ======== constraint checks ======== ##
            # Check if any parameter is negative:
            if(min(optimpars,na.rm = TRUE) < 0) return(error)
            #
            if (optimpars["omega"] <= 0) return(error)
            if (optimpars["alpha"]+optimpars["beta"] >= 1) return(error)

            ## ======== calculate loglikelihood ======== ##
            this$Estimated$pars <- .parsVecToMatrix(this,optimpars)
            this <- .calculate_h(this,e)

            #Return the LogLiklihood value:
            sum( -0.5*log(2*pi) - 0.5*log(this@h) - 0.5*(e*e)/this@h ) 
            
          })

## --- Override Methods --- ####

## -- plot() ####
setMethod("plot",signature = c(x="garch_class",y="missing"),
          function(x, y, ...){
            plot.default(x=x@h, type='l', ylab = "Cond.Variance", ...)
          })


## -- summary() ####
setMethod("summary",signature="garch_class",
          function(object,...){
            this <- object
            
            TypeNames <- c("No Garch","General","GJR Garch")
            
            if(this$type == garchtype$noGarch){
              cat("\nGARCH OBJECT\n")
              cat("\nType:",TypeNames[this$type+1])
              cat("\nCannot be estimated - this type only exists to support `mtvgarch` objects")
              return()
            } 
            
            if(is.null(this$Estimated)){
              cat("\nGARCH OBJECT\n")
              cat("\nType:",TypeNames[this$type+1])
              cat("\nThis Garch object has not been estimated yet.")
              return()
            } else {

              parsVec <-  round(as.vector(this$Estimated$pars),6)
              parsRows <- NROW(this$Estimated$pars)

              if(!is.null(this$Estimated$se) ){
                seVec <- round(as.vector(this$Estimated$se),6)
                seVecSig <- vector("character", length(seVec))

                for(n in seq_along(parsVec)){
                  if(is.nan(seVec[n])) {
                    seVecSig[n] <- "   "
                  } else {
                    # Calculate a significance indicator
                    if(seVec[n]*2 < abs((parsVec[n]/100)) ) { (seVecSig[n] <- "***") } 
                    else if(seVec[n]*2 < abs((parsVec[n]/10)) ) { (seVecSig[n] <- "** ") } 
                    else if(seVec[n]*2 < abs((parsVec[n])) ) { (seVecSig[n] <- "*  ") } 
                  }
                } 
              } else {
                seVec <- rep(NaN,length(this$pars))
                seVecSig <- rep("   ", length(seVec))
              }
                
              seMat <- matrix(seVec,nrow=parsRows)
              colnames(seMat) <- paste("se" ,1:max(this@order),sep = "")
              # Build Results table and insert the significance indicators
              results <- data.frame(NA,stringsAsFactors = FALSE)
              for (n in 1:NCOL(this$Estimated$pars)){
                sig <- matrix(seVecSig[1:parsRows],nrow=parsRows)
                results <- cbind(results,round(this$Estimated$pars[,n,drop=F],6),seMat[,n,drop=F],sig)
                seVecSig <- tail(seVecSig,-parsRows)
              }
            }
            
            cat("\nGARCH OBJECT\n")
            cat("\nOrder: (",this@order[1],",",this@order[2],")")
            cat("\nEstimation Results:\n")
            print(results[,-1])
            cat("\nLog-liklihood value: ",this$Estimated$value)

          })


