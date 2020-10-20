## --- Correlation Structures --- ##
##
## -- We have a number of Correlation objects, with many common properties
## -- This class file maintains the Structure for STCC1 (STCC with One Transition)


## --- stcc1_class Definition --- ####
stcc1 <- setClass(Class = "stcc1_class",
               slots = c(st="numeric",nr.covPars="integer",nr.trPars="integer",Tobs="integer",N="integer"),
               contains = c("namedList")
               )

## --- Initialise --- ####
setMethod("initialize","stcc1_class",
          function(.Object,...){
            .Object <- callNextMethod(.Object,...)

            # Default initial values
            .Object$mtvgarch <- list()
            .Object@N <- as.integer(0)
            .Object@Tobs <- as.integer(0)
            .Object$shape <- corrshape$single
            .Object$speedopt <- corrspeedopt$eta
            .Object$optimcontrol <- list(fnscale = -1, maxit = 1000, reltol = 1e-5)

            # Return:
            .Object
          })

## -- Constructor:stcc1 -- ####
setGeneric(name="stcc1",
           valueClass = "stcc1_class",
           signature = c("mtvgarchObj"),
           def = function(mtvgarchObj){
             this <- new("stcc1_class")

             ## -- Do validation checks -- ####
             objType <- class(mtvgarchObj)
             if(objType[1] != "mtvgarch_class"){
               warning("a valid instance of the mtvgarch_class is required to create an STCC1 model")
               return(this)
             }
             # End validation

             # Add the Estimated components from the mtvgarch
             for(n in 1:mtvgarchObj@N){
               this$mtvgarch[[n]] <- list()
               this$mtvgarch[[n]]$tv <- mtvgarchObj[[n]]$Estimated$tv
               this$mtvgarch[[n]]$garch <- mtvgarchObj[[n]]$Estimated$garch
             }
             names(this$mtvgarch) <- names(mtvgarchObj)

             # Set Default Values:
             this@N <- mtvgarchObj@N
             this@st <- (1:mtvgarchObj@Tobs)/mtvgarchObj@Tobs
             this@Tobs <- mtvgarchObj@Tobs

             if(this$shape==corrshape$double) {
               this@nr.trPars <- as.integer(3)
               this$pars <- c(2.5,0.33,0.66)
             }else {
               this@nr.trPars <- as.integer(2)
               this$pars <- c(2.5,0.5,NA)
             }
             names(this$pars) <- c("speed","loc1","loc2")

             N <- this@N
             this@nr.covPars <- as.integer((N^2-N)/2)
             this$P1 <- matrix(0.3,N,N)
             diag(this$P1) <- 1
             this$P2 <- matrix(0.7,N,N)
             diag(this$P2) <- 1

             return(this)
           }
)

## -- calc.Gt -- ####
setGeneric(name="calc.Gt",
           valueClass = "matrix",
           signature = c("stcc1Obj"),
           def = function(stcc1Obj){
             this <- stcc1Obj

             if(is.null(this$Estimated$pars)){
               speed <- this$pars["speed"]
               loc1 <- this$pars["loc1"]
               loc2 <- this$pars["loc2"]
             } else {
               speed <- this$Estimated$pars["speed"]
               loc1 <- this$Estimated$pars["loc1"]
               loc2 <- this$Estimated$pars["loc2"]
             }

             st_c <- 0
             if(this$shape == corrshape$single) { st_c <- this@st - loc1 }
             if(this$shape == corrshape$double) { st_c <- (this@st - loc1)*(this@st - loc2) }
             if(this$shape == corrshape$double1loc) { st_c <- (this@st - loc1)^2 }

             G <- 0
             if(this$speedopt == corrspeedopt$gamma) { G <- 1/(1+exp(-speed*st_c)) }
             if(this$speedopt == corrspeedopt$gamma_std) { G <- 1/(1+exp(-speed/sd(this@st)*st_c)) }
             if(this$speedopt == corrspeedopt$eta) { G <- 1/(1+exp(-exp(speed)*st_c)) }

             return(matrix(G,nrow = this@Tobs,ncol = 1))
           }
)

## -- calc.Pt -- ####
setGeneric(name=".calc.Pt",
           valueClass = "matrix",
           signature = c("stcc1Obj"),
           def =   function(stcc1Obj){
             this <- stcc1Obj

             if(is.null(this$Estimated$P1)){
               vP1 <- .vecL(this$P1)
               vP2 <- .vecL(this$P2)
             } else {
               vP1 <- .vecL(this$Estimated$P1)
               vP2 <- .vecL(this$Estimated$P2)
             }

             Gt <- calc.Gt(this)
             Pt <- apply(Gt,MARGIN = 1,FUN = function(X,P1,P2) ((1-X)*P1 + X*P2), P1=vP1, P2=vP2)

             if(is.vector(Pt)) Pt <- matrix(Pt,ncol = 1)

             return(Pt)

           }
)


## -- .dG_dtr(STCC) ####
setGeneric(name=".dG_dtr",
           valueClass = "matrix",
           signature = c("stccObj","trNum"),
           def =  function(stccObj,trNum){

           this <- stccObj

           rtn <- matrix(nrow=this@Tobs,ncol=this@nr.trPars)

           G <- calc.Gt(this)  # T x 1

           loc1 <- loc2 <- speed <- 0
           if(trNum == 1){
             loc1 <- this$Estimated$pars["loc1"]
             loc2 <- this$Estimated$pars["loc2"]
             speed <- this$Estimated$pars["speed"]
           }

            # corrshape$single,
            if(this$speedopt == corrspeedopt$gamma) {
               col_idx <- 1
               rtn[,col_idx] <- G * (1-G) * (this@st - loc1)
               col_idx <- 2
               rtn[,col_idx] <- -1 * G * (1-G) * speed
            }

           if(this$speedopt == corrspeedopt$gamma_std) {
             col_idx <- 1
             rtn[,col_idx] <- G * (1-G) * (this@st - loc1) / sd(this@st)
             col_idx <- 2
             rtn[,col_idx] <- -1 * G * (1-G) * speed / sd(this@st)
           }

           if(this$speedopt == corrspeedopt$eta) {
             col_idx <- 1
             rtn[,col_idx] <- G * (1-G) * (this@st - loc1) * exp(speed)
             col_idx <- 2
             rtn[,col_idx] <- -1 * G * (1-G) * exp(speed)
           }

           # # TODO:  Complete for other corrshapes & other corr speedOpt's:

           return(rtn)

           }
)


## -- loglik.stcc1() --####
setGeneric(name=".loglik.stcc1",
           valueClass = "numeric",
           signature = c("optimpars","z","stcc1Obj"),
           def = function(optimpars,z,stcc1Obj){

             err_output <- -1e10
             this <- stcc1Obj

             #### ======== constraint checks ======== ####
             this$Estimated$pars <- tail(optimpars,this@nr.trPars)

             tmp.par <- optimpars

             vP1 <- tmp.par[1:this@nr.covPars]
             mP <- .unVecL(vP1)
             eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
             # Check for SPD - positive-definite check:
             if (min(eig$values) <= 0) return(err_output)

             #Remove the P1 covPars, then extract the P2 covPars
             tmp.par <- tail(tmp.par,-this@nr.covPars)
             vP2 <- tmp.par[1:this@nr.covPars]
             mP <- .unVecL(vP2)
             eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
             # Check for SPD - positive-definite check:
             if (min(eig$values) <= 0) return(err_output)


             # Check 2: Check the boundary values for speed params:
             speed <- this$Estimated$pars["speed"]
             maxSpeed <- switch(this$speedopt,1000,(1000/sd(this@st)),7.0,0.30)
             if (speed > maxSpeed) return(err_output)
             if (speed < 0) return(err_output)

             # Check 3: Check the locations fall within min-max values of st
             loc1 <- this$Estimated$pars["loc1"]
             if(this$shape == corrshape$double) loc2 <- this$Estimated$pars["loc2"] else loc2 <- NA
             if (loc1 < min(this@st)) return(err_output)
             if (loc1 > max(this@st)) return(err_output)
             if (!is.na(loc2)) {
               if (loc2 < min(this@st)) return(err_output)
               if (loc2 > max(this@st)) return(err_output)
             }


             #### ======== calculate loglikelihood ======== ####
             Pt <- .calc.Pt(this)  # T x nr.covPars

             llt <- vector("numeric")
             for(t in 1:this@Tobs) {
               mPt <- .unVecL(Pt[t,,drop=FALSE])
               llt[t] <- -0.5*log(det(mPt)) -0.5*( z[t,,drop=FALSE] %*% (qr.solve(mPt)) %*% t(z[t,,drop=FALSE]) )
             }
             # Return:
             return(sum(llt))

           }
)
## --- estimateSTCC1 --- ####
setGeneric(name="estimateSTCC1",
           valueClass = "stcc1_class",
           signature = c("z","stcc1Obj","estimationCtrl"),
           def = function(z,stcc1Obj,estimationCtrl){
             this <- stcc1Obj

             calcSE <- estimationCtrl$calcSE
             verbose <- estimationCtrl$verbose

             this$Estimated <- list()

             optimpars <- c( .vecL(this$P1), .vecL(this$P2), this$pars )
             optimpars <- optimpars[!is.na(optimpars)]

             ### ---  Call optim to calculate the estimate --- ###
             if (verbose) this$optimcontrol$trace <- 10

             tmp <- NULL
             try(tmp <- optim(optimpars,.loglik.stcc1,z,this,gr=NULL,method="BFGS",control=this$optimcontrol,hessian=calcSE))

             ### ---  Interpret the response from optim --- ###
             # An unhandled error could result in a NULL being returned by optim()
             if (is.null(tmp)) {
               this$Estimated$value <- -1e10
               this$Estimated$error <- TRUE
               return(this)
             }

             #Optim converged successfully => we expect tmp$par to have good estimates!
             if (tmp$convergence==0) {
               this$Estimated$error <- FALSE
               this$Estimated$value <- tmp$value

               tmp.par <- tmp$par
               this$Estimated$P1 <- .unVecL(tmp.par[1:this@nr.covPars])
               tmp.par <- tail(tmp.par,-this@nr.covPars)
               this$Estimated$P2 <- .unVecL(tmp.par[1:this@nr.covPars])
               tmp.par <- tail(tmp.par,-this@nr.covPars)
               this$Estimated$pars <- tmp.par
               if(this$shape != corrshape$double) this$Estimated$pars <- c(this$Estimated$pars,NA)
               names(this$Estimated$pars) <- names(this$pars)

               if (calcSE) {
                 this$Estimated$hessian <- tmp$hessian
                 vecSE <- vector("numeric")
                 try(vecSE <- sqrt(-diag(qr.solve(tmp$hessian))))

                 if(length(vecSE) > 0) {
                   this$Estimated$P1.se <- .unVecL(vecSE[1:this@nr.covPars])
                   vecSE <- tail(vecSE,-this@nr.covPars)
                   this$Estimated$P2.se <- .unVecL(vecSE[1:this@nr.covPars])
                   vecSE <- tail(vecSE,-this@nr.covPars)

                   this$Estimated$pars.se <- vecSE
                   if(this$shape != corrshape$double) this$Estimated$pars.se <- c(this$Estimated$pars.se,NA)
                   names(this$Estimated$pars.se) <- names(this$pars)
                 }
               }
               this$Estimated$Pt <- .calc.Pt(this)

             } else {
               #Failed to converge
               this$Estimated$error <- TRUE
               this$Estimated$value <- -1e10
               this$Estimated$optimoutput <- tmp
             }

             if (verbose) this$Estimated$optimoutput <- tmp
             #Return:
             return(this)
           }
)
## --- estimateSTCC1 --- ####
setMethod("estimateSTCC1",signature = c("matrix","stcc1_class","missing"),
          function(z,stcc1Obj){
            estimationControl <- list(calcSE = TRUE,verbose = TRUE)
            estimateSTCC1(z,stcc1Obj,estimationControl)
          })

## --- unCorrelateData --- ####
setGeneric(name="unCorrelateData",
           valueClass = "matrix",
           signature = c("e","stcc1Obj"),
           def = function(e,stcc1Obj){
             this <- stcc1Obj
             Pt <- this$Estimated$Pt
             # Return matrix 'u' of uncorrelated data:
             u <- matrix(0,nrow = NROW(e),ncol = NCOL(e))

             for(t in 1:this@Tobs){
               Pt_t <- .unVecL(Pt[t,])
               Pt_t_inv <- solve(Pt_t)
               u[t,] <- sqrt_mat1(Pt_t_inv) %*% e[t,]
             }

             return(u)
           }
)


## --- PRIVATE METHODS --- ####
