## -- The MTVGARCH package supports a number of Correlation objects
## -- This class file maintains the structure for STCC1 (STCC with One Transition) & STCC2 (STCC with Two Transitions)
## -- Both classes only support single transitions at this time


## --- stcc1_class Definition --- ####
stcc1 <- setClass(Class = "stcc1_class",
               slots = c(ntvg="ntvgarch_class",nr.corPars="integer",nr.trPars="integer",z="matrix"),
               contains = c("namedList")
               )

## --- Initialise --- ####
setMethod("initialize","stcc1_class",
          function(.Object,...){
            .Object <- callNextMethod(.Object,...)

            # Default initial values
            .Object$N <- 0
            .Object$e <- matrix("numeric")
            .Object$Tobs <- 0

            .Object$shape <- corrshape$single
            .Object$speedopt <- corrspeedopt$eta
            .Object$optimcontrol <- list(fnscale = -1, reltol = 1e-5, trace = 0)

            .Object$P1 <- matrix("numeric")
            .Object$P2 <- matrix("numeric")
            .Object$pars <- c(3,0.5)
            names(.Object$pars) <- c("speed","loc")


            # Return:
            .Object
          })

## -- Constructor:stcc1 -- ####
setGeneric(name="stcc1",
           valueClass = "stcc1_class",
           signature = c("ntvgarchObj"),
           def = function(ntvgarchObj){
             this <- new("stcc1_class")

             ## -- Do validation checks -- ####

             if(class(ntvgarchObj)[1] != "ntvgarch_class"){
               warning("a valid instance of the ntvgarch_class is required to create an STCC1 model")
               return(this)
             }
             # End validation

             # Set Default Values:
             N <- this$N <- ntvgarchObj$N
             this$Tobs <- ntvgarchObj$Tobs
             this$e <- ntvgarchObj$e
             z <- this@z <- ntvgarchObj$z

             # Extract the Data & Estimated components from the ntvgarch
             this@ntvg <- ntvgarchObj

             this$st <- (1:ntvgarchObj$Tobs)/ntvgarchObj$Tobs
             this@nr.corPars <- as.integer((N^2-N)/2)
             this@nr.trPars <- as.integer(2)

             # # Extract the Data & Estimated components from the ntvgarch
             # this@ntvg <- list()
             # for(n in 1:ntvgarchObj@N){
             #   this@ntvg[[n]] <- list()
             #   this@ntvg[[n]]$tv <- ntvgarchObj[[n]]$Estimated$tv
             #   this@ntvg[[n]]$garch <- ntvgarchObj[[n]]$Estimated$garch
             #   this$e[,n] <- ntvgarchObj[[n]]@e
             # }
             # names(this@ntvg) <- names(ntvgarchObj)


             # Use the unconditional correlation of the first & last 1/3 of the data as starting values
             zDiv3 <- round(this$Tobs/3)
             z.start <- 1
             z.end <- zDiv3
             this$P1 <- cor(z[(z.start:z.end),])

             z.start <- this$Tobs - zDiv3
             z.end <- this$Tobs
             this$P2 <- cor(z[(z.start:z.end),])

             return(this)
           }
)


## --- stcc2_class Definition --- ####
stcc2 <- setClass(Class = "stcc2_class",
                  slots = c(ntvg="ntvgarch_class",nr.corPars="integer",nr.trPars="integer",z="matrix"),
                  contains = c("namedList")
)

## --- Initialise --- ####
setMethod("initialize","stcc2_class",
          function(.Object,...){
            .Object <- callNextMethod(.Object,...)

            # Default initial values
            .Object$N <- 0
            .Object$e <- matrix("numeric")
            .Object$Tobs <- 0

            .Object$shape <- corrshape$single
            .Object$speedopt <- corrspeedopt$eta
            .Object$optimcontrol <- list(fnscale = -1, reltol = 1e-5, trace = 0)

            .Object$P1 <- matrix("numeric")
            .Object$P2 <- matrix("numeric")
            .Object$P3 <- matrix("numeric")
            .Object$pars <- c(3,0.33,3,0.66)
            names(.Object$pars) <- c("speed1","loc1","speed2","loc2")

            # TODO: Update parscale, ndeps
            .Object$optimcontrol$parscale

            # Return:
            .Object
          })

## -- Constructor:stcc2 -- ####
setGeneric(name="stcc2",
           valueClass = "stcc2_class",
           signature = c("ntvgarchObj"),
           def = function(ntvgarchObj){
             this <- new("stcc2_class")

             ## -- Do validation checks -- ####

             if(class(ntvgarchObj)[1] != "ntvgarch_class"){
               warning("a valid instance of the ntvgarch_class is required to create an STCC1 model")
               return(this)
             }
             # End validation

             # Set Default Values:
             N <- this$N <- ntvgarchObj$N
             this$Tobs <- ntvgarchObj$Tobs
             this$e <- ntvgarchObj$e
             z <- this@z <- ntvgarchObj$z

             # Extract the Data & Estimated components from the ntvgarch
             this@ntvg <- ntvgarchObj

             this$st <- (1:ntvgarchObj$Tobs)/ntvgarchObj$Tobs
             this@nr.corPars <- as.integer((N^2-N)/2)
             this@nr.trPars <- as.integer(2)

             # Use the correlation of the first & middle & last 1/3 of the data as starting values
             zDiv3 <- round(this$Tobs/3)
             z.start <- 1
             z.end <- zDiv3
             this$P1 <- cor(z[(z.start:z.end),])

             z.start <- z.end+1
             z.end <- this$Tobs - zDiv3 - 1
             this$P2 <- cor(z[(z.start:z.end),])

             z.start <- this$Tobs - zDiv3
             z.end <- this$Tobs
             this$P3 <- cor(z[(z.start:z.end),])

             return(this)
           }
)

## -- .calc.Gt -- ####
setGeneric(name=".calc.Gt",
           valueClass = "matrix",
           signature = c("stcc1Obj"),
           def = function(stcc1Obj){
             this <- stcc1Obj

             speed <- this$Estimated$pars["speed"]
             loc1 <- this$Estimated$pars["loc1"]

             st_c <- 0
             if(this$shape == corrshape$single) { st_c <- this$st - loc1 }

             G <- 0
             if(this$speedopt == corrspeedopt$gamma) { G <- 1/(1+exp(-speed*st_c)) }
             if(this$speedopt == corrspeedopt$gamma_std) { G <- 1/(1+exp(-speed/sd(this$st)*st_c)) }
             if(this$speedopt == corrspeedopt$eta) { G <- 1/(1+exp(-exp(speed)*st_c)) }

             return(matrix(G,nrow = this$Tobs,ncol = 1))
           }
)

## -- .calc.Pt -- ####
setGeneric(name=".calc.Pt",
           valueClass = "matrix",
           signature = c("stcc1Obj"),
           def =   function(stcc1Obj){
             this <- stcc1Obj

             vP1 <- vecL(this$Estimated$P1)
             vP2 <- vecL(this$Estimated$P2)

             Gt <- .calc.Gt(this)
             Pt <- apply(Gt,MARGIN = 1,FUN = function(X,P1,P2) ((1-X)*P1 + X*P2), P1=vP1, P2=vP2)

             if(is.vector(Pt)) Pt <- matrix(Pt,ncol = 1) else Pt <- t(Pt)

             return(Pt)

           }
)

## -- .calc.Gt2 -- ####
setGeneric(name=".calc.Gt2",
           valueClass = "matrix",
           signature = c("stcc2Obj"),
           def = function(stcc2Obj){
             this <- stcc2Obj

             speed1 <- this$Estimated$pars["speed1"]
             loc1 <- this$Estimated$pars["loc1"]
             speed2 <- this$Estimated$pars["speed2"]
             loc2 <- this$Estimated$pars["loc2"]

             st_c_1 <- this$st - loc1
             st_c_2 <- this$st - loc2

             G <- matrix(0,nrow = this$Tobs, ncol = 2)
             if(this$speedopt == corrspeedopt$gamma) { G[,1] <- 1/(1+exp(-speed1*st_c_1)) }
             if(this$speedopt == corrspeedopt$gamma_std) { G[,1] <- 1/(1+exp(-speed1/sd(this$st)*st_c_1)) }
             if(this$speedopt == corrspeedopt$eta) { G[,1] <- 1/(1+exp(-exp(speed1)*st_c_1)) }
             #
             if(this$speedopt == corrspeedopt$gamma) { G[,2] <- 1/(1+exp(-speed2*st_c_2)) }
             if(this$speedopt == corrspeedopt$gamma_std) { G[,2] <- 1/(1+exp(-speed2/sd(this$st)*st_c_2)) }
             if(this$speedopt == corrspeedopt$eta) { G[,2] <- 1/(1+exp(-exp(speed2)*st_c_2)) }

             return(matrix(G,nrow = this$Tobs,ncol = 2))
           }
)



## -- .calc.Pt2 -- ####
setGeneric(name=".calc.Pt2",
           valueClass = "matrix",
           signature = c("stcc2Obj"),
           def =   function(stcc2Obj){
             this <- stcc2Obj

             vP1 <- matrix(vecL(this$Estimated$P1),nrow=1)
             vP2 <- matrix(vecL(this$Estimated$P2),nrow=1)
             vP3 <- matrix(vecL(this$Estimated$P3),nrow=1)

             Gt <- .calc.Gt2(this) # T x 2
             Pt <- ((1-Gt[,2])*(1-Gt[,1]))%*%vP1 + ((1-Gt[,2])*Gt[,1])%*%vP2 + Gt[,2]%*%vP3 # T x N(N-1)/2
             if(is.vector(Pt)) Pt <- matrix(Pt,ncol = 1)

             return(Pt)
           }
)


## -- .loglik.stcc1() --####
setGeneric(name=".loglik.stcc1",
           valueClass = "numeric",
           signature = c("optimpars","z","stcc1Obj"),
           def = function(optimpars,z,stcc1Obj){

             err_output <- -1e10
             this <- stcc1Obj

             this$Estimated$pars <- tail(optimpars,this@nr.trPars)
             tmp.par <- optimpars

             vP1 <- tmp.par[1:this@nr.corPars]
             mP <- unVecL(vP1)
             eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
             # Check for SPD - positive-definite check:
             if (min(eig$values) <= 0) return(err_output)
             this$Estimated$P1 <- mP

             #Remove the P1 corPars, then extract the P2 corPars
             tmp.par <- tail(tmp.par,-this@nr.corPars)
             vP2 <- tmp.par[1:this@nr.corPars]
             mP <- unVecL(vP2)
             eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
             # Check for SPD - positive-definite check:
             if (min(eig$values) <= 0) return(err_output)
             this$Estimated$P2 <- mP

             #### ======== constraint checks ======== ####

             # Check 2: Check the boundary values for speed params:
             speed <- this$Estimated$pars[1]
             maxSpeed <- switch(this$speedopt,1000,(1000/sd(this$st)),7.0,0.30)
             if (speed > maxSpeed) return(err_output)
             if (speed < 0) return(err_output)

             # Check 3: Check the locations fall within min-max values of st
             loc <- this$Estimated$pars[2]
             if (loc < min(this$st)) return(err_output)
             if (loc > max(this$st)) return(err_output)


             #### ======== calculate loglikelihood ======== ####
             Pt <- .calc.Pt(this)  # T x nr.corPars

             llt <- vector("numeric")
             for(t in 1:this$Tobs) {
               mPt <- unVecL(Pt[t,,drop=FALSE])
               llt[t] <- -0.5*log(det(mPt)) -0.5*( z[t,,drop=FALSE] %*% (solve(mPt)) %*% t(z[t,,drop=FALSE]) )
             }
             # Return:
             return(sum(llt))

           }
)

## --- estimateSTCC1 --- ####
setGeneric(name="estimateSTCC1",
           valueClass = "stcc1_class",
           signature = c("stcc1Obj","estimationCtrl"),
           def = function(stcc1Obj,estimationCtrl){
             this <- stcc1Obj
             e <- this$e

             calcSE <- estimationCtrl$calcSE
             verbose <- estimationCtrl$verbose

             this$Estimated <- list()

             optimpars <- c( vecL(this$P1), vecL(this$P2), this$pars )
             optimpars <- optimpars[!is.na(optimpars)]

             ### ---  Call optim to calculate the estimate --- ###
             if (verbose) this$optimcontrol$trace <- 10

             # Filter the data:
             z <- this$z
             tmp <- NULL
             try(tmp <- optim(optimpars,.loglik.stcc1,z,this,gr=NULL,method="BFGS",control=this$optimcontrol))

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
               this$Estimated$P1 <- unVecL(tmp.par[1:this@nr.corPars])
               tmp.par <- tail(tmp.par,-this@nr.corPars)
               this$Estimated$P2 <- unVecL(tmp.par[1:this@nr.corPars])
               tmp.par <- tail(tmp.par,-this@nr.corPars)
               this$Estimated$pars <- tmp.par
               names(this$Estimated$pars) <- names(this$pars)

               if (calcSE) {
                 cat("\nCalculating STCC standard errors...\n")
                 this$Estimated$hessian <- NULL
                 try(this$Estimated$hessian <- optimHess(tmp$par,.loglik.stcc1,z,this,gr=NULL,control=this$optimcontrol))
                 # Handle optimHess returns non-matrix
                 vecSE <- vector("numeric")
                 try(vecSE <- sqrt(-diag(invertHess(this$Estimated$hessian))))
                 if(length(vecSE) > 0) {
                   this$Estimated$P1.se <- unVecL(vecSE[1:this@nr.corPars])
                   vecSE <- tail(vecSE,-this@nr.corPars)
                   this$Estimated$P2.se <- unVecL(vecSE[1:this@nr.corPars])
                   vecSE <- tail(vecSE,-this@nr.corPars)

                   this$Estimated$pars.se <- vecSE
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
setMethod("estimateSTCC1",signature = c("stcc1_class","missing"),
          function(stcc1Obj){
            estimationControl <- list(calcSE = TRUE,verbose = TRUE)
            estimateSTCC1(stcc1Obj,estimationControl)
          })



## -- .loglik.stcc2() --####
setGeneric(name=".loglik.stcc2",
           valueClass = "numeric",
           signature = c("optimpars","z","stcc2Obj"),
           def = function(optimpars,z,stcc2Obj){

             err_output <- -1e10
             this <- stcc2Obj

             this$Estimated$pars <- tail(optimpars,2*this@nr.trPars)
             tmp.par <- optimpars

             # P1 always full of parameters:
             vP1 <- tmp.par[1:this@nr.corPars]
             mP <- unVecL(vP1)
             eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
             # Check for SPD - positive-definite check:
             if (min(eig$values) <= 0) return(err_output)
             this$Estimated$P1 <- mP

             #Remove the P1 corPars, then extract the P2 corPars
             tmp.par <- tail(tmp.par,-this@nr.corPars)
             nr.freepars2 <- sum(this$sel.vec2,na.rm=TRUE) # count of number of free pars
             freepars <- tmp.par[1:nr.freepars2]
             vP2<-vP1 # copy of P1
             # replace vP2 elements that have 1's in selvec by free pars
             vP2[which(as.logical(this$sel.vec2))] <- freepars
             mP <- unVecL(vP2)
             eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
             # Check for SPD - positive-definite check:
             if (min(eig$values) <= 0) return(err_output)
             this$Estimated$P2 <- mP

             #Remove the P2 corPars, then extract the P3 corPars
             tmp.par <- tail(tmp.par,-nr.freepars2)
             nr.freepars3 <- sum(this$sel.vec3,na.rm=TRUE) # count of number of free pars
             freepars <- tmp.par[1:nr.freepars3]
             vP3<-vP2 # copy of P2
             # replace vP3 elements that have 1's in selvec by free pars
             vP3[which(as.logical(this$sel.vec3))] <- freepars
             mP <- unVecL(vP3)
             eig <- eigen(mP,symmetric=TRUE,only.values=TRUE)
             # Check for SPD - positive-definite check:
             if (min(eig$values) <= 0) return(err_output)
             this$Estimated$P3 <- mP


             #### ======== constraint checks ======== ####

             # Check 2.1: Check the boundary values for speed params:
             pos = 1
             speed <- this$Estimated$pars[pos]
             maxSpeed <- switch(this$speedopt,1000,(1000/sd(this$st)),7.0,0.30)
             if (speed > maxSpeed) return(err_output)
             if (speed < 0) return(err_output)
             pos = pos+1

             # Check 3.1: Check the locations fall within min-max values of st
             loc1 <- this$Estimated$pars[pos]
             if (loc1 < min(this$st)) return(err_output)
             if (loc1 > max(this$st)) return(err_output)
             pos = pos+1

             # Check 2.2: Check the boundary values for speed params:
             speed <- this$Estimated$pars[pos]
             maxSpeed <- switch(this$speedopt,1000,(1000/sd(this$st)),7.0,0.30)
             if (speed > maxSpeed) return(err_output)
             if (speed < 0) return(err_output)
             pos = pos+1

             # Check 3.2: Check the locations fall within min-max values of st
             loc2 <- this$Estimated$pars[pos]
             if (loc2 < min(this$st)) return(err_output)
             if (loc2 > max(this$st)) return(err_output)

             if (loc2 < loc1) return(err_output)
             #if (loc2 - loc1 < 0.2) return(err_output)


             #### ======== calculate loglikelihood ======== ####
             Pt <- .calc.Pt2(this)  # T x nr.corPars

             llt <- vector("numeric")
             for(t in 1:this$Tobs) {
               mPt <- unVecL(Pt[t,,drop=FALSE])
               llt[t] <- -0.5*log(det(mPt)) -0.5*( z[t,,drop=FALSE] %*% (solve(mPt)) %*% t(z[t,,drop=FALSE]) )
             }
             # Return:
             return(sum(llt))

           }
)


## --- estimateSTCC2 --- ####
setGeneric(name="estimateSTCC2",
           valueClass = "stcc2_class",
           signature = c("stcc2Obj","estimationCtrl"),
           def = function(stcc2Obj,estimationCtrl){

             this <- stcc2Obj
             this$Estimated <- list()
             e <- this$e
             z <- this@z

             calcSE <- estimationCtrl$calcSE
             verbose <- estimationCtrl$verbose

             veclP1 <- vecL(this$P1)
             veclP2 <- vecL(this$P2)
             veclP3 <- vecL(this$P3)

             # Set the selection vectors to manage free/restricted parameters
             this$sel.vec2 = as.numeric(!is.na(veclP2))
             this$sel.vec3 = as.numeric(!is.na(veclP3))

             optimpars <- c( veclP1, veclP2, veclP3, this$pars )
             optimpars <- optimpars[!is.na(optimpars)]

             # Manage the optimControl params

             ### ---  Call optim to calculate the estimate --- ###
             if (verbose) this$optimcontrol$trace <- 10

             tmp <- NULL
             try(tmp <- optim(optimpars,.loglik.stcc2,z,this,gr=NULL,method="BFGS",control=this$optimcontrol,hessian=calcSE))

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
               vP1 <- tmp.par[1:this@nr.corPars]
               this$Estimated$P1 <- unVecL(vP1)
               tmp.par <- tail(tmp.par,-this@nr.corPars)

               nr.freepars2 <- sum(this$sel.vec2,na.rm=TRUE) # count of number of free pars
               freepars <- tmp.par[1:nr.freepars2]
               vP2 <- vP1 # copy of P1
               # replace vP2 elements that have 1's in selvec by free pars
               vP2[which(as.logical(this$sel.vec2))] <- freepars
               this$Estimated$P2 <- unVecL(vP2)
               tmp.par <- tail(tmp.par,-nr.freepars2)

               nr.freepars3 <- sum(this$sel.vec3,na.rm=TRUE) # count of number of free pars
               freepars <- tmp.par[1:nr.freepars3]
               vP3 <- vP2 # copy of P2
               # replace vP3 elements that have 1's in selvec by free pars
               vP3[which(as.logical(this$sel.vec3))] <- freepars
               this$Estimated$P3 <- unVecL(vP3)
               tmp.par <- tail(tmp.par,-nr.freepars3)

               # TO DO : probably going to fall over but likely never to be generlised to double
               trPars <- tail(tmp$par,(2*this@nr.trPars))
               pars1 <- trPars[1:this@nr.trPars]
               pars2 <- tail(trPars,-this@nr.trPars)
               #
               this$Estimated$pars <- c(pars1,pars2)
               names(this$Estimated$pars) <- names(this$pars)

               if (calcSE) cat("\nCalculating STCC standard errors...\n")
               if (calcSE) {
                 this$Estimated$hessian <- tmp$hessian
                 vecSE <- vector("numeric")
                 try(vecSE <- sqrt(-diag(qr.solve(tmp$hessian))))

                 if(length(vecSE) > 0) {
                   vecSE1 <- vecSE[1:this@nr.corPars]
                   this$Estimated$P1.se <- unVecL(vecSE1)
                   vecSE <- tail(vecSE,-this@nr.corPars)

                   vecSE2 <- vecSE1
                   vecSE2[which(as.logical(this$sel.vec2))] <- vecSE[1:nr.freepars2]
                   this$Estimated$P2.se <- unVecL(vecSE2)
                   vecSE <- tail(vecSE,-nr.freepars2)

                   vecSE3 <- vecSE2
                   vecSE3[which(as.logical(this$sel.vec3))] <- vecSE[1:nr.freepars3]
                   this$Estimated$P3.se <- unVecL(vecSE3)
                   vecSE <- tail(vecSE,-nr.freepars3)

                   vecSE1 <- vecSE[1:this@nr.trPars]
                   vecSE2 <- tail(vecSE, -this@nr.trPars)
                   this$Estimated$pars.se <- c(vecSE1,vecSE2)
                   #
                   names(this$Estimated$pars.se) <- names(this$pars)
                 }
               }
               this$Estimated$Pt <- .calc.Pt2(this)

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



##   Might be needed later:
## --- unCorrelateData --- ####
setGeneric(name="unCorrelateData",
           valueClass = "matrix",
           signature = c("e","stcc1Obj"),
           def = function(e,stcc1Obj){
             this <- stcc1Obj

             # Return matrix 'u' of uncorrelated data:
             u <- matrix(0,nrow = NROW(e),ncol = NCOL(e))

             Pt <- this$Estimated$Pt
             for(t in 1:this$Tobs){
               P <- unVecL(Pt[t,,drop=FALSE])
               P.sqrt <- sqrt_mat1(P)
               u[t,] <- t(solve(P.sqrt) %*% t(e[t,,drop=FALSE]))
             }

             return(u)
           }
)

## -- .dG_dtr(STCC) ####
setGeneric(name=".dG_dtr",
           valueClass = "matrix",
           signature = c("stccObj","trNum"),
           def =  function(stccObj,trNum){
             this <- stccObj

             rtn <- matrix(nrow=this$Tobs,ncol=this@nr.trPars)

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
               rtn[,col_idx] <- G * (1-G) * (this$st - loc1)
               col_idx <- 2
               rtn[,col_idx] <- -1 * G * (1-G) * speed
             }

             if(this$speedopt == corrspeedopt$gamma_std) {
               col_idx <- 1
               rtn[,col_idx] <- G * (1-G) * (this$st - loc1) / sd(this$st)
               col_idx <- 2
               rtn[,col_idx] <- -1 * G * (1-G) * speed / sd(this$st)
             }

             if(this$speedopt == corrspeedopt$eta) {
               col_idx <- 1
               rtn[,col_idx] <- G * (1-G) * (this$st - loc1) * exp(speed)
               col_idx <- 2
               rtn[,col_idx] <- -1 * G * (1-G) * exp(speed)
             }

             # # TODO:  Complete for other corrshapes & other corr speedOpt's:

             return(rtn)

           }
)

