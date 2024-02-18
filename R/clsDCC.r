## -- The MTVGARCH package supports a number of Correlation objects
## -- This class file maintains the structure for DCC (Dynamic Conditional Correlation)


dcctype <- list(General=1,Dynamic=2)

## --- dcc_class Definition --- ####
dcc <- setClass(Class = "dcc_class",
               slots = c(st="numeric",nr.corPars="integer",nr.trPars="integer",Tobs="integer",N="integer",e="matrix"),
               contains = c("namedList")
               )

## --- Initialise --- ####
setMethod("initialize","dcc_class",
          function(.Object,...){
            .Object <- callNextMethod(.Object,...)

            # Default initial values
            .Object@N <- as.integer(0)
            .Object@Tobs <- as.integer(0)
            .Object@e <- matrix("numeric")
            .Object$ntvgarch <- list()
            .Object$speedopt <- corrspeedopt$eta
            .Object$optimcontrol <- list(fnscale = -1, reltol = 1e-5)
            .Object$pars <- c(0.05,0.8)
            names(.Object$pars) <- c("alpha","beta")

            # Return:
            .Object
          })

## -- Constructor:dcc -- ####
setGeneric(name="dcc",
           valueClass = "dcc_class",
           signature = c("ntvgarchObj","dcctype"),
           def = function(ntvgarchObj,dcctype){
             this <- new("dcc_class")

             ## -- Do validation checks -- ####

             if(class(ntvgarchObj)[1] != "ntvgarch_class"){
               warning("a valid instance of the ntvgarch_class is required to create an DCC model")
               return(this)
             }
             # End validation

             # Set Default Values based on type:
             # Common:
             N <- this@N <- ntvgarchObj@N
             this@Tobs <- ntvgarchObj@Tobs
             # TODO: Should we make this a variable?  Or inherit from NTVGarch?
             this@st <- (1:ntvgarchObj@Tobs)/ntvgarchObj@Tobs
             this@nr.corPars <- as.integer((N^2-N)/2)
             this@e <- matrix(nrow = this@Tobs,ncol = N)

             # Extract the Data & Estimated components from the ntvgarch
             this$ntvgarch <- list()
             for(n in 1:ntvgarchObj@N){
               this$ntvgarch[[n]] <- list()
               this$ntvgarch[[n]]$tv <- ntvgarchObj[[n]]$Estimated$tv
               this$ntvgarch[[n]]$garch <- ntvgarchObj[[n]]$Estimated$garch
               this@e[,n] <- ntvgarchObj[[n]]@e
             }
             names(this$ntvgarch) <- names(ntvgarchObj)

             # Filter the data:
             z <- w <- e <- this@e
             for(n in 1:this@N){
               w[,n] <- e[,n]/sqrt(this$ntvgarch[[n]]$tv@g)
               z[,n] <- w[,n]/sqrt(this$ntvgarch[[n]]$garch@h)
             }

             if(dcctype == dcctype$General){
               this$Qbar <- toeplitz(0.5^seq.int(0, N-1))
             }else if (dcctype == dcctype$Dynamic){
               this$Q1 <- diag(N)
               this$Q2 <- toeplitz(0.5^seq.int(0, N-1))
               this$speed <- 3
               this$loc <- 0.5

             }else{
               # Invalid type passed
               warning("dcctype not recognised.  Try using the dcctype$...  as a second parameter.")
               return(this)
             }

             return(this)
           }
)

## -- calc.Gt.eta(speed,loc,Tobs) -- ##
calc.Gt.eta = function(speed,loc,Tobs){
  st <- (1:Tobs)/Tobs
  st_c <- st - loc
  G <- 1/(1+exp(-exp(speed)*st_c))
  return(matrix(G,nrow = Tobs,ncol = 1))
}

## -- calc.Pt -- ####
setGeneric(name=".calc.Pt",
           valueClass = "matrix",
           signature = c("dccObj"),
           def =   function(dccObj){
             this <- dccObj

             vP1 <- vecL(this$Estimated$P1)
             vP2 <- vecL(this$Estimated$P2)

             Gt <- calc.Gt(this)
             Pt <- apply(Gt,MARGIN = 1,FUN = function(X,P1,P2) ((1-X)*P1 + X*P2), P1=vP1, P2=vP2)

             if(is.vector(Pt)) Pt <- matrix(Pt,ncol = 1) else Pt <- t(Pt)

             return(Pt)

           }
)


## -- loglik.dcc() --####
setGeneric(name=".loglik.dcc",
           valueClass = "numeric",
           signature = c("optimpars","z","dccObj"),
           def = function(optimpars,z,dccObj){

             err_output <- -1e10
             this <- dccObj

             # optim pars for DCC are:
             # alpha
             # beta
             # +
             # dcctype=General:
             # Qbar_pars
             #  -- or --
             # dcctype=Dynamic:
             # Q1_pars
             # Q2_pars
             # tr.pars = (speed,loc)
             #


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
             maxSpeed <- switch(this$speedopt,1000,(1000/sd(this@st)),7.0,0.30)
             if (speed > maxSpeed) return(err_output)
             if (speed < 0) return(err_output)

             # Check 3: Check the locations fall within min-max values of st
             loc1 <- this$Estimated$pars[2]
             if(this$shape == corrshape$double) loc2 <- this$Estimated$pars[3] else loc2 <- NA
             if (loc1 < min(this@st)) return(err_output)
             if (loc1 > max(this@st)) return(err_output)
             if (!is.na(loc2)) {
               if (loc2 < min(this@st)) return(err_output)
               if (loc2 > max(this@st)) return(err_output)
             }


             #### ======== calculate loglikelihood ======== ####
             Pt <- .calc.Pt(this)  # T x nr.corPars

             llt <- vector("numeric")
             for(t in 1:this@Tobs) {
               mPt <- unVecL(Pt[t,,drop=FALSE])
               llt[t] <- -0.5*log(det(mPt)) -0.5*( z[t,,drop=FALSE] %*% (solve(mPt)) %*% t(z[t,,drop=FALSE]) )
             }
             # Return:
             return(sum(llt))

           }
)

## --- estimateDCC --- ####
setGeneric(name="estimateDCC",
           valueClass = "dcc_class",
           signature = c("dccObj","estimationCtrl"),
           def = function(dccObj,estimationCtrl){
             this <- dccObj
             e <- this@e

             calcSE <- estimationCtrl$calcSE
             verbose <- estimationCtrl$verbose

             this$Estimated <- list()

             optimpars <- c( vecL(this$P1), vecL(this$P2), this$pars )
             optimpars <- optimpars[!is.na(optimpars)]

             ### ---  Call optim to calculate the estimate --- ###
             if (verbose) this$optimcontrol$trace <- 10

             # Filter the data:
             z <- w <- e <- this@e
             for(n in 1:this@N){
               w[,n] <- e[,n]/sqrt(this$ntvgarch[[n]]$tv@g)
               z[,n] <- w[,n]/sqrt(this$ntvgarch[[n]]$garch@h)
             }

             tmp <- NULL
             try(tmp <- optim(optimpars,.loglik.dcc,z,this,gr=NULL,method="BFGS",control=this$optimcontrol,hessian=calcSE))

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
               if(this$shape != corrshape$double) this$Estimated$pars <- c(this$Estimated$pars,NA)
               names(this$Estimated$pars) <- names(this$pars)

               if (calcSE) {
                 cat("\nCalculating STCC standard errors...\n")
                 this$Estimated$hessian <- tmp$hessian
                 vecSE <- vector("numeric")
                 try(vecSE <- sqrt(-diag(qr.solve(tmp$hessian))))

                 if(length(vecSE) > 0) {
                   this$Estimated$P1.se <- unVecL(vecSE[1:this@nr.corPars])
                   vecSE <- tail(vecSE,-this@nr.corPars)
                   this$Estimated$P2.se <- unVecL(vecSE[1:this@nr.corPars])
                   vecSE <- tail(vecSE,-this@nr.corPars)

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
## --- estimateDCC --- ####
setMethod("estimateDCC",signature = c("dcc_class","missing"),
          function(dccObj){
            estimationControl <- list(calcSE = TRUE,verbose = TRUE)
            estimateDCC(dccObj,estimationControl)
          })

## --- unCorrelateData --- ####
setGeneric(name="unCorrelateData",
           valueClass = "matrix",
           signature = c("e","dccObj"),
           def = function(e,dccObj){
             this <- dccObj

             # Return matrix 'u' of uncorrelated data:
             u <- matrix(0,nrow = NROW(e),ncol = NCOL(e))

             Pt <- this$Estimated$Pt
             for(t in 1:this@Tobs){
               P <- unVecL(Pt[t,,drop=FALSE])
               P.sqrt <- sqrt_mat1(P)
               u[t,] <- t(solve(P.sqrt) %*% t(e[t,,drop=FALSE]))
             }

             return(u)
           }
)



