
## -- generateRefData -- ####
## generates Reference Data with TV, GARCH processes, optionally with Correlation
setGeneric(name="generateRefData",
           valueClass = "matrix",
           signature = c("nr.series","nr.obs","tvObj","garchObj","corrObj","noiseDist"),
           def =  function(nr.series,nr.obs,tvObj,garchObj,corrObj,noiseDist)
           {
             ## TODO: 1. Override function to handle TV & GARCH as optional params

             ## Note:
             ## noiseDist is a named-list describing the error-distribution and parameters
             ## e.g. noiseDist$name = 'Normal'     noiseDist$mean = 0  noiseDist$sd = 1
             ## or   noiseDist$name = 'Student-t'  noiseDist$df = 6    noiseDist$ncp = 0


             # # used for debugging
             # nr.series=N; nr.obs=Tobs; tvObj=simTv
             # garchObj=simGarch; corrObj=simCorr
             # noiseDist=errDist;seed=1


             # Generate Noise Data:
             #set.seed(NULL)     # Reset the RNG process
             df=0
             if(isTRUE(toupper(substr(trimws(noiseDist$name),1,1)) =="S")) {
               # Standardised Student-t Error/Noise Distribution
               if(!is.null(noiseDist$df)) {
                 df <- noiseDist$df
               } else {
                 df <- 6
                 message("param noiseDist did not contain a valid $df value.  df set = 6 ")
               }
               u <- matrix(rt(nr.obs * nr.series,df),nrow=nr.obs, ncol=nr.series)
               u <- u*sqrt((df-2)/df)
             }else{
               # Normal Error/Noise Distribution (Default is Standard-Normal)
               u <- matrix(rnorm(nr.obs * nr.series),nrow=nr.obs, ncol=nr.series)
             }

             e <- u

             # Step1: Create Correlation (if required)

             if (!is.null(corrObj)){

               # - - - CCC - - -
               if (isTRUE(class(corrObj) == "ccc_class")){
                 if(is.null(corrObj$Estimated)) {
                   P <- corrObj$P
                 } else P <- corrObj$Estimated$P
                 P.sqrt <-  sqrt_mat1(P)
                 e <- t(P.sqrt %*% t(u))
               }

               if (isTRUE(class(corrObj) == "stcc1_class")){
                 if(is.null(corrObj$Estimated)) {
                   corrObj$Estimated$P1 <- corrObj$P1
                   corrObj$Estimated$P2 <- corrObj$P2
                   corrObj$Estimated$pars <- corrObj$pars
                 }
                 Pt <- .calc.Pt(corrObj)
                 for (t in 1:corrObj@Tobs){
                   mPt <- unVecL(Pt[t,,drop=FALSE])
                   mPt.sqrt <- sqrt_mat1(mPt)
                   e[t,] <- t( mPt.sqrt %*% t(u[t,,drop=FALSE]) )
                 }
               }

               #End: Generate Correlated Data

             }

             # Step2: Inject GARCH into Data

             if (!is.null(garchObj)){
               # Generate Discard Data:
               discardObs <- 1500
               if(isTRUE(toupper(substr(trimws(noiseDist$name),1,1)=="S")) ){
                 # Standardised Student-t Error/Noise Distribution
                 discardData <- matrix(rt(discardObs*nr.series,df),nrow=discardObs, ncol=nr.series)
                 discardData <- discardData*sqrt((df-2)/df)
               } else {
                 # Normal Error/Noise Distribution (Default is Standard-Normal)
                 discardData <- matrix(rnorm(discardObs*nr.series),nrow=discardObs, ncol=nr.series)
               }

               garchObj$pars["omega",1] <- ( 1 - garchObj$pars["alpha",1] - garchObj$pars["beta",1] )
               e <- rbind(discardData,e)
               endRow <- discardObs + nr.obs

               for (b in 1:nr.series){
                 w <- z <- e[,b]
                 ht_1 <- 1
                 w[1] <- z[1]
                 for (t in 2:endRow) {
                   ht <- garchObj$pars["omega",1] + garchObj$pars["alpha",1]*(w[t-1])^2 + garchObj$pars["beta",1]*ht_1
                   if(garchObj$type == garchtype$gjr) { ht <- ht + garchObj$pars["gamma",1]*(min(w[t-1],0)^2) }
                   ht_1 <- ht
                   w[t] <- sqrt(ht)*z[t]
                 }
                 e[,b] <- as.numeric(w)
               }

               # Discard the first 2000
               startRow <- discardObs + 1
               e <- e[(startRow:endRow), ]

             }

             if (!is.null(tvObj)){
               # Step3: Inject TV into Data
               gt <- get_g(tvObj)
               e <- e*sqrt(gt)

             }

             #Return:
             e

           }
)

setMethod("generateRefData",
          signature = c(nr.series="numeric",nr.obs="numeric",tvObj="tv_class",garchObj="garch_class",corrObj="ccc_class",noiseDist="missing"),
          function(nr.series,nr.obs,tvObj,garchObj,corrObj){
            # Default Std Normal Noise Dist:
            noiseDist = list(name='Normal', mean=0, sd=1)
            generateRefData(nr.series,nr.obs,tvObj,garchObj,corrObj,noiseDist)
          }

)

setMethod("generateRefData",
          signature = c(nr.series="numeric",nr.obs="numeric",tvObj="tv_class",garchObj="garch_class",corrObj="stcc1_class",noiseDist="missing"),
          function(nr.series,nr.obs,tvObj,garchObj,corrObj){
            # Default Std Normal Noise Dist:
            noiseDist = list(name='Normal', mean=0, sd=1)
            generateRefData(nr.series,nr.obs,tvObj,garchObj,corrObj,noiseDist)
          }

)

setMethod("generateRefData",
          signature = c(nr.series="numeric",nr.obs="numeric",tvObj="tv_class",garchObj="garch_class",corrObj="missing",noiseDist="missing"),
          function(nr.series,nr.obs,tvObj,garchObj){
            # Default Std Normal Noise Dist:
            noiseDist = list(name='Normal', mean=0, sd=1)
            generateRefData(nr.series,nr.obs,tvObj,garchObj,NULL,noiseDist)
          }

)

setMethod("generateRefData",
          signature = c(nr.series="numeric",nr.obs="numeric",tvObj="tv_class",garchObj="missing",corrObj="missing",noiseDist="missing"),
          function(nr.series,nr.obs,tvObj){
            # Default Std Normal Noise Dist:
            noiseDist = list(name='Normal', mean=0, sd=1)
            generateRefData(nr.series,nr.obs,tvObj,NULL,NULL,noiseDist)
          }

)

setMethod("generateRefData",
          signature = c(nr.series="numeric",nr.obs="numeric",tvObj="missing",garchObj="garch_class",corrObj="missing",noiseDist="missing"),
          function(nr.series,nr.obs,garchObj){
            # Default Std Normal Noise Dist:
            noiseDist = list(name='Normal', mean=0, sd=1)
            generateRefData(nr.series,nr.obs,NULL,garchObj,NULL,noiseDist)
          }

)


## -- generateDCCRefData -- ####
generateDCCRefData=function(nr.series,nr.obs,Qbar,a,b)
{
  # Step1: Generate iid Data & create Correlation
  # u = un-correlated data
  u <- matrix(rnorm(nr.obs * nr.series),nrow=nr.obs, ncol=nr.series)

  # e = correlated data (defaults to 'u' if there is no correlation object)
  e <- u

  # - - - DCC - - -
  # pass in Qbar matrix and "a" alpha & "b" beta for DCC
  # need to create more than just nr.obs, add discard amount then drop the discard part before returning
  discard <- 2000
  discardData <- matrix(rnorm(discard * nr.series),nrow=discard, ncol=nr.series)
  u <- rbind(discardData,u)
  e <- u
  endRow <- discard + nr.obs
  # starting from t=2! (set Qt[1]=Qbar)
  Qt_1 <- Qbar
  for (t in 2:endRow){
    Qt <- (1-a-b)*Qbar + a*t(e[t-1,,drop=FALSE]) %*% e[t-1,,drop=FALSE] + b*Qt_1 # N x N
    #scale Qt by inverse of its sqrt diagonals (from front and back) to make it correlation matrix
    Pt <- diag(sqrt(diag(Qt))^(-1),nrow=nr.series, ncol=nr.series) %*% Qt %*% diag(sqrt(diag(Qt))^(-1),nrow=nr.series, ncol=nr.series)  # N x N
    # create DCC correlated data
    Pt.sqrt <- sqrt_mat1(Pt)
    e[t,] <- t( Pt.sqrt %*% t(u[t,,drop=FALSE]) )
    # for the next loop round, store the "previous" Qt
    Qt_1 <- Qt
  }

  startRow <- discard + 1

  #Return e:
  e[(startRow:endRow), ]

}


## -- calc.Gt(spd,loc,Tobs) -- ##
 .calc.Gt = function(speed,loc,Tobs){
   st <- (1:Tobs)/Tobs
   st_c <- st - loc
   G <- 1/(1+exp(-exp(speed)*st_c))
   return(matrix(G,nrow = Tobs,ncol = 1))
 }

## -- generateDynDCCRefData -- ####
generateDynDCCRefData=function(nr.series,nr.obs,Q1,Q2,speed,loc,a,b)
{

  # Step1: Generate iid Data & create Correlation
  # u = un-correlated data
  u <- matrix(rnorm(nr.obs * nr.series),nrow=nr.obs, ncol=nr.series)

  # e = correlated data (defaults to 'u' if there is no correlation object)
  e <- u

  # - - - dynDCC - - -
  # pass in Q1 and Q2 matrices, speed and loc, and "a" alpha & "b" beta for DCC
  # need to create more than just nr.obs, add discard amount then drop the discard part before returning
  discard <- 2000
  discardData <- matrix(rnorm(discard * nr.series),nrow=discard, ncol=nr.series)
  u <- rbind(discardData,u)
  e <- u
  endRow <- discard + nr.obs
  # starting from t=2! (set Qt[1]=Qbar)
  Qt_1 <- Q1
  for (t in 2:endRow){
    if (t > discard){
      Gt <- .calc.Gt(speed,loc,Tobs) # st=timetrend, loc=0.5, gamma s.t. transition takes place between 0.25 and 0.75 of the sample
      Qbar <- (1-Gt) %*% Q1+Gt %*%Q2
    } else Qbar<-Q1
    Qt <- (1-a-b)*Qbar + a*t(e[t-1,,drop=FALSE]) %*% e[t-1,,drop=FALSE] + b*Qt_1 # N x N
    #scale Qt by inverse of its sqrt diagonals (from front and back) to make it correlation matrix
    Pt <- diag(sqrt(diag(Qt))^(-1),nrow=nr.series, ncol=nr.series) %*% Qt %*% diag(sqrt(diag(Qt))^(-1),nrow=nr.series, ncol=nr.series)  # N x N
    # create DCC correlated data
    Pt.sqrt <- sqrt_mat1(Pt)
    e[t,] <- t( Pt.sqrt %*% t(u[t,,drop=FALSE]) )
    # for the next loop round, store the "previous" Qt
    Qt_1 <- Qt
  }

  startRow <- discard + 1

  #Return e:
  e[(startRow:endRow), ]

}


