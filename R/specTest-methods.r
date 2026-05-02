# *************************************************************************************************************** #
#                        Methods for the TV-Specification / MisSpecification Tests                             ####
# *************************************************************************************************************** #


# ---  TESTS ---  ####

{# Help doco: test.LM.TR2(e,tv,ord) ####
  #' @noRd
#' @title
  #' LM TR-Squared Test
  #'
  #' @description
  #' `test.LM.TR2` is a
  #'
  #' @usage test.LM.TR2(e,tvObj,testOrder)
  #'
  #' @param e An estimated TV model
  #' @param tvObj Use Enum, garchtype$...
  #' @param testOrder An integer
  #'
  #' @details
  #' This object
  #'
  #' ```
  #'  testOrder = 1
  #'  testResult = test.LM.TR2(e,myTv,testOrder)
  #' ```
  #'
  #'
  #' @returns A numeric test statistic
  #'
  #' @note
  #' I am a note
  #'
  #'
}

getTestStats = function(e,tvObj,testOrder){
  this <- list()
  this$TR2 <- test.LM.TR2(e,tvObj,testOrder)
  this$Robust <- test.LM.Robust(e,tvObj,testOrder)
  return(this)
}

test.LM.TR2=function(e,tvObj,testOrder){

  # debug:
  # this <- obj_single
  # testOrder = 1

  this <- tvObj

  if(testOrder <= 0){
    message("Cannot execute test with no alternate hypothesis. Please set a valid testOrder")
    return(NaN)
  }

  # Test Method: Regress psi2_1 on 1/gt*(dgdt and dgdt2)
  # 1. Calc derivatives of params dgdt = Tx1 or Tx4 or Tx7 or...
  #    NCOL(dgdt) increases with the order of TV function.
  dgdt <- .dg_dt(this)

  # 2. Calc derivatives of taylor pars (linearised component) under the null
  dgdt2 <- .dg_dt2(this@st,testOrder)

  g <- calculate_g(this)
  X <- cbind(dgdt,dgdt2)/g

  # 3. Invert crossprod(X) to calculate SSR1
  Xinv <- NULL
  try(Xinv <- solve(crossprod(X)))
  if (is.null(Xinv)){
    rm(g,dgdt,dgdt2)
    return(NaN)
  }

  # 4. Calculate psi2_1 to calculate SSR0
  psi2_1 <- matrix(data=(e^2/g-1),nrow = this@Tobs,ncol = 1)

  # 5. Calc the TestStat:
  SSR0 <- sum(psi2_1*psi2_1)    # Scalar
  SSR1 <- sum((psi2_1-X%*%Xinv%*%t(X)%*%psi2_1)^2)

  Result <- this@Tobs*(SSR0-SSR1)/SSR0

  # Tidy up & release memory before returning:
  rm(this,psi2_1,g,X,Xinv,dgdt,dgdt2)

  # Return:
  Result

}


{# Help doco: test.LM.Robust(e,tv,ord) ####
  #' @noRd
#' @title
  #' LM TR-Squared Test - Robust Version
  #'
  #' @description
  #' `test.LM.Robust` is a
  #'
  #' @usage test.LM.Robust(e,tvObj,testOrder)
  #'
  #' @param e An estimated TV model
  #' @param tvObj Use Enum, garchtype$...
  #' @param testOrder An integer
  #'
  #' @details
  #' This object
  #'
  #' ```
  #'  testOrder = 1
  #'  testResult = test.LM.Robust(e,myTv,testOrder)
  #' ```
  #'
  #'
  #' @returns A numeric test statistic
  #'
  #' @note
  #' I am a note
  #'
  #'
}

test.LM.Robust = function(e,tvObj,testOrder){
  this <- tvObj

  if(testOrder <= 0){
    message("Cannot execute test with no alternate hypothesis. Please set a valid testOrder")
    return(NaN)
  }
  # 1. Calc derivatives of params dgdt = Tx1 or Tx4 or Tx7 or...
  #    NCOL(dgdt) increases with the order of TV function.
  dgdt <- .dg_dt(this)

  # 2. Calc derivatives of taylor pars (linearised component) under the null
  dgdt2 <- .dg_dt2(this@st,testOrder)

  g <- calculate_g(this)
  X <- dgdt/g

  # 3. Invert crossprod(X) to calculate SSR1
  Xinv <- NULL
  try(Xinv <- solve(crossprod(X)))
  if (is.null(Xinv)){
    message("error")
    rm(g,X,dgdt,dgdt2)
    return(NaN)
  }

  XXXX <- X%*%Xinv%*%t(X)
  Y <- as.matrix(dgdt2/g)
  W <- as.matrix(Y-XXXX%*%Y)

  #4. Regress 1 on (psi2-1)*w, and compute SSR
  psi2_1 <- as.vector(e^2/g - 1)
  X <- psi2_1*W  #psi2_1 must be a vector for this!!

  #5. Compute test statistic:
  Xinv <- NULL
  try(Xinv <- solve(crossprod(X)))
  if(is.null(Xinv)) {
    message("error")
    rm(psi2_1,g,W,Y,X,XXXXdgdt,dgdt2)
    return(NaN)
  }

  Result <- this@Tobs-sum(diag(this@Tobs)-(X%*%Xinv%*%t(X)))

  # Tidy up & release memory before returning:
  rm(this,psi2_1,g,W,Y,X,XXXX,Xinv,dgdt,dgdt2)

  # Return:
  Result

}

.test.misSpec.Robust =  function(z,r1,r2){

  rtn <- list()
  z2_1 <- z
  Tobs <- NROW(z)

  SSR0 <- t(z2_1) %*% z2_1

  # 2: regress z2_1 on r1~r2, get SSR
  X <- cbind(r1,r2) # T x ncol(r1)+ncol(r2)
  Y <- z2_1  # T x 1
  #tXX <- solve(t(X),X,tol=1e-10)
  #b <- tXX %*% t(X) %*% Y  # vector len=ncol(X)
  b <- solve(t(X) %*% X) %*% t(X) %*% Y  # vector len=ncol(X)
  resid <- Y - (X %*% b)   # T x 1
  SSR1 <- t(resid) %*% resid # 1x1
  # LM stat
  rtn$LM <- as.numeric(Tobs * (SSR0-SSR1)/SSR0)
  rtn$pVal <- as.numeric(pchisq(rtn$LM,df=NCOL(r2),lower.tail=FALSE))

  # ROBUST
  # 1 as above
  # 2 regress r2 on r1, get residual vectors
  resid <- matrix(0,Tobs,NCOL(r2))
  X <- r1
  for (i in 1:NCOL(r2)){
    Y <- r2[,i,drop=FALSE]
    #tXX <- solve(t(X),X,tol=1e-10)
    #b <- tXX %*% t(X) %*% Y  # vector len=ncol(X)
    b <- solve(t(X) %*% X) %*% t(X) %*% Y  # vector len=ncol(X)
    resid[,i] <- Y-(X %*% b)   # T x 1
  }

  # regress 1 on (z2_1)resid, get SSR
  Y <- matrix(1,Tobs,1)
  X <- as.vector(z2_1) * resid
  #tXX <- solve(t(X),X,tol=1e-10)
  #b <- tXX %*% t(X) %*% Y  # vector len=ncol(X)
  b <- solve(t(X) %*% X) %*% t(X) %*% Y  # vector len=ncol(X)
  resid <- Y - (X %*% b)   # T x 1
  SSR <- t(resid) %*% resid # 1x1
  # LM Robust
  rtn$LMrob <- as.numeric(Tobs-SSR)
  rtn$pValrob <- as.numeric(pchisq(rtn$LMrob,df=NCOL(r2),lower.tail=FALSE))

  return(rtn)
}


{# Help Doco: testStatDist() ####
  #' @noRd
#' @title
  #' Generates a test statistic distribution...
  #'
  #' @description
  #' `testStatDist` is a
  #'
  #' @usage testStatDist(refdata,tvObj,reftests,simcontrol)
  #'
  #' @param refdata A matrix
  #' @param tvObj Use Enum, garchtype$...
  #' @param reftests An integer
  #' @param simcontrol An integer
  #'
  #' @details
  #' This object
  #'
  #' ```
  #'  simcontrol$maxTestorder = 2
  #'  testResult = testStatDist(refdata,tvObj,reftests,simcontrol)
  #' ```
  #'
  #'
  #' @returns A list containing
  #'
  #' @note
  #' I am a note
  #'
  #'
}

testStatDist = function(refdata,tvObj,reftests,simcontrol){
  this <- tvObj

  if(is.null(simcontrol)){
    simControl <- list()
    simControl$saveAs <- paste("TestStatDist-",strftime(Sys.time(),format="%Y%m%d-%H%M%S",usetz = FALSE))
    simControl$numLoops <- 1100
    simControl$numCores <- parallel::detectCores() - 1
    simcontrol$maxTestorder <- 3
  }

  # Validate SimControl:
  if(is.null(simcontrol$maxTestorder)){
    warning("\nMaximum Test Order must be a valid number between 1 - 4")
    return(list())
  }else if((simcontrol$maxTestorder < 1) || (simcontrol$maxTestorder > 4) ){
    warning("\nMaximum Test Order must be a valid number between 1 - 4")
    return(list())
  }else if(length(reftests) != simcontrol$maxTestorder) {
    warning("\nMaximum Test Order is mis-matched with the number of Reference Tests provided")
    return(list())
  }

  # 1. Setup the default params
  if(!is.null(simcontrol$numLoops)) numLoops <- simcontrol$numLoops else numLoops <- 1100
  if(!is.null(simcontrol$numCores)) numCores <- simcontrol$numCores else numCores <- detectCores() - 1

  # 2. Load the generated data with Garch and add the 'g' from our TV object
  refdata <- refdata[1:this@Tobs,]*sqrt(this@g)

  # 3. Setup the parallel backend
  Sys.setenv("MC_CORES" = numCores)
  cl <- makeCluster(numCores)
  registerDoParallel(cl, cores = numCores)

  # 4. Set the estimation controls to suppress SE & console output
  estCtrl <- list(calcSE = FALSE, verbose = FALSE)

  # 5. Setup the timer to provide duration feedback
  tmr <- proc.time()
  timestamp(prefix = "\nStarting to build Test Stat Distribution - ",suffix = "\nPlease be patient as this may take a while...\n")

  # 6. Calculate Results for all test orders required (optimised for parallel processing)

  # First estimate all TV objects in parallel

  noGarch <- garch(garchtype$noGarch,1)  # Prevent the estimateTV method from creating a Garch object every time
  listTV <- foreach(a = 1:numLoops, .inorder=FALSE, .packages = "MTVGARCH", .verbose = FALSE) %dopar%{
    estimateTV(refdata[,a],this,noGarch,estCtrl)    # Note: The tv params don't change, only the data changes
  }

  testStats <- foreach(b = 1:numLoops, .inorder=FALSE, .combine=rbind, .verbose = FALSE) %dopar% {
    if (isFALSE(listTV[[b]]$Estimated$error)) {

      foreach(testOrder = 1:simcontrol$maxTestorder, .inorder=FALSE, .combine=c, .verbose = FALSE) %do% {
        if(is.nan(reftests[[testOrder]]$TR2)) simTEST1 <- NA else simTEST1 <- test.LM.TR2(refdata[,b],listTV[[b]],testOrder)
        if(is.nan(reftests[[testOrder]]$Robust)) simTEST2 <- NA else simTEST2 <- test.LM.Robust(refdata[,b],listTV[[b]],testOrder)
        c(b,reftests[[testOrder]]$TR2,simTEST1,as.integer(simTEST1 > reftests[[testOrder]]$TR2),reftests[[testOrder]]$Robust,simTEST2,as.integer(simTEST2 > reftests[[testOrder]]$Robust),listTV[[b]]$Estimated$value)
        #runSimrow <- c(runSimrow,b,reftests[[testOrder]]$TR2,simTEST1,as.integer(simTEST1 > reftests[[testOrder]]$TR2),reftests[[testOrder]]$Robust,simTEST2,as.integer(simTEST2 > reftests[[testOrder]]$Robust),listTV[[b]]$Estimated$value)
      }
      #gc()  # Release memory locked by tests
    }
  }


  # Debug: Internal Parallel Result:
  #  runSimrow
  # End: testStats <- foreach(b = 1:numloops,...

  # Extract Test P_Values from Results & express as decimal
  colnamesResults <- vector("character")
  if(simcontrol$maxTestorder >= 1) colnamesResults <- c(colnamesResults,"TestOrd_1","Ref$LMTR2","Stat_TR2.1","Pval_TR2.1","Ref$LMRobust","Stat_Robust.1","Pval_Robust.1","Estimated_LL")
  if(simcontrol$maxTestorder >= 2) colnamesResults <- c(colnamesResults,"TestOrd_2","Ref$LMTR2","Stat_TR2.2","Pval_TR2.2","Ref$LMRobust","Stat_Robust.2","Pval_Robust.2","Estimated_LL")
  if(simcontrol$maxTestorder >= 3) colnamesResults <- c(colnamesResults,"TestOrd_3","Ref$LMTR2","Stat_TR2.3","Pval_TR2.3","Ref$LMRobust","Stat_Robust.3","Pval_Robust.3","Estimated_LL")
  if(simcontrol$maxTestorder >= 4) colnamesResults <- c(colnamesResults,"TestOrd_4","Ref$LMTR2","Stat_TR2.4","Pval_TR2.4","Ref$LMRobust","Stat_Robust.4","Pval_Robust.4","Estimated_LL")
  #
  colnames(testStats) <- colnamesResults

  Results <- list()
  if(simcontrol$maxTestorder >= 1){
    Results$Order1 <- list()
    Results$Order1$pVal_TR2 <- mean(testStats[,"Pval_TR2.1"],na.rm = TRUE)
    Results$Order1$pVal_ROB <- mean(testStats[,"Pval_Robust.1"],na.rm = TRUE)
  }
  #
  if(simcontrol$maxTestorder >= 2){
    Results$Order2 <- list()
    Results$Order2$pVal_TR2 <- mean(testStats[,"Pval_TR2.2"],na.rm = TRUE)
    Results$Order2$pVal_ROB <- mean(testStats[,"Pval_Robust.2"],na.rm = TRUE)
  }
  #
  if(simcontrol$maxTestorder >= 3){
    Results$Order3 <- list()
    Results$Order3$pVal_TR2 <- mean(testStats[,"Pval_TR2.3"],na.rm = TRUE)
    Results$Order3$pVal_ROB <- mean(testStats[,"Pval_Robust.3"],na.rm = TRUE)
  }
  #
  if(simcontrol$maxTestorder >= 4){
    Results$Order4 <- list()
    Results$Order4$pVal_TR2 <- mean(testStats[,"Pval_TR2.4"],na.rm = TRUE)
    Results$Order4$pVal_ROB <- mean(testStats[,"Pval_Robust.4"],na.rm = TRUE)
  }
  Results$TestStatDist <- testStats

  # 7. Save the distribution
  if(!is.null(simcontrol$saveAs)) {
    # Create SimDist folder (if not there) & set Save filename
    if (!dir.exists(file.path(getwd(),"SimDist"))) dir.create(file.path(getwd(),"SimDist"))
    saveAs <- file.path("SimDist",simcontrol$saveAs)
    try(saveRDS(Results,saveAs))
  }

  # 8. Stop the parallel cluster
  stopCluster(cl)
  rm(cl)

  # 9. Print the time taken to the console:
  cat("\nTest Stat Distribution Completed \nRuntime:",(proc.time()-tmr)[3],"seconds\n")

  # 10. Attempt to release memory:
  rm(refdata,listTV,testStats)

  # Return:
  return(Results)

}

.dg_dt =  function(tvObj){
  this <- tvObj

  rtn <- matrix(nrow=this@Tobs,ncol=this@nr.pars)
  col_idx <- 0

  if(isTRUE(this@delta0free)){
    col_idx <- col_idx + 1
    rtn[,col_idx] <- 1  # derivative of delta0
  }

  if (this@nr.transitions > 0) {
    # initialise some variables
    stdev_st <- sd(this@st)
    st_c <- speed_transf <- Gi <- 0

    for (i in 1:this@nr.transitions) {

      if(this$shape[i] == tvshape$single) st_c <- this@st - this$Estimated$pars["locN1",i]
      if(this$shape[i] == tvshape$double) st_c <- (this@st - this$Estimated$pars["locN1",i]) * (this@st - this$Estimated$pars["locN2",i])
      if(this$shape[i] == tvshape$double1loc) st_c <- (this@st - this$Estimated$pars["locN1",i])^2

      if(this$speedopt == speedopt$gamma) {
        speed_transf <- this$Estimated$pars["speedN",i]
        Gi <- 1/(1+exp(-this$Estimated$pars["speedN",i] * st_c))
      }
      if(this$speedopt == speedopt$gamma_std) {
        speed_transf <- this$Estimated$pars["speedN",i]/stdev_st
        Gi <- 1/(1+exp(-this$Estimated$pars["speedN",i] * st_c/stdev_st))
      }
      if(this$speedopt == speedopt$eta) {
        speed_transf <- exp(this$Estimated$pars["speedN",i])
        Gi <- 1/(1+exp(-exp(this$Estimated$pars["speedN",i]) * st_c))
      }

      deriv_const <- this$Estimated$pars["deltaN",i]*speed_transf*Gi*(1-Gi)

      col_idx <- col_idx + 1
      rtn[,col_idx] <- Gi    # derivative of delta1..n
      col_idx <- col_idx + 1
      rtn[,col_idx] <- deriv_const*st_c    # derivative of speed1..n

      if(this$shape[i] == tvshape$single){
        col_idx <- col_idx + 1
        rtn[,col_idx] <- -deriv_const    # derivative of loc1..n (shape=TVshape$single)
      }
      if(this$shape[i] == tvshape$double){
        col_idx <- col_idx + 1
        rtn[,col_idx] <- -deriv_const*(this@st-this$Estimated$pars["locN1",i])  # derivative of loc1..n (shape=TVshape$double)
        col_idx <- col_idx + 1
        rtn[,col_idx] <- -deriv_const*(this@st-this$Estimated$pars["locN2",i])  # derivative of loc2..n (shape=TVshape$double)
      }
      if(this$shape[i] == tvshape$double1loc){
        col_idx <- col_idx + 1
        rtn[,col_idx] <- -deriv_const*2*(this@st-this$Estimated$pars["locN1",i])    # derivative of loc1..n (shape=TVshape$double1loc)
      }

    } # End: for loop

  } # End: if (this@nr.transitions > 0)

  return(rtn)

}

.dg_dt2 =  function(st,testOrder){

  rtn <- matrix(nrow=NROW(st),ncol=testOrder)
  for(n in 1:testOrder){
    rtn[,n] <- st^n
  }
  # Return:
  return(rtn)
}

.dh_dg =  function(e,tvgarchObj){
  this <- tvgarchObj
  Tobs <- this@Tobs
  w <- e/sqrt(this$Estimated$tv$g)

  dhdg <- matrix(NA,nrow = Tobs,ncol = this@tvObj@nr.pars)  # Initialise return matrix = T x TV@nr.pars, value = 1

  if (this$garch$type == garchtype$noGarch){
    return(dhdg)
  } else if (this$garch$type == garchtype$general){
    v_garch <- cbind(c(0,rep(1,(Tobs-1))),c(0,w[1:(Tobs-1)]^2),c(0,this$Estimated$garch$h[1:(Tobs-1)])) # T x Num_garch_pars, each row = "1~w(i,t-1)^2~h(i,t-1)", i=1,...,N
  } else if (this$garch$type == garchtype$gjr){
    v_garch <- cbind(c(0,rep(1,(Tobs-1))),c(0,w[1:(Tobs-1)]^2),c(0,this$Estimated$garch$h[1:(Tobs-1)]),c(0,(min(w[1:(Tobs-1)],0))^2)) # T x Num_garch_pars, each row = "1~w(i,t-1)^2~h(i,t-1)", i=1,...,N
  }

  beta <- rep(this$Estimated$garch$Estimated$pars["beta",1],NCOL(v_garch))
  dhdg <- .ar1.Filter(v_garch,beta) # T x Num_garch_pars, each row = dh(t).dgarchpars
  return(dhdg)

}

.dh_dt =  function(e,tvgarchObj){
  this <- tvgarchObj
  Tobs <- this@tvObj@Tobs
  w <- e/sqrt(this$Estimated$tv$g)

  dhdt <- matrix(NA,nrow = Tobs,ncol = this@tvObj@nr.pars)  # Initialise return matrix = T x TV@nr.pars, value = 1

  dgdt <- .dg_dt(this$Estimated$tv)  # T x TV@nr.pars (includes d0's derivative if it is a free param)
  dgdt_l1 <- rbind(rep(0,NCOL(dgdt)),dgdt[(1:Tobs-1),])
  g <- this$Estimated$tv$g

  if (this$garch$type != garchtype$noGarch){
    v_tv <- (-this$Estimated$garch$Estimated$pars["alpha",1] * c(0,1/g[1:(Tobs-1)])*c(0,w[1:(Tobs-1)]^2)) * dgdt_l1
    beta <- rep(this$Estimated$garch$Estimated$pars["beta",1],NCOL(v_tv))
    dhdt <- .ar1.Filter(v_tv,beta) # T x Num_tv_pars, each row = dh(t).dtvpar
  }
  return(dhdt)
}

.df_dli =  function(tvObj,testOrder){
  this <- tvObj

  st <- this@st

  if(isTRUE(this@delta0free)){
    # We already have a derivative of the constant
    nrDerivs <- testOrder
    ret <- matrix(nrow=NROW(st),ncol=nrDerivs)
    for(n in 1:nrDerivs){ ret[,n] <- st^(n) }

  }else{

    nrDerivs <- testOrder + 1
    ret <- matrix(nrow=NROW(st),ncol=nrDerivs)
    for(n in 1:nrDerivs){ ret[,n] <- st^(n-1) }
  }
  return(ret)
}


{# -- test.misSpec1(e,tvgarch,testOrd) ####
#' @noRd
#' @title
#' Tests for mis-specification of...
#'
#' @description
#' `test.misSpec1` is a
#'
#' @usage test.misSpec1(e,tvgarch,testOrd)
#'
#' @param e A numeric vector
#' @param tvgarch An estimated TvGarch model
#' @param testOrd Use Enum, garchtype$...
#'
#' @details
#' This object
#'
#'
#'
#' @returns A list...
#'
#' @note
#' I am a note
#'
#'
#'
}

## Tests for mis-specification of the tv component of the model, e.g. is there another transition?
test.misSpec1 =  function(e,tvgarchObj,testOrder){
  this <- tvgarchObj
  rtn <- list()

  Tobs <- this@Tobs
  g <- this$Estimated$tv$g
  h <- this$Estimated$garch$h

  # derivatives:
  dg_dtv <- .dg_dt(this$Estimated$tv)        # T x nr.tv.pars
  dh_dtv <- .dh_dt(e,this)                   # T x nr.tv.pars
  dh_dga <- .dh_dg(e,this)                   # T x nr.garch.pars
  df_dli <- .df_dli(this@tvObj,testOrder)    # T x (testorder+1)

  z <- e/sqrt(g*h)
  z2 <- z^2
  z2_1 <- z^2 - 1

  # matrices
  u <- g^(-1) * dg_dtv     # (1/g)*(dg/dtvpars); T x nr.tv.pars
  x1 <- (h^(-1)) * dh_dtv  # (1/ht)*(dh/dtvpars); T x nr.tv.pars
  x2 <- (h^(-1)) * dh_dga  # (1/ht)*(dh/dgarchpars); T x nr.garch.pars
  #
  if(this$Estimated$garch$type == garchtype$noGarch){
    r1 <- u     # T x nr.tv.pars
  }else{
    r1 <- cbind(u+x1,x2)     # T x (tot.tv+garch.pars)
  }
  r2 <- (g^(-1)) * df_dli  # (1/gt)*(df/dlinpars); T x nr.lin.pars (=testorder+1)
  rtn <- .test.misSpec.Robust(z2_1,r1,r2)
  return(rtn)
}


{# -- test.misSpec2(e,tvgarch,type) ####
#' @noRd
#' @title
#' Tests for mis-specification of...
#'
#' @description
#' `test.misSpec2` is a
#'
#' @usage test.misSpec2(e,tvgarch,type)
#'
#' @param e A numeric vector
#' @param tvgarch An estimated TvGarch model
#' @param type Use Enum, garchtype$...
#'
#' @details
#' This object
#'
#'
#'
#' @returns A list...
#'
#' @note
#' I am a note
#'
#'
}

## Tests for mis-specification of the Garch-order of the model
test.misSpec2 =  function(e,tvgarchObj,type){
  # type 1: GARCH(1,1) vs GARCH(1,2)
  # type 2: GARCH(1,1) vs GARCH(2,1)
  # type 3: GARCH vs STGARCH

  this <- tvgarchObj
  rtn <- list()

  Tobs <- this@Tobs
  g <- this$Estimated$tv$g
  h <- this$Estimated$garch$h
  w <- e/sqrt(g)
  z <- e/sqrt(g*h)
  z2 <- z^2
  z2_1 <- z^2 - 1

  # derivatives:
  dg_dtv <- .dg_dt(this$Estimated$tv)     # T x nr.tv.pars
  dh_dtv <- .dh_dt(e,this)                # T x nr.tv.pars
  dh_dga <- .dh_dg(e,this)                # T x nr.garch.pars

  # matrices
  u <- g^(-1) * dg_dtv     # (1/g)*(dg/dtvpars); T x nr.tv.pars
  x1 <- (h^(-1)) * dh_dtv  # (1/ht)*(dh/dtvpars); T x nr.tv.pars
  x2 <- (h^(-1)) * dh_dga  # (1/ht)*(dh/dgarchpars); T x nr.garch.pars

  r1 <- cbind(u+x1,x2)     # T x (tot.tv+garch.pars)
  if (type==1){ r2 <- (h^(-1)) * lag0(w^2,2)  }                 # T x 1
  if (type==2){ r2 <- (h^(-1)) * lag0(h,2)    }                 # T x 1
  if (type==3){ r2 <- (h^(-1)) * cbind(lag0(w,1),lag0(w^3,1)) } # T x 2

  rtn <- .test.misSpec.Robust(z2_1,r1,r2)
  return(rtn)
}


{# -- test.misSpec3(e,tvgarch,maxLag) ####
#' @noRd
#' @title
#' Tests for mis-specification of...
#'
#' @description
#' `test.misSpec3` is a
#'
#' @usage test.misSpec3(e,tvgarch,maxLag)
#'
#' @param e A numeric vector
#' @param tvgarch An estimated TvGarch model
#' @param maxLag Use Enum, garchtype$...
#'
#' @details
#' This object
#'
#'
#'
#' @returns A list...
#'
#' @note
#' I am a note
#'
#'
}

## Tests for mis-specification of the remaining-ARCH in the model

test.misSpec3 =  function(e,tvgarchObj,maxLag){
  this <- tvgarchObj
  rtn <- list()

  Tobs <- this@tvObj@Tobs
  g <- this$Estimated$tv$g
  h <- this$Estimated$garch$h
  z <- e/sqrt(g*h)
  z2 <- z^2
  z2_1 <- z^2 - 1

  # derivatives:
  dg_dtv <- .dg_dt(this$Estimated$tv)     # T x nr.tv.pars
  dh_dtv <- .dh_dt(e,this)                # T x nr.tv.pars
  dh_dga <- .dh_dg(e,this)                # T x nr.garch.pars

  # matrices
  u <- g^(-1) * dg_dtv         # (1/g)*(dg/dtvpars); T x nr.tv.pars
  x1 <- (h^(-1)) * dh_dtv      # (1/ht)*(dh/dtvpars); T x nr.tv.pars
  x2 <- (h^(-1)) * dh_dga      # (1/ht)*(dh/dgarchpars); T x nr.garch.pars
  #
  if(this$Estimated$garch$type == garchtype$noGarch){
    r1 <- u     # T x nr.tv.pars
  }else{
    r1 <- cbind(u+x1,x2)     # T x (tot.tv+garch.pars)
  }
  #
  r2 <- lag0(z2,(1:maxLag))    # T x maxlag

  rtn <- .test.misSpec.Robust(z2_1,r1,r2)
  return(rtn)
}


{# -- test.misSpec4(e,tvgarch,exogVar) ####
#' @noRd
#' @title
#' Tests for mis-specification of...
#'
#' @description
#' `test.misSpec4` is a
#'
#' @usage test.misSpec4(e,tvgarch,exogVar)
#'
#' @param e A numeric vector
#' @param tvgarch An estimated TvGarch model
#' @param exogVar Use Enum, garchtype$...
#'
#' @details
#' This object
#'
#'
#'
#' @returns A list...
#'
#' @note
#' I am a note
#'
#'
}

## Tests for mis-specification of the additive exogenous variable(s) in the model
test.misSpec4 =  function(e,tvgarchObj,exogVar){
  this <- tvgarchObj
  rtn <- list()

  Tobs <- this@Tobs
  g <- this$Estimated$tv$g
  h <- this$Estimated$garch$h

  # derivatives:
  dg_dtv <- .dg_dt(this$Estimated$tv)        # T x nr.tv.pars
  dh_dtv <- .dh_dt(e,this)                   # T x nr.tv.pars
  dh_dga <- .dh_dg(e,this)                   # T x nr.garch.pars

  z <- e/sqrt(g*h)
  z2 <- z^2
  z2_1 <- z^2 - 1

  # matrices
  u <- g^(-1) * dg_dtv     # (1/g)*(dg/dtvpars); T x nr.tv.pars
  x1 <- (h^(-1)) * dh_dtv  # (1/ht)*(dh/dtvpars); T x nr.tv.pars
  x2 <- (h^(-1)) * dh_dga  # (1/ht)*(dh/dgarchpars); T x nr.garch.pars
  #
  if(this$Estimated$garch$type == garchtype$noGarch){
    r1 <- u     # T x nr.tv.pars
  }else{
    r1 <- cbind(u+x1,x2)     # T x (tot.tv+garch.pars)
  }
  #
  r2 <- (g^(-1)) * exogVar  # (1/gt)*(exogVar); T x nr. Exogenous Vars

  rtn <- .test.misSpec.Robust(z2_1,r1,r2)
  return(rtn)
}








