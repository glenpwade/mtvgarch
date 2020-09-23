## --- ccc_class Definition --- ####

ccc <- setClass(Class = "ccc_class",
                  slots = c(nr.covPars="integer",Tobs="integer",N="integer"),
                  contains = c("namedList")
)

## --- Initialise --- ####
setMethod("initialize","ccc_class",
          function(.Object,...){
            .Object <- callNextMethod(.Object,...)

            # Default initial values
            .Object@N <- as.integer(0)
            .Object@Tobs <- as.integer(0)
            .Object@nr.covPars <- as.integer(0)

            # Return:
            .Object
          })

## -- Constructor:ccc -- ####
setGeneric(name="ccc",
           valueClass = "ccc_class",
           signature = c("mtvgarchObj"),
           def = function(mtvgarchObj){
             this <- new("ccc_class")

             ## -- Do validation checks -- ####
             objType <- class(mtvgarchObj)
             if(objType[1] != "mtvgarch_class"){
               warning("a valid instance of the mtvgarch_class is required to create a ccc model")
               return(this)
             }
             # End validation

             # Set Default Values:
             this@N <- mtvgarchObj@N
             this@Tobs <- mtvgarchObj@Tobs

             N <- this@N
             this@nr.covPars <- as.integer((N^2-N)/2)
             this$P <- matrix(0,N,N)
             diag(this$P) <- 1

             return(this)
           }
)


setGeneric(name=".loglik.ccc",
           valueClass = "list",
           signature = c("optimpars","z","cccObj"),
           def = function(optimpars,z,cccObj){

             err_output <- -1e10
             this <- list()

             #### ======== constraint checks ======== ####

             # # Check 1: Confirm we have

             vP <- optimpars
             mP <- .unVecl(vP)
             eig <- eigen(mP,symmetric=TRUE,only.values = TRUE)
             if (min(eig$values) <= 0) return(err_output)

             # - - - P(t) and loglik-value
             Tobs <- NROW(z)
             llt <- rep(0,Tobs)
             mPinv <- solve(mP)

             for(t in seq(1,Tobs)) llt[t] <- -0.5*log(det(mP)) - 0.5*(t(z[t,])%*%(mPinv)%*%z[t,])
             Pt <- matrix(vP,nrow=Tobs,ncol=length(vP),byrow=TRUE)

             this$value <- sum(llt)
             this$P <- Pt

             return(this)

           }
)

## --- estimateCCC --- ####
setGeneric(name="estimateCCC",
           valueClass = "ccc_class",
           signature = c("z","cccObj","estimationCtrl"),
           def = function(z,cccObj,estimationCtrl){
             this <- cccObj

             calcSE <- estimationCtrl$calcSE
             verbose <- estimationCtrl$verbose

             this$Estimated <- list()

             optimpars <- c( this$P[lower.tri(this$P)] )

             ### ---  Call optim to calculate the estimate --- ###
             if (verbose) this$optimcontrol$trace <- 10

             this$Estimated <- .loglik.ccc(optimpars,z,this)

             return(this)
           }
)

##====  TESTS  ====####


test.CCCParsim.LM <- function(z,H0,H1,testorder=1) {

  ####--- Initialise ---####
  if(T){

    nGARCH <- H0$nGARCH
    nTV <- H0$nTV
    CCC <- H0$CCC
    ALT <- H1  # the alternative is any CC model, which has a transition, i.e. an st-variable

    if (is.null(CCC)) stop("LM Test: need CCC as input (H0)")
    if (is.null(ALT$st)) stop("LM Test: need transition variable for the test as input (H1)")

    N <- NCOL(z)  # Number of series
    Tobs <- NROW(z)
    st <- ALT$st
    P <- CCC$P
    tmp <- eigen(P)
    Q <- tmp$vectors
    L <- tmp$values # vector
    L.inv <- diag(L^(-1)) # matrix
    L2.inv <- diag(L^(-2)) # matrix
    L <- diag(L) # matrix
    I <- diag(nrow = N,ncol = N) # identity matrix
    One_N1 <- matrix(1,nrow=N,ncol=1)
    One_NN <- matrix(1,nrow=N,ncol=N)
    One_1.N_1 <- matrix(1,nrow=1,ncol=(N-1))
    One_N_1.1 <- matrix(1,nrow=(N-1),ncol=1)
    One_N_1.N_1 <- matrix(1,nrow=(N-1),ncol=(N-1))
    P.inv <- try(solve(P))
    if (!is.matrix(P.inv)) stop("LM Test: Can't invert the correlation matrix P")

    Num_garch_pars <- 0
    Num_tv_pars <- 0

    #Construct the K matrix: N^2 x N^2
    K <- NULL
    for (i in 1:N) {
      # block rows
      Krow <- NULL
      for (j in 1:N) {
        # block columns
        block <- matrix(0,N,N)
        block[j,i] <- 1
        Krow <- cbind(Krow,block)
      }
      K <- rbind(K,Krow)
    }
  }

  # 2. Create the Information Matrix for the Correlation, and extract the SE corner of its inverse:

  ####--- Initialise IM Matrix variables ---####
  if(T){
    g <- matrix(1,Tobs,N) # for now, no TV estimated, so these are 1, T x N
    h <- matrix(1,Tobs,N) # place for estimated GARCH variances, ones at this stage, T x N
    beta <- rep(0,N) # estimated beta coefficients from GARCH equations
    # Extract the number of TV & Garch pars:
    for (n in 1:N) {
      Num_tv_pars <- Num_tv_pars + nTV[[n]]@nr.pars
      g[,n] <- nTV[[n]]@g # univariate TV volatilities, T x N
      if(nGARCH[[n]]$type != garchtype$noGarch) {
        Num_garch_pars <- Num_garch_pars + nGARCH[[n]]@nr.pars
        h[,n] <- nGARCH[[n]]@h        # univariate GARCH volatilities, T x N
        beta[n] <- nGARCH[[n]]$Estimated$pars["beta",1]
      }
    }
    # TODO: Convert for() loops to lapply functions
    #Num_garch_pars <- lapply(nGARCH,sum())
    IM_cor       <- matrix(0, nrow=(testorder+1)*(N-1), ncol=(testorder+1)*(N-1))
    IM_garch     <- matrix(0, nrow=Num_garch_pars,      ncol=Num_garch_pars)
    IM_garch_cor <- matrix(0, nrow=Num_garch_pars,      ncol=(testorder+1)*(N-1))
    IM_tv        <- matrix(0, nrow=Num_tv_pars,         ncol=Num_tv_pars)
    IM_tv_cor    <- matrix(0, nrow=Num_tv_pars,         ncol=(testorder+1)*(N-1))
    IM_tv_garch  <- matrix(0, nrow=Num_tv_pars,         ncol=Num_garch_pars)

    w <- e/sqrt(g) # returns standardised by g(t), T x N
    z <- w/sqrt(h) # returns standardised by GARCH volatilities and TV g(t), T x N
  }

  # partial derivatives of P w.r.t correlation parameters, T x 1 or T x 2
  # these are for alternative parameters, later add column of ones at the front

  ####--- Set x_tau ---####
  if(T){
    if (testorder==1) x_tau <- (-0.5)*cbind(st)
    if (testorder==2) x_tau <- (-0.5)*cbind(st,st^2)
    if (testorder==3) x_tau <- (-0.5)*cbind(st,st^2,st^3)

    # score for rho
    qxq <- matrix(0,nrow=N,ncol=(N^2))  # each row = qi' %x% qi', 1 x N^2, N rows
    for (n in 1:N){
      qxq[n,] <- t(Q[,n,drop=FALSE]%x%Q[,n,drop=FALSE])
    }
    mHelp0 <- matrix(0,nrow=Tobs,ncol=(N-1))
    for (t in 1:Tobs){
      vec_zz_t <- vec(t(z[t,,drop=FALSE])%*%(z[t,,drop=FALSE]))  # N^2 x 1
      mHelp0[t,] <-  L.inv[1:(N-1),1:(N-1)]%*%(One_N_1.1-L.inv[1:(N-1),1:(N-1)]%*%qxq[1:(N-1),,drop=FALSE]%*%vec_zz_t) - One_N_1.1%*%L.inv[N,N]%*%(1-L.inv[N,N]*(qxq[N,,drop=FALSE]%*%vec_zz_t))  # 1 x N-1
    }
    # mHelp0_v3 <- matrix(0,nrow=Tobs,ncol=N-1)
    # for (t in 1:Tobs){
    #   for (n in 1:(N-1)){
    #     mHelp0_v3[t,n] <- 1/L[n,n]-(1/L[n,n]^2)%*%(z[t,,drop=FALSE]%*%Q[,n,drop=FALSE])^2-1/L[N,N]+(1/L[N,N]^2)%*%(z[t,,drop=FALSE]%*%Q[,N,drop=FALSE])^2
    #   }
    # }
    # mHelp0_v2 <- matrix(0,nrow=Tobs,ncol=N-1)
    # for (t in 1:Tobs){
    #    for (n in 1:(N-1)){
    #      mHelp0_v2[t,n] <- (1/L[n,n])-((1/L[n,n])^2)*(t(Q[,n])%*%z[t,])^2-(1/L[N,N])+((1/L[N,N])^2)*t(Q[,N]%*%z[t,])^2
    #    }
    #  }
    # tmp<-rep(0,((testorder)*(N-1)))
    # for (t in 1:Tobs){
    #   tmp <- tmp + t(x_tau[t,])%x%t(mHelp0[t,])
    # }

    dlldrho_A <- t(mHelp0)%*%x_tau   # (N-1)xT %*% (T x testorder) = N-1 x testorder, SUM OVER TIME
    dlldrho_A <- vec(dlldrho_A)  # testorder*(N-1) x 1, SUM OVER TIME
    x_tau <- cbind(rep(-0.5,Tobs),x_tau) # T x 2 or T x 3, now add column of ones at the front (Note:xtau includes -0.5)
    # Check the size:
    if(NCOL(x_tau) != testorder+1) stop("x_tau is the wrong size")
    if(NROW(dlldrho_A) != testorder*(N-1)) stop("dlldrho_A is the wrong size")
  }

  ####--- partial derivatives for GARCH and "x" matrices for GARCH parts ---####
  if(T){
    # partial derivatives of h1,h2,...,hN w.r.t garch_pars
    v_garch <- NULL
    h_scale <- NULL
    beta_scale <- NULL
    for (n in 1:N){
      if(nGARCH[[n]]$type==garchtype$noGarch ) {
        # Do nothing, move on to next series
      } else if(nGARCH[[n]]$type==garchtype$general) {
        # standard GARCH(1,1) case
        v_garch <- cbind(v_garch,c(0,rep(1,(Tobs-1))),c(0,w[1:(Tobs-1),n]^2),c(0,h[1:(Tobs-1),n])) # T x Num_garch_pars, each row = "1~w(i,t-1)^2~h(i,t-1)", i=1,...,N
        beta_scale <- beta_scale <- c(beta_scale,(beta[n] %x% rep(1,3)))    # vector, length 3Num_garch_pars
        scaleFactor <- matrix(1,nrow=1,ncol=3)
        h_scale <- cbind(h_scale,h[,n,drop=FALSE] %x% scaleFactor)    # T x Num_garch_pars
      } else if (nGARCH[[n]]$type > garchtype$general) {
        # GJR or other forms not done yet
        print("ARRGH!")
      }
    }
    if(!is.null(v_garch) && !is.null(beta_scale)){
      dhdt <- myFilter(v_garch,beta_scale) # T x Num_garch_pars, each row = "(dh(i,t)/dtheta(i))' ", i=1,...,N
      # matrix containing the x_t's (garch part)
      x_garch <- -0.5*dhdt/(h_scale) # T x Num_garch_pars, each row = "x_it'", i=1,...,N
      # Check the size:
      if(NCOL(dhdt) != Num_garch_pars) stop("dhdt is the wrong size")
      if(NCOL(x_garch) != Num_garch_pars) stop("x_garch is the wrong size")
    } else {
      x_garch <- matrix(nrow=Tobs,ncol=0)
    }
  }

  ####--- partial derivatives for TV and "x" matrices for TV parts ---####

  if(T){
    v_tv <- NULL
    h_scale <- NULL
    g_scale <- NULL
    beta_scale <- NULL
    dgdt <- NULL
    dhdt <- NULL

    # partial derivatives of g1,g2,...,gN w.r.t tv_pars
    for (n in 1:N){
      #TODO: Confirm with Anna - are the bug-fixes below ok?
      dgdt_n <- dg_dt(nTV[[n]])  # T x #of TVpars in TV[n], includes d0's derivative if it is a free param
      dgdt <- cbind(dgdt,dgdt_n)
      scaleFactor <- matrix(1,nrow=1,ncol=NCOL(dgdt_n))
      g_scale <- cbind( g_scale,g[,n,drop=FALSE] %x% scaleFactor ) # cbinds g(nt) as many times as g(n) has tv parameters

      if (nGARCH[[n]]$type==garchtype$noGarch){
        # do nothing
      } else if (nGARCH[[n]]$type==garchtype$general) {
        # standard GARCH(1,1) case
        v_tv_n <- ((-nGARCH[[n]]$Estimated$pars["alpha",1]*c(0,1/g[1:(Tobs-1),n])*c(0,w[1:(Tobs-1),n]^2))%x% scaleFactor) *dgdt_n
        v_tv <- cbind(v_tv,v_tv_n)
        beta_scale_n <- (beta[n] %x% scaleFactor)
        beta_scale <- c(beta_scale,beta_scale_n)
        h_scale_n <- h[,n,drop=FALSE] %x% scaleFactor
        h_scale <- cbind(h_scale,h_scale_n)
      } else if(nGARCH[[n]]$type > garchtype$general) {
        # GJR or other type not done yet!
        stop("GARCH Type not supported by this Test yet")
      }
    } # End: for loop

    if(Num_tv_pars > 0 && Num_garch_pars > 0){
      dhdt <- myFilter(v_tv,beta_scale) # T x Num_tv_pars, each row = dh(i,t).dtvpar(i), i=1...N
      x_tv <- -0.5*dhdt/h_scale -0.5*dgdt/g_scale  # T x Num_tv_pars, each row = x_it, i=1...N
    } else if(Num_tv_pars > 0 && Num_garch_pars == 0){
      x_tv <- -0.5*dgdt/g_scale  # T x Num_tv_pars, each row = x_it, i=1...N
    } else {
      x_tv <- matrix(1,nrow=Tobs,ncol=0)
    }

    # Check Size:
    if(NCOL(x_tv) != Num_tv_pars) stop("x_tv is the wrong size")

  }

  # Information matrix:
  ####--- Blocks involving GARCH only or GARCH and Correlation ---####
  # IM_garch (Num_garch_pars x Num_garch_pars), SUM OVER TIME
  One_31 <- matrix(1,3,1)
  One_33 <- matrix(1,3,3)
  I.P.Pinv <- I + P*P.inv
  I.P.Pinv_scale <- NULL
  scaleFactor <- NULL
  for (i in 1:N) {
    # i = row index
    I.P.Pinv_scale_row <- NULL
    if (nGARCH[[i]]$type == garchtype$noGarch){
      # do nothing
    } else {
      for (j in 1:N) {
        if(nGARCH[[j]]$type == garchtype$noGarch){
          # do nothing
        } else {
          # j = col index
          if(nGARCH[[i]]$type==garchtype$general  && nGARCH[[j]]$type==garchtype$general) scaleFactor <- One_33
          else stop("Garch type is not supported")
          I.P.Pinv_scale_row <- cbind(I.P.Pinv_scale_row, (I.P.Pinv[i,j] %x% scaleFactor))
        }
      } # End: for "j" loop
    }
    I.P.Pinv_scale <- rbind(I.P.Pinv_scale,I.P.Pinv_scale_row)
  } # End: for "i" loop

  if(!is.null(I.P.Pinv_scale)){
    IM_garch <- ((t(x_garch)%*%x_garch) * I.P.Pinv_scale) / Tobs
    # Check Size:
    if(NCOL(IM_garch) != Num_garch_pars) stop("IM_garch is the wrong size")
    if(NROW(IM_garch) != Num_garch_pars) stop("IM_garch is the wrong size")

    # IM_garch_cor, will be 3N x (testorder+1)*(N-1)
    mHelp1 <- matrix(0,nrow=0,ncol=N-1)

    for (n in 1:N) {
      if (nGARCH[[n]]$type==garchtype$general) {
        scaleFactor <- matrix(1,nrow=3,ncol=1)
        mHelp1 <- rbind(mHelp1, scaleFactor %x% (2*Q[n,1:(N-1),drop=FALSE]^2%*%L.inv[1:(N-1),1:(N-1)] - 2*Q[n,N]^2*L.inv[N,N]*One_1.N_1))
      }
    }
    mHelp2 <- t(x_garch) %*% x_tau
    mHelp3 <- matrix(0,nrow=NROW(mHelp2),ncol=(N-1)*NCOL(mHelp2))
    for (i in 1:NROW(mHelp2)){
      mHelp3[i,] <- mHelp2[i,,drop=FALSE]%x%mHelp1[i,,drop=FALSE]
    }

    IM_garch_cor <- mHelp3/Tobs # (Num_Garch_Pars x (testorder+1)*N), SUM OVER TIME
    # Check Size:
    if(NCOL(IM_garch_cor) != (N-1)*(testorder+1)) stop("IM_garch_cor is the wrong size - COLs")
    if(NROW(IM_garch_cor) != Num_garch_pars) stop("IM_garch_cor is the wrong size - ROWs")

  } # End: if(!is.null(I.P.Pinv_scale))  => Only process when we have Garch

  ####--- Blocks involving TV only or TV and Correlation ---####
  #if (nTV$Type != TRshape$none){
  #IM_tv <- NULL
  I.P.Pinv <- I + P*P.inv
  I.P.Pinv_scale <- NULL
  scaleFactor <- NULL
  for (i in 1:N) {
    # i = row index
    ii <- nTV[[i]]@nr.pars
    if (ii==0){
      # do nothing
    } else {
      I.P.Pinv_scale_row <- NULL
      for (j in 1:N) {
        # j = col index
        jj <- nTV[[j]]@nr.pars
        if (jj==0){
          # do nothing
        } else {
          scaleFactor <- matrix(1,ii,jj)
          I.P.Pinv_scale_row <- cbind(I.P.Pinv_scale_row, (I.P.Pinv[i,j] %x% scaleFactor))
        } # End: if (jj==0)
      } # End: for (j in 1:N)
      I.P.Pinv_scale <- rbind(I.P.Pinv_scale,I.P.Pinv_scale_row)
    } # End: if (ii==0)
  } # End: for (i in 1:N)
  IM_tv <- ((t(x_tv)%*%x_tv) * I.P.Pinv_scale) / Tobs  # sum(num_tv_pars) x sum(num_tv_pars)

  # IM_tv_cor, (total #of tv pars in model) x (testorder+1)*(N-1)
  mHelp4 <- matrix(0,nrow=0,ncol=N-1)
  for (n in 1:N) {
    ii <- nTV[[n]]@nr.pars
    if (ii==0){
      # do nothing
    } else {
      scaleFactor <- matrix(1,nrow=ii,ncol=1)
      mHelp4 <- rbind(mHelp4, scaleFactor %x% (2*Q[n,1:(N-1),drop=FALSE]^2%*%L.inv[1:(N-1),1:(N-1)] - 2*Q[n,N]^2*L.inv[N,N]*One_1.N_1) )
    } # End: if (ii==0)
  } # End: for (n in 1:N)
  mHelp5 <- t(x_tv)%*%x_tau
  mHelp6 <- matrix(0,nrow=NROW(mHelp5),ncol=(N-1)*NCOL(mHelp5))
  if (NROW(mHelp6)>0){
    for (i in 1:NROW(mHelp6)){
      mHelp6[i,] <- mHelp5[i,,drop=FALSE] %x% mHelp4[i,,drop=FALSE]
    }
  }

  IM_tv_cor <- mHelp6 / Tobs  # (Num_tv_Pars x (testorder+1)*N),  SUM OVER TIME

  # Check Size:
  if(NCOL(IM_tv) != Num_tv_pars) stop("IM_tv is the wrong size")
  if(NCOL(IM_tv_cor) != (testorder+1)*(N-1)) stop("IM_tv_cor is the wrong size")
  if(NROW(IM_tv_cor) != Num_tv_pars) stop("IM_tv_cor is the wrong size")


  ####--- Blocks involving TV and GARCH ---####
  #if (nGARCH$Type != garchtype$noGarch && nTV$Type != TRshape$none){
  I.P.Pinv <- I + P*P.inv  # N x N
  I.P.Pinv_scale <- NULL
  for (i in 1:N){
    # i = row index
    ii <- nTV[[i]]@nr.pars
    if (ii==0){
      # do nothing
    } else {
      I.P.Pinv_scale_row <- NULL
      for (j in 1:N){
        # j = col index
        jj <- nGARCH[[j]]@nr.pars
        if (jj==0){
          # do nothing
        } else {
          scaleFactor <- matrix(1,ii,jj)
          I.P.Pinv_scale_row <- cbind(I.P.Pinv_scale_row, (I.P.Pinv[i,j] %x% scaleFactor))
        } # End if (jj==0)
      } # End for (j in 1:N)
      I.P.Pinv_scale <- rbind(I.P.Pinv_scale,I.P.Pinv_scale_row)  # (total of tvpars) x (total of garchpars)
    } # End if (ii==0)
  }
  #IM_tv_garch <- NULL
  IM_tv_garch <- (t(x_tv)%*%x_garch)*(I.P.Pinv_scale)/Tobs  #  (total # of tvpars)x(total # of garchpars)

  # Check Size:
  if(NCOL(IM_tv_garch) != Num_garch_pars) stop("IM_tv_garch is the wrong size")
  if(NROW(IM_tv_garch) != Num_tv_pars) stop("IM_tv_garch is the wrong size")

  #}

  ####--- Block involving Correlation only ---####
  if(T){
    # IM_cor (testorder+1)*(N-1) x (testorder+1)*(N-1), SUM OVER TIME
    mHelp7 <- matrix(0,nrow=(testorder+1),ncol=(testorder+1))
    if (testorder==1){
      mHelp7 <- matrix(c(1,1/2,1/2,1/3),nrow=2,ncol=2)
    }
    if (testorder==2){
      mHelp7 <- matrix(c(1,1/2,1/3,1/2,1/3,1/4,1/3,1/4,1/5),nrow=3,ncol=3)
    }
    if (testorder==3){
      mHelp7 <- matrix(c(1,1/2,1/3,1/4,1/2,1/3,1/4,1/5,1/3,1/4,1/5,1/6,1/4,1/5,1/6,1/7),nrow=4,ncol=4)
    }
    IM_cor <- ( t(x_tau)%*%x_tau ) %x% (2*L2.inv[1:(N-1),1:(N-1)] + 2*L2.inv[N,N]*One_N_1.N_1) / Tobs # (testorder+1)*N x (testorder+1)*N
    #IM_cor <- (0.25*mHelp7) %x% (2*L2.inv[1:(N-1),1:(N-1)]  + 2*L2.inv[N,N]*One_N_1.N_1) # (testorder+1)*N x (testorder+1)*N
    # Check Size:
    if(NCOL(IM_cor) != (testorder+1)*(N-1)) stop("IM_cor is the wrong size")
    if(NROW(IM_cor) != (testorder+1)*(N-1)) stop("IM_cor is the wrong size")
  }


  # IM
  # IM <- rbind(cbind(IM_tv,IM_tv_garch,IM_tv_cor),cbind(t(IM_tv_garch),IM_garch,IM_garch_cor),cbind(t(IM_tv_cor),t(IM_garch_cor),IM_cor)) # the whole IM, not really needed...
  ####--- Block corresponding to correlations of the inverse of the IM matrix ---####

  IM <- rbind(cbind(IM_tv,IM_tv_garch,IM_tv_cor),cbind(t(IM_tv_garch),IM_garch,IM_garch_cor),cbind(t(IM_tv_cor),t(IM_garch_cor),IM_cor)) # the whole IM, not really needed...
  IM_inv <- solve(IM)

  ## The above is not necessarily very efficient - could use block inversion methods - TO DO LATER

  if (!is.matrix(IM_inv)) stop("LM Test: Can't invert the information matrix IM_inv")


  ####--- Block corresponding to the corr.parameters that are set to zero under null ---####
  SE_dim <- testorder*(N-1)
  block_start <- NCOL(IM_inv)-SE_dim+1
  block_end <- NCOL(IM_inv)
  IM_inv_SE <- IM_inv[(block_start:block_end),(block_start:block_end)]

  ##--- Return LM ---####
  LM <- (1/Tobs)*t(dlldrho_A)%*%IM_inv_SE%*%dlldrho_A


}  # End:  Parsim



myTest.CCCvSTCC.LM.new <- function(e,H0,H1,testorder=1) {

  ####--- Initialise ---####
  if(T){

    nGARCH <- H0$nGARCH
    nTV <- H0$nTV
    CCC <- H0$CCC
    STCC <- H1
    #NEW <- list()

    if (is.null(CCC)) stop("LM Test: need CCC as input (H0)")
    if (is.null(STCC$st)) stop("LM Test: need transition variable for the test as input (H1)")

    N <- NCOL(e)  # Number of series
    Tobs <- NROW(e)
    P <- CCC$P
    st <- STCC$st
    I <- diag(nrow = N,ncol = N) # identity matrix
    Pinv <- try(solve(P))
    if (!is.matrix(Pinv)) stop("LM Test: Can't invert the correlation matrix P")

    Num_garch_pars <- 0
    Num_tv_pars <- 0

    #Construct the K matrix: N^2 x N^2
    K <- NULL
    for (i in 1:N) {
      # block rows
      Krow <- NULL
      for (j in 1:N) {
        # block columns
        block <- matrix(0,N,N)
        block[j,i] <- 1
        Krow <- cbind(Krow,block)
      }
      K <- rbind(K,Krow)
    }

    #Construct the U matrix:Dimensions = N^2 x N*(N-1)/2
    U <- NULL
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        block <- matrix(0,N,N)
        block[i,j] <- block[j,i] <- 1
        Ucol <- as.vector(block)
        U <- cbind(U,Ucol)
      }
    }

  }

  # 2. Create the Information Matrix for the Correlation, and extract the SE corner of its inverse:

  ####--- Initialise IM Matrix variables ---####
  if(T){
    IM_cor <- matrix(0,(testorder+1)*N*(N-1)/2,(testorder+1)*N*(N-1)/2)
    g <- matrix(1,Tobs,N) # for now, no TV estimated, so these are 1, T x N
    h <- matrix(1,Tobs,N) # place for estimated GARCH variances, ones at this stage, T x N
    beta <- rep(0,N) # estimated beta coefficients from GARCH equations
    # Extract the number of TV & Garch pars:
    for (n in 1:N) {
      Num_tv_pars <- Num_tv_pars + nTV[[n]]@nr.pars
      g[,n] <- nTV[[n]]@g # univariate TV volatilities, T x N
      if(nGARCH[[n]]$type != garchtype$noGarch) {
        Num_garch_pars <- Num_garch_pars + nGARCH[[n]]@nr.pars
        h[,n] <- nGARCH[[n]]@h        # univariate GARCH volatilities, T x N
        beta[n] <- nGARCH[[n]]$Estimated$pars["beta",1]
      }
    }
    # TODO: Convert for() loops to lapply functions
    #Num_garch_pars <- lapply(nGARCH,sum())
    IM_garch <- matrix(0,Num_garch_pars,Num_garch_pars)
    IM_garch_cor <- matrix(0,Num_garch_pars,(testorder+1)*N*(N-1)/2)

    w <- e/sqrt(g) # returns standardised by g(t), T x N
    z <- w/sqrt(h) # returns standardised by GARCH volatilities and TV g(t), T x N
  }

  # partial derivatives of P w.r.t correlation parameters, T x 1 or T x 2
  # these are for alternative parameters, later add column of ones at the front

  ####--- Set v_rho ---####
  if(T){
    if (testorder==1) v_rho <- cbind(st)
    if (testorder==2) v_rho <- cbind(st,st^2)
    if (testorder==3) v_rho <- cbind(st,st^2,st^3)

    # score for rho
    zKRONz <- matrix(0,nrow=N^2,ncol=Tobs)
    for (t in 1:Tobs) {
      zKRONz[,t] <- t(z[t,,drop=FALSE] %x% z[t,,drop=FALSE]) # (N^2 x T), each col = "z_t kron z_t"
    }
    scaleFactor <- matrix(1,nrow=1,ncol=Tobs)
    dlldrho_A <- -0.5*t(U)%*%( vec(Pinv)%x%scaleFactor-(Pinv%x%Pinv)%*%(zKRONz) )%*%v_rho # N*(N-1)/2 x testorder
    dlldrho_A <- vec(dlldrho_A) # testorder*N*(N-1)/2 x 1, SUM OVER TIME
    v_rho <- cbind(1,v_rho) # T x 2 or T x 3, now add column of ones at the front

    # Check the size:
    if(NROW(dlldrho_A) != testorder*N*(N-1)/2) stop("dlldrho_A is the wrong size")

  }

  ####--- partial derivatives for GARCH and "x" matrices for GARCH parts ---####
  if(T){
    # partial derivatives of h1,h2,...,hN w.r.t garch_pars
    v_garch <- NULL
    h_scale <- NULL
    beta_scale <- NULL
    for (n in 1:N){
      if(nGARCH[[n]]$type==garchtype$noGarch ) {
        # Do nothing, move on to next series
      } else if(nGARCH[[n]]$type==garchtype$general) {
        # standard GARCH(1,1) case
        v_garch <- cbind(v_garch,c(0,rep(1,(Tobs-1))),c(0,w[1:(Tobs-1),n]^2),c(0,h[1:(Tobs-1),n])) # T x Num_garch_pars, each row = "1~w(i,t-1)^2~h(i,t-1)", i=1,...,N
        beta_scale <- c(beta_scale,(beta[n] %x% rep(1,3)))    # vector, length 3Num_garch_pars
        scaleFactor <- matrix(1,nrow=1,ncol=3)
        h_scale <- cbind(h_scale,h[,n,drop=FALSE] %x% scaleFactor)                     # T x Num_garch_pars
      } else if (nGARCH[[n]]$type > garchtype$general) {
        # GJR or other forms not done yet
        print("ARRGH!")
      }
    }

    if(Num_garch_pars > 0){
      dhdt <- myFilter(v_garch,beta_scale) # T x Num_garch_pars, each row = "(dh(i,t)/dtheta(i))' ", i=1,...,N
      # matrix containing the x_t's (garch part)
      x_garch <- -0.5*dhdt/h_scale  # T x Num_garch_pars, each row = "x_it'", i=1,...,N
      # Check the size:
      if(NCOL(dhdt) != Num_garch_pars) stop("dhdt is the wrong size")
      if(NCOL(x_garch) != Num_garch_pars) stop("x_garch is the wrong size")
    } else {
      x_garch <- matrix(nrow=Tobs,ncol=0)
    }

  }
  # NEW[[1]] <- IM_garch
  # NEW[[2]] <- IM_garch_cor
  # NEW[[3]] <- zKRONz
  # NEW[[4]] <- dlldrho_A
  # NEW[[5]] <- v_rho
  # NEW[[6]] <- dhdt
  # NEW[[7]] <- x_garch

  ####--- partial derivatives for TV and "x" matrices for TV parts ---####

  if(T){
    v_tv <- NULL
    h_scale <- matrix(NA,nrow = Tobs,ncol = 0)
    g_scale <- matrix(NA,nrow = Tobs,ncol = 0)
    beta_scale <- NULL
    dgdt_n <- NULL
    dgdt <- NULL
    dhdt <- NULL

    # partial derivatives of g1,g2,...,gN w.r.t tv_pars
    for (n in 1:N){
      dgdt_n <- dg_dt(nTV[[n]])  # T x #of TVpars in TV[n], includes d0's derivative only if d0 is free

      dgdt <- cbind(dgdt,dgdt_n) # T x #of TVpars
      scaleFactor <- matrix(1,nrow=1,ncol=NCOL(dgdt_n))
      g_scale <- cbind( g_scale,g[,n,drop=FALSE] %x% scaleFactor ) # cbinds g(nt) as many times as g(n) has tv parameters

      if (nGARCH[[n]]$type==garchtype$noGarch){
        # do nothing
      } else if (nGARCH[[n]]$type==garchtype$general) {
        # standard GARCH(1,1) case
        v_tv_n <- ((-nGARCH[[n]]$Estimated$pars["alpha",1]*c(0,1/g[1:(Tobs-1),n])*c(0,w[1:(Tobs-1),n]^2))%x% scaleFactor) *dgdt_n
        v_tv <- cbind(v_tv,v_tv_n)
        beta_scale_n <- (beta[n] %x% scaleFactor)
        beta_scale <- c(beta_scale,beta_scale_n)
        h_scale_n <- h[,n,drop=FALSE] %x% scaleFactor
        h_scale <- cbind(h_scale,h_scale_n)
      } else if(nGARCH[[n]]$type > garchtype$general) {
        # GJR or other type not done yet!
        stop("GARCH Type not supported by this Test yet")
      }
    } # End: for loop

    if(Num_tv_pars > 0 && Num_garch_pars > 0){
      dhdt <- myFilter(v_tv,beta_scale) # T x Num_tv_pars, each row = dh(i,t).dtvpar(i), i=1...N
      x_tv <- -0.5*dhdt/h_scale -0.5*dgdt/g_scale  # T x Num_tv_pars, each row = x_it, i=1...N
    } else if(Num_tv_pars > 0 && Num_garch_pars == 0){
      x_tv <- -0.5*dgdt/g_scale  # T x Num_tv_pars, each row = x_it, i=1...N
    } else {
      x_tv <- matrix(1,nrow=Tobs,ncol=0)
    }

    # Check Size:
    if(NCOL(x_tv) != Num_tv_pars) stop("x_tv is the wrong size")

  }

  #NEW[[8]] <- x_tv
  #saveRDS(NEW,"NEW.RDS")

  # Information matrix:
  ####--- Blocks involving GARCH only or GARCH and Correlation ---####
  # IM_garch (Num_garch_pars x Num_garch_pars), SUM OVER TIME
  One_31 <- matrix(1,3,1)
  One_33 <- matrix(1,3,3)
  I.P.Pinv <- I + P*Pinv
  I.P.Pinv_scale <- NULL
  scaleFactor <- NULL
  for (i in 1:N) {
    # i = row index
    I.P.Pinv_scale_row <- NULL
    if (nGARCH[[i]]$type == garchtype$noGarch){
      # do nothing
    } else {
      for (j in 1:N) {
        if(nGARCH[[j]]$type == garchtype$noGarch){
          # do nothing
        } else {
          # j = col index
          if(nGARCH[[i]]$type==garchtype$general && nGARCH[[j]]$type==garchtype$general) scaleFactor <- One_33
          else stop("Garch type is not supported")

          I.P.Pinv_scale_row <- cbind(I.P.Pinv_scale_row, (I.P.Pinv[i,j] %x% scaleFactor))
        }
      } # End: for "j" loop
    }
    I.P.Pinv_scale <- rbind(I.P.Pinv_scale,I.P.Pinv_scale_row)
  } # End: for "i" loop

  if(!is.null(I.P.Pinv_scale)){
    IM_garch <- ((t(x_garch)%*%x_garch) * I.P.Pinv_scale) / Tobs

    #NEW[[9]] <- IM_garch

    # Check Size:
    if(NCOL(IM_garch) != Num_garch_pars) stop("IM_garch is the wrong size")
    if(NROW(IM_garch) != Num_garch_pars) stop("IM_garch is the wrong size")

    # IM_garch_cor, 3N x (testorder+1)*N*(N-1)/2
    Mhelp1 <- matrix(0,nrow = 0,ncol=N*N)  # N x N^2
    for (n in 1:N) Mhelp1 <- rbind(Mhelp1, (Pinv[n,]%x%I[n,]+I[n,]%x%Pinv[n,])) # N x N^2

    Mhelp1 <- -0.5*(Mhelp1%*%U)
    Mhelp1_scale <- matrix(0,nrow = 0,ncol=N*(N-1)/2)  # (Num_Garch_Pars x N*(N-1)/2)

    for (n in 1:N) {
      if (nGARCH[[n]]$type==garchtype$general) Mhelp1_scale <- rbind(Mhelp1_scale, (Mhelp1[n, ,drop=FALSE] %x% One_31))
    }
    #NEW[[10]] <- Mhelp1_scale

    Mhelp2 <- t(t(v_rho)%*%x_garch)/Tobs # Num_garch_pars x 2 (or x3 if TestOrder=2), SUM OVER TIME
    # Check Size:
    if(NCOL(Mhelp2) != testorder+1) stop("Mhelp2 is the wrong size - COLs")
    if(NROW(Mhelp2) != Num_garch_pars) stop("Mhelp2 is the wrong size - ROWs")

    IM_garch_cor <- matrix(0,nrow = Num_garch_pars,ncol=0)  # (Num_Garch_Pars x (testorder+1)*N*(N-1)/2)
    for (i in 1:NCOL(v_rho)){
      IM_garch_cor <- cbind(IM_garch_cor,((Mhelp2[,i,drop=FALSE]%x%t(rep(1,(N*(N-1))/2))) * Mhelp1_scale))  # (Num_Garch_Pars x (testorder+1)*N*(N-1)/2),  SUM OVER TIME
    }
    # Check Size:
    if(NCOL(IM_garch_cor) != (testorder+1)*N*(N-1)/2) stop("IM_garch_cor is the wrong size")
  } # End: if(!is.null(I.P.Pinv_scale))  => Only process when we have Garch

  #NEW[[11]] <- IM_garch_cor

  ####--- Blocks involving TV only or TV and Correlation ---####
  #if (nTV$Type != TRshape$none){
  IM_tv <- NULL
  I.P.Pinv <- I + P*Pinv
  I.P.Pinv_scale <- NULL
  scaleFactor <- NULL
  for (i in 1:N) {
    # i = row index
    ii <- nTV[[i]]@nr.pars
    if (ii==0){
      # do nothing
    } else {
      I.P.Pinv_scale_row <- NULL
      for (j in 1:N) {
        # j = col index
        jj <- nTV[[j]]@nr.pars
        if (jj==0){
          # do nothing
        } else {
          scaleFactor <- matrix(1,ii,jj)
          I.P.Pinv_scale_row <- cbind(I.P.Pinv_scale_row, (I.P.Pinv[i,j] %x% scaleFactor))
        } # End: if (jj==0)
      } # End: for (j in 1:N)
      I.P.Pinv_scale <- rbind(I.P.Pinv_scale,I.P.Pinv_scale_row)
    } # End: if (ii==0)
  } # End: for (i in 1:N)
  IM_tv <- ((t(x_tv)%*%x_tv) * I.P.Pinv_scale) / Tobs  # sum(num_tv_pars) x sum(num_tv_pars)

  #NEW[[12]] <- IM_tv

  # IM_tv_cor, (#of tv pars in n,n=1...N)N x (testorder+1)*N*(N-1)/2
  Mhelp3 <- matrix(0,nrow = 0,ncol = N*N)
  for (n in 1:N){
    Mhelp3 <- rbind(Mhelp3, (Pinv[n,,drop=FALSE]%x%I[n,,drop=FALSE]+I[n,,drop=FALSE]%x%Pinv[n,,drop=FALSE])) # N x N^2
  }
  # Check Size:
  if(NCOL(Mhelp3) != N^2) stop("Mhelp3 is the wrong size - COLs")
  if(NROW(Mhelp3) != N) stop("Mhelp3 is the wrong size - ROWs")

  Mhelp3 <- -0.5*(Mhelp3%*%U)  # N x N*(N-1)/2
  # Check Size:
  if(NCOL(Mhelp3) != N*(N-1)/2) stop("Mhelp3 is the wrong size - COLs")
  if(NROW(Mhelp3) != N) stop("Mhelp3 is the wrong size - ROWs")

  Mhelp3_scale <- matrix(0,nrow=0,ncol=(N*(N-1)/2))  # (Total Num_tv_Pars x N*(N-1)/2)
  for (n in 1:N) {
    ii <- nTV[[i]]@nr.pars
    if (ii==0){
      # do nothing
    } else {
      #if (nGARCH[[n]]$type==garchtype$noGarch) ii <- ii+1   #add delta0 to the count
      scaleFactor <- matrix(1,nrow=ii,ncol=1)
      Mhelp3_scale <- rbind(Mhelp3_scale, Mhelp3[n,,drop=FALSE] %x% scaleFactor)
    } # End: if (ii==0)
  } # End: for (n in 1:N)

  # Check Size:
  if(NCOL(Mhelp3_scale) != N*(N-1)/2) stop("Mhelp3_scale is the wrong size - COLs")
  if(NROW(Mhelp3_scale) != Num_tv_pars) stop("Mhelp3_scale is the wrong size - ROWs")

  Mhelp4 <- t(t(v_rho)%*%x_tv)/Tobs # Num_tv_pars x 2 (or x3 if TestOrder=2), SUM OVER TIME
  if(NCOL(Mhelp4) != testorder+1) stop("Mhelp4 is the wrong size")
  if(NROW(Mhelp4) != Num_tv_pars) stop("Mhelp4 is the wrong size")

  IM_tv_cor <- matrix(0,nrow=Num_tv_pars,ncol=0)
  scaleFactor <- matrix(1,nrow=1,ncol=N*(N-1)/2)
  for (i in 1:NCOL(v_rho)){
    IM_tv_cor <- cbind(IM_tv_cor,(Mhelp4[,i,drop=FALSE]%x%scaleFactor) * Mhelp3_scale)  # (Num_tv_Pars x (testorder+1)*N*(N-1)/2),  SUM OVER TIME
  }

  # Check Size:
  if(NCOL(IM_tv) != Num_tv_pars) stop("IM_tv is the wrong size")
  if(NCOL(IM_tv_cor) != (testorder+1)*N*(N-1)/2) stop("IM_tv_cor is the wrong size")
  if(NROW(IM_tv_cor) != Num_tv_pars) stop("IM_tv_cor is the wrong size")

  #NEW[[13]] <- IM_tv_cor
  #saveRDS(NEW,"NEW.RDS")

  ####--- Blocks involving TV and GARCH ---####
  #if (nGARCH$Type != garchtype$noGarch && nTV$Type != TRshape$none){
  I.P.Pinv <- I + P*Pinv  # N x N
  I.P.Pinv_scale <- NULL
  for (i in 1:N){
    # i = row index
    ii <- nTV[[i]]@nr.pars
    if (ii==0){
      # do nothing
    } else {
      I.P.Pinv_scale_row <- NULL
      for (j in 1:N){
        # j = col index
        jj <- nGARCH[[j]]@nr.pars
        if (jj==0){
          # do nothing
        } else {
          scaleFactor <- matrix(1,ii,jj)
          I.P.Pinv_scale_row <- cbind(I.P.Pinv_scale_row, (I.P.Pinv[i,j] %x% scaleFactor))
        } # End if (jj==0)
      } # End for (j in 1:N)
      I.P.Pinv_scale <- rbind(I.P.Pinv_scale,I.P.Pinv_scale_row)  # (total of tvpars) x (total of garchpars)
    } # End if (ii==0)
  }
  IM_tv_garch <- NULL
  IM_tv_garch <- (t(x_tv)%*%x_garch)*(I.P.Pinv_scale)/Tobs  #  (total # of tvpars)x(total # of garchpars)

  # Check Size:
  if(NCOL(IM_tv_garch) != Num_garch_pars) stop("IM_tv_garch is the wrong size")
  if(NROW(IM_tv_garch) != Num_tv_pars) stop("IM_tv_garch is the wrong size")

  #NEW[[14]] <- I.P.Pinv_scale
  #NEW[[15]] <- IM_tv_garch
  #saveRDS(NEW,"NEW.RDS")

  #}

  ####--- Block involving Correlation only ---####
  if(T){
    # IM_cor (testorder+1)*N*(N-1)/2, SUM OVER TIME
    Mhelp5 <- t(U)%*%(Pinv%x%Pinv+(Pinv%x%I)%*%K%*%(Pinv%x%I))%*%U # N*(N-1)/2 x N*(N-1)/2
    Mhelp6 <- t(v_rho)%*%v_rho/Tobs # 2x2 or 3x3, SUM OVER TIME
    IM_cor <- 0.25*(Mhelp6%x%Mhelp5) # (testorder+1)*N*(N-1)/2 x (testorder+1)*N*(N-1)/2
    # Check Size:
    if(NCOL(Mhelp5) != N*(N-1)/2) stop("Mhelp5 is the wrong size")
    if(NROW(Mhelp5) != N*(N-1)/2) stop("Mhelp5 is the wrong size")
    #if(NCOL(Mhelp6) != NCOL(v_rho)) stop("Mhelp6 is the wrong size")
    #if(NROW(Mhelp6) != NROW(v_rho)) stop("Mhelp6 is the wrong size")
    if(NCOL(IM_cor) != (testorder+1)*N*(N-1)/2) stop("IM_cor is the wrong size")
    if(NROW(IM_cor) != (testorder+1)*N*(N-1)/2) stop("IM_cor is the wrong size")

  }

  #NEW[[16]] <- IM_cor
  #saveRDS(NEW,"NEW.RDS")

  # IM
  # IM <- rbind(cbind(IM_tv,IM_tv_garch,IM_tv_cor),cbind(t(IM_tv_garch),IM_garch,IM_garch_cor),cbind(t(IM_tv_cor),t(IM_garch_cor),IM_cor)) # the whole IM, not really needed...
  ####--- Block corresponding to correlations of the inverse of the IM matrix ---####

  IM <- rbind(cbind(IM_tv,IM_tv_garch,IM_tv_cor),cbind(t(IM_tv_garch),IM_garch,IM_garch_cor),cbind(t(IM_tv_cor),t(IM_garch_cor),IM_cor)) # the whole IM, not really needed...
  IM_inv <- solve(IM)

  ## The above is not necessarily very efficient - could use block inversion methods - TO DO LATER

  if (!is.matrix(IM_inv)) stop("LM Test: Can't invert the information matrix IM_inv")

  ####--- Block corresponding to the corr.parameters that are set to zero under null ---####
  SE_dim <- testorder*(N*(N-1))/2
  block_start <- NCOL(IM_inv)-SE_dim+1
  block_end <- NCOL(IM_inv)
  IM_inv_SE <- IM_inv[(block_start:block_end),(block_start:block_end)]



  # Clean up:
  #rm(e,w,z,g,h,H0,H1,nGARCH,nTV,STCC,K,U,IM_tv,IM_tv_garch,IM_tv_cor,IM_garch,IM_garch_cor,IM_cor,IM,IM_inv)
  #rm(zKRONz,scaleFactor,v_tv,h_scale,g_scale,beta_scale,dgdt_n,dgdt,dhdt,v_garch,h_scale,beta_scale)
  #rm(Mhelp1,Mhelp2,Mhelp3,Mhelp4,Mhelp5,Mhelp6)
  #rm(I.P.Pinv,I.P.Pinv_scale)

  ##--- Return LM ---####
  LM <- (1/Tobs)*t(dlldrho_A)%*%IM_inv_SE%*%dlldrho_A

}  # End: myTest.CCCvSTCC.LM.new()



## ===== Test Sub Functions =====####


setGeneric(name=".x_tau",
           valueClass = "list",
           signature = c("z","H0","testOrder"),
           def = function(z,H0,testOrder){
             if (testOrder==1) x_tau <- (-0.5)*cbind(st)
             if (testOrder==2) x_tau <- (-0.5)*cbind(st,st^2)
             if (testOrder==3) x_tau <- (-0.5)*cbind(st,st^2,st^3)

             N <- H0@N
             P <- H0$Estimated$P
             tmp <- eigen(P)
             Q <- tmp$vectors
             L.inv <- diag(tmp$values^(-1)) # matrix

             One_N_1.1 <- matrix(1,nrow=(N-1),ncol=1)
             I <- diag(nrow = N,ncol = N) # NxN Identity matrix
             I.P.Pinv <- I + P*solve(P)
             I.P.Pinv_scale <- NULL
             scaleFactor <- NULL

             # score for rho
             qxq <- matrix(0,nrow=N,ncol=(N^2))  # each row = qi' %x% qi', 1 x N^2, N rows
             for (n in 1:N){
               qxq[n,] <- t(Q[,n,drop=FALSE]%x%Q[,n,drop=FALSE])
             }
             mHelp0 <- matrix(0,nrow=Tobs,ncol=(N-1))
             for (t in 1:H0@Tobs){
               vec_zz_t <- vec(t(z[t,,drop=FALSE]) %*% (z[t,,drop=FALSE]))  # N^2 x 1
               mHelp0[t,] <-  L.inv[1:(N-1),1:(N-1)] %*% (One_N_1.1 - L.inv[1:(N-1),1:(N-1)] %*% qxq[1:(N-1),,drop=FALSE] %*% vec_zz_t) - One_N_1.1 %*% L.inv[N,N] %*% (1-L.inv[N,N] * (qxq[N,,drop=FALSE] %*% vec_zz_t))  # 1 x N-1
             }

             dlldrho_A <- t(mHelp0) %*% x_tau   # (N-1)xT %*% (T x testorder) = N-1 x testorder, SUM OVER TIME
             dlldrho_A <- vec(dlldrho_A)  # testorder*(N-1) x 1, SUM OVER TIME
             x_tau <- cbind(rep(-0.5,Tobs),x_tau) # T x 2 or T x 3, now add column of ones at the front (Note:xtau includes -0.5)

             rtn <- list()
             rtn$x_tau <- x_tau
             rtn$dlldrho_A <- dlldrho_A
             return(rtn)

           }
)


setGeneric(name=".x_garch",
           valueClass = "matrix",
           signature = c("w","H0"),
           def = function(w,H0){

             # partial derivatives of h1,h2,...,hN w.r.t garch_pars
             v_garch <- NULL
             h_scale <- NULL
             beta_scale <- NULL
             # TODO: Do we need to ensure the matrix is padded to 4 columns?
             x_garch <- matrix(NA,H0@Tobs,4)  # Initialise return matrix to be able to contain largest Garch Specification (GJR - 4 params)

             # Loop over every series:
             for (n in 1:H0@N){
               h <- as.matrix(H0[[n]]$garchObj@h)
               beta <- as.matrix(H0[[n]]$garchObj$Estimated$pars["beta",1])

               if(H0[[n]]$garchObj$type==garchtype$noGarch ) {
                 # Do nothing, move on to next series
                 x_garch <- matrix(NA,nrow=H0@Tobs,ncol=0)

               } else if(H0[[n]]$garchObj$type==garchtype$general) {
                 # standard GARCH(1,1) case
                 v_garch <- cbind(v_garch,c(0,rep(1,(H0@Tobs-1))),c(0,w[1:(H0@Tobs-1),n]^2),c(0,h[1:(H0@Tobs-1),1])) # T x Num_garch_pars, each row = "1~w(i,t-1)^2~h(i,t-1)", i=1,...,N
                 beta_scale <- beta_scale <- c(beta_scale,(beta %x% c(1,1,1)))    # vector, length 3
                 scaleFactor <- matrix(1,nrow=1,ncol=3)
                 h_scale <- cbind(h_scale,h %x% scaleFactor)    # T x Num_garch_pars
                 #
                 dhdt <- .ar1.Filter(v_garch,beta_scale) # T x Num_garch_pars, each row = "(dh(i,t)/dtheta(i))' ", i=1,...,N
                 x_garch <- -0.5*dhdt/h_scale # T x Num_garch_pars, each row = "x_it'", i=1,...,N

               } else if (H0[[n]]$garchObj$type == garchtype$gjr) {
                 v_garch <- cbind(v_garch,c(0,rep(1,(H0@Tobs-1))),c(0,w[1:(H0@Tobs-1),n]^2),c(0,h[1:(H0@Tobs-1),1])) # T x Num_garch_pars, each row = "1~w(i,t-1)^2~h(i,t-1)", i=1,...,N
                 beta_scale <- beta_scale <- c(beta_scale,(beta %x% rc(1,1,1,1)))    # vector, length 4
                 scaleFactor <- matrix(1,nrow=1,ncol=4)
                 h_scale <- cbind(h_scale,h %x% scaleFactor)    # T x Num_garch_pars
                 #
                 dhdt <- .ar1.Filter(v_garch,beta_scale) # T x Num_garch_pars, each row = "(dh(i,t)/dtheta(i))' ", i=1,...,N
                 x_garch <- -0.5*dhdt/h_scale # T x Num_garch_pars, each row = "x_it'", i=1,...,N
                 warning("Annastiina - please confirm this matrix is correct for GJR garch!")
               }

             }
             return(x_garch)

           }
)

setGeneric(name=".x_tv",
           valueClass = "matrix",
           signature = c("w","H0"),
           def = function(w,H0){

             v_tv <- NULL
             h_scale <- NULL
             g_scale <- NULL
             beta_scale <- NULL
             dgdt <- NULL
             dhdt <- NULL

             # partial derivatives of g1,g2,...,gN w.r.t tv_pars
             for (n in 1:H0@N){

               if(H0[[n]]$tvObj@nr.pars == 0){
                 x_tv <- matrix(1,nrow=H0@Tobs,ncol=0)
                 break  #Goto next 'n' in for..loop
               }

               ## All code below is assured that there are tv parameters
               g <- as.matrix(H0[[n]]$tvObj@g)
               h <- as.matrix(H0[[n]]$garchObj@h)
               beta <- as.matrix(H0[[n]]$garchObj$Estimated$pars["beta",1])

               dgdt_n <- .dg_dt(H0[[n]]$tvObj)  # T x #of TVpars in TV[n], includes d0's derivative if it is a free param
               dgdt <- cbind(dgdt,dgdt_n)
               scaleFactor <- matrix(1,nrow=1,ncol=NCOL(dgdt_n))
               g_scale <- cbind( g_scale,g %x% scaleFactor ) # cbinds g(nt) as many times as g(n) has tv parameters

               if (H0[[n]]$garchObj$type==garchtype$noGarch){
                 x_tv <- -0.5*dgdt/g_scale  # T x Num_tv_pars, each row = x_it, i=1...N

               } else if (H0[[n]]$garchObj$type==garchtype$general) {
                 # standard GARCH(1,1) case
                 v_tv_n <- ((-H0[[n]]$garchObj$Estimated$pars["alpha",1]*c(0,1/g[1:(Tobs-1),1])*c(0,w[1:(Tobs-1),n]^2))%x% scaleFactor) * dgdt_n
                 v_tv <- cbind(v_tv,v_tv_n)
                 beta_scale_n <- (beta %x% scaleFactor)
                 beta_scale <- c(beta_scale,beta_scale_n)
                 h_scale_n <- h %x% scaleFactor
                 h_scale <- cbind(h_scale,h_scale_n)
                 dhdt <- .ar1.Filter(v_tv,beta_scale) # T x Num_tv_pars, each row = dh(i,t).dtvpar(i), i=1...N
                 x_tv <- -0.5*dhdt/h_scale -0.5*dgdt/g_scale  # T x Num_tv_pars, each row = x_it, i=1...N

               } else if(H0[[n]]$garchObj$type == garchtype$gjr) {
                 v_tv_n <- ((-H0[[n]]$garchObj$Estimated$pars["alpha",1]*c(0,1/g[1:(Tobs-1),1])*c(0,w[1:(Tobs-1),n]^2))%x% scaleFactor) * dgdt_n
                 v_tv <- cbind(v_tv,v_tv_n)
                 beta_scale_n <- (beta %x% scaleFactor)
                 beta_scale <- c(beta_scale,beta_scale_n)
                 h_scale_n <- h %x% scaleFactor
                 h_scale <- cbind(h_scale,h_scale_n)
                 dhdt <- .ar1.Filter(v_tv,beta_scale) # T x Num_tv_pars, each row = dh(i,t).dtvpar(i), i=1...N
                 x_tv <- -0.5*dhdt/h_scale -0.5*dgdt/g_scale  # T x Num_tv_pars, each row = x_it, i=1...N
               }
             } # End: for loop

             return(x_tv)

           }
)

setGeneric(name=".im_garch",
           valueClass = "list",
           signature = c("H0","x_garch","x_tau"),
           def = function(H0,x_garch,x_tau){
             # IM_garch (Num_garch_pars x Num_garch_pars), SUM OVER TIME
             N <- H0@N
             P <- H0$Estimated$P
             tmp <- eigen(P)
             Q <- tmp$vectors
             L.inv <- diag(tmp$values^(-1)) # matrix
             L2.inv <- diag(tmp$values^(-2)) # matrix
             L <- diag(tmp$values) # matrix

             One_31 <- matrix(1,3,1)
             One_33 <- matrix(1,3,3)
             I <- diag(nrow = N,ncol = N) # NxN Identity matrix
             I.P.Pinv <- I + P*solve(P)
             I.P.Pinv_scale <- NULL
             scaleFactor <- NULL

             for (i in 1:N) {
               # i = row index
               I.P.Pinv_scale_row <- NULL
               if (H0[[i]]$garchObj$type == garchtype$noGarch){
                 # do nothing
               } else {
                 for (j in 1:N) {
                   if(H0[[j]]$garchObj$type == garchtype$noGarch){
                     # do nothing
                   } else {
                     # j = col index
                     if(H0[[i]]$garchObj$type==garchtype$general && H0[[j]]$garchObj$type==garchtype$general) scaleFactor <- One_33
                     else stop("Garch type is not supported")
                     I.P.Pinv_scale_row <- cbind(I.P.Pinv_scale_row, (I.P.Pinv[i,j] %x% scaleFactor))
                   }
                 } # End: for "j" loop
               }
               I.P.Pinv_scale <- rbind(I.P.Pinv_scale,I.P.Pinv_scale_row)
             } # End: for "i" loop


             if(!is.null(I.P.Pinv_scale)){
               IM_garch <- ((t(x_garch)%*%x_garch) * I.P.Pinv_scale) / H0@Tobs
               # IM_garch_cor, will be 3N x (testorder+1)*(N-1)
               mHelp1 <- matrix(0,nrow=0,ncol=(N - 1) )

               for (n in 1:H0@N) {
                 ## TODO: What about GJR??
                 if (H0[[n]]$garchObj$type == garchtype$general) {
                   scaleFactor <- matrix(1,nrow=3,ncol=1)
                   mHelp1 <- rbind(mHelp1, scaleFactor %x% (2*Q[n,1:(N-1),drop=FALSE]^2 %*% L.inv[1:(N-1),1:(N-1)] - 2*Q[n,N]^2 * L.inv[N,N] * One_1.N_1))
                 }
               }
               mHelp2 <- t(x_garch) %*% x_tau
               mHelp3 <- matrix(0,nrow=NROW(mHelp2),ncol=(N-1)*NCOL(mHelp2))
               for (i in 1:NROW(mHelp2)){
                 mHelp3[i,] <- mHelp2[i,,drop=FALSE] %x% mHelp1[i,,drop=FALSE]
               }
               IM_garch_cor <- mHelp3/H0@Tobs # (Num_Garch_Pars x (testorder+1)*N), SUM OVER TIME

             } # End: if(!is.null(I.P.Pinv_scale))  => Only process when we have 'some' Garch

             rtn <- list()
             rtn$IM_garch <- IM_garch
             rtn$IM_garch_cor <- IM_garch_cor
             return(rtn)

           }
)

setGeneric(name=".im_tv",
           valueClass = "list",
           signature = c("H0","x_tv","x_tau"),
           def = function(H0,x_tv,x_tau){

             N <- H0@N
             P <- H0$Estimated$P
             tmp <- eigen(P)
             Q <- tmp$vectors
             L.inv <- diag(tmp$values^(-1)) # matrix
             L2.inv <- diag(tmp$values^(-2)) # matrix
             L <- diag(tmp$values) # matrix

             I <- diag(nrow = N,ncol = N) # NxN Identity matrix
             I.P.Pinv <- I + P * solve(P)
             I.P.Pinv_scale <- NULL
             scaleFactor <- NULL

             for (i in 1:N) {
               # i = row index
               if(H0[[i]]$tvObj@nr.pars > 0) {
                I.P.Pinv_scale_row <- NULL

                for (j in 1:N) {
                   # j = col index
                   if (H0[[j]]$tvObj@nr.pars > 0){
                     scaleFactor <- matrix(1,H0[[i]]$tvObj@nr.pars,H0[[j]]$tvObj@nr.pars)
                     I.P.Pinv_scale_row <- cbind(I.P.Pinv_scale_row, (I.P.Pinv[i,j] %x% scaleFactor))
                   } # End: if (jj==0)
                 } # End: for (j in 1:N)
                 I.P.Pinv_scale <- rbind(I.P.Pinv_scale,I.P.Pinv_scale_row)

               } # End: if (ii==0)
             } # End: for (i in 1:N)

             IM_tv <- ((t(x_tv) %*% x_tv) * I.P.Pinv_scale) / H0@Tobs  # sum(num_tv_pars) x sum(num_tv_pars)


             # IM_tv_cor, (total #of tv pars in model) x (testorder+1)*(N-1)
             mHelp4 <- matrix(0,nrow=0,ncol=(N-1) )
             for (n in 1:N) {
               if (H0[[n]]$tvObj@nr.pars > 0){
                 scaleFactor <- matrix(1,nrow=H0[[n]]$tvObj@nr.pars,ncol=1)
                 mHelp4 <- rbind(mHelp4, scaleFactor %x% (2*Q[n,1:(N-1),drop=FALSE]^2 %*% L.inv[1:(N-1),1:(N-1)] - 2*Q[n,N]^2 * L.inv[N,N] * One_1.N_1) )
               } # End: if (ii==0)
             } # End: for (n in 1:N)
             mHelp5 <- t(x_tv) %*% x_tau
             mHelp6 <- matrix(0,nrow=NROW(mHelp5),ncol=(N-1)*NCOL(mHelp5))
             if (NROW(mHelp6)>0){
               for (i in 1:NROW(mHelp6)){
                 mHelp6[i,] <- mHelp5[i,,drop=FALSE] %x% mHelp4[i,,drop=FALSE]
               }
             }

             IM_tv_cor <- mHelp6 / H0@Tobs  # (Num_tv_Pars x (testorder+1)*N),  SUM OVER TIME

             rtn <- list()
             rtn$IM_tv <- IM_tv
             rtn$IM_tv_cor <- IM_tv_cor
             return(rtn)


           }
)

setGeneric(name=".im_tv_garch",
           valueClass = "matrix",
           signature = c("H0","x_tv","x_garch"),
           def = function(H0,x_tv,x_garch){

             N <- H0@N
             P <- H0$Estimated$P
             I <- diag(nrow = N,ncol = N) # NxN Identity matrix
             I.P.Pinv <- I + P * solve(P)
             I.P.Pinv_scale <- NULL

             for (i in 1:N){
               # i = row index
               if (H0[[i]]$tvObj@nr.pars > 0){
                 I.P.Pinv_scale_row <- NULL
                 for (j in 1:N){
                   # j = col index
                   if (H0[[j]]$garchObj@nr.pars > 0){
                     scaleFactor <- matrix(1,H0[[i]]$tvObj@nr.pars,H0[[j]]$garchObj@nr.pars)
                     I.P.Pinv_scale_row <- cbind(I.P.Pinv_scale_row, (I.P.Pinv[i,j] %x% scaleFactor))
                   } # End if (jj==0)
                 } # End for (j in 1:N)
                 I.P.Pinv_scale <- rbind(I.P.Pinv_scale,I.P.Pinv_scale_row)  # (total of tvpars) x (total of garchpars)
               } # End if (ii==0)
             }
             #IM_tv_garch <- NULL
             IM_tv_garch <- (t(x_tv) %*% x_garch)*I.P.Pinv_scale/H0@Tobs  #  (total # of tvpars)x(total # of garchpars)

             return(IM_tv_garch)

           }
)

setGeneric(name=".im_cor",
           valueClass = "matrix",
           signature = c("H0","x_tau","testOrder"),
           def = function(H0,x_tau,testOrder){

             # IM_cor (testorder+1)*(N-1) x (testorder+1)*(N-1), SUM OVER TIME
             N <- H0@N
             P <- H0$Estimated$P
             tmp <- eigen(P)
             L2.inv <- diag(tmp$values^(-2)) # matrix
             One_N_1.N_1 <- matrix(1,nrow=(N-1),ncol=(N-1))

             if (testorder==1){
               mHelp7 <- matrix(c(1,1/2,1/2,1/3),nrow=2,ncol=2)
             }
             if (testorder==2){
               mHelp7 <- matrix(c(1,1/2,1/3,1/2,1/3,1/4,1/3,1/4,1/5),nrow=3,ncol=3)
             }
             if (testorder==3){
               mHelp7 <- matrix(c(1,1/2,1/3,1/4,1/2,1/3,1/4,1/5,1/3,1/4,1/5,1/6,1/4,1/5,1/6,1/7),nrow=4,ncol=4)
             }

             IM_cor <- ( t(x_tau)%*%x_tau ) %x% (2*L2.inv[1:(N-1),1:(N-1)] + 2*L2.inv[N,N]*One_N_1.N_1) / H0@Tobs # (testorder+1)*N x (testorder+1)*N

             return(IM_cor)

           }
)

setGeneric(name=".LM",
           valueClass = "matrix",
           signature = c("H0","IM_list","dlldrho_A","testOrder"),
           def = function(H0,IM_list,dlldrho_A,testOrder){

             N <- H0@N
             P <- H0$Estimated$P
             Pinv <- solve(P)
             I <- diag(nrow = N,ncol = N) # NxN Identity matrix

             IM_tv <- IM_list$IM_tv
             IM_tv_cor <- IM_list$IM_tv_cor
             IM_garch <- IM_list$IM_garch
             IM_garch_cor <- IM_list$IM_garch_cor
             IM_tv_garch <- IM_list$IM_tv_garch
             IM_cor <- IM_list$IM_cor

             IM <- rbind(cbind(IM_tv,IM_tv_garch,IM_tv_cor),cbind(t(IM_tv_garch),IM_garch,IM_garch_cor),cbind(t(IM_tv_cor),t(IM_garch_cor),IM_cor)) # the whole IM, not really needed...
             IM_inv <- solve(IM)
             ## The above is not necessarily very efficient - could use block inversion methods - TO DO LATER

             if (!is.matrix(IM_inv)) stop("LM Test: Can't invert the information matrix IM_inv")

             ####--- Block corresponding to the corr.parameters that are set to zero under null ---####
             block_start <- NCOL(IM_inv) - (testOrder*(N-1) + 1)
             block_end <- NCOL(IM_inv)
             IM_inv_SE <- IM_inv[(block_start:block_end),(block_start:block_end)]

             ##--- Return LM ---####
             return( (1/H0@Tobs)*t(dlldrho_A) %*% IM_inv_SE %*% dlldrho_A )

           }
)

