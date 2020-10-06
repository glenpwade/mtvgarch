## --- ccc_class Definition --- ####

ccc <- setClass(Class = "ccc_class",
                  slots = c(Tobs="integer",N="integer"),
                  contains = c("namedList")
)

## --- Initialise --- ####
setMethod("initialize","ccc_class",
          function(.Object,...){
            .Object <- callNextMethod(.Object,...)

            # Default initial values
            .Object@N <- as.integer(0)
            .Object@Tobs <- as.integer(0)

            # Return:
            .Object
          })

## -- Constructor:ccc -- ####
setGeneric(name="ccc",
           valueClass = "ccc_class",
           signature = c("mtvgarchObj"),
           def = function(mtvgarchObj){
             this <- new("ccc_class")

             ## -- Do validation checks -- ##
             objType <- class(mtvgarchObj)
             if(objType[1] != "mtvgarch_class"){
               warning("a valid instance of the mtvgarch_class is required to create a ccc model")
               return(this)
             }
             # End validation

             # Add the Estimated components from the mtvgarch
             this$mtvgarch <- list()
             for(n in 1:mtvgarchObj@N){
               this$mtvgarch[[n]] <- list()
               this$mtvgarch[[n]]$tv <- mtvgarchObj[[n]]$Estimated$tv
               this$mtvgarch[[n]]$garch <- mtvgarchObj[[n]]$Estimated$garch
             }
             names(this$mtvgarch) <- names(mtvgarchObj)

             # Set Default Values:
             this@N <- mtvgarchObj@N
             this@Tobs <- mtvgarchObj@Tobs

             return(this)
           }
)

## -- .loglik.ccc -- ####
setGeneric(name=".loglik.ccc",
           valueClass = "numeric",
           signature = c("optimpars","z","cccObj"),
           def = function(optimpars,z,cccObj){

             err_output <- -1e10

             vP <- optimpars
             mP <- .unVecL(vP)
             eig <- eigen(mP,symmetric=TRUE,only.values = TRUE)
             if (min(eig$values) <= 0) return(err_output)

             # - - - calc loglik-value
             Tobs <- NROW(z)
             llt <- NULL
             mPinv <- solve(mP)

             for(t in seq(1,Tobs)) llt[t] <- -0.5*log(det(mP)) - 0.5*(t(z[t,]) %*% mPinv %*% z[t,])

             return(sum(llt))

           }
)

## --- estimateCCC --- ####
setGeneric(name="estimateCCC",
           valueClass = "ccc_class",
           signature = c("z","cccObj","estimationCtrl"),
           def = function(z,cccObj,estimationCtrl){
             this <- cccObj

             this$Estimated <- list()

             this$Estimated$P <- cor(z)
             optimpars <- this$Estimated$P[lower.tri(this$Estimated$P)]
             this$Estimated$value <- .loglik.ccc(optimpars,z,this)

             return(this)
           }
)



##====  TESTS  ====####

##===  test.CCCParsim ===####
setGeneric(name="test.CCCParsim",
           valueClass = "numeric",
           signature = c("e","H0","H1","testOrder"),
           def = function(e,H0,H1,testOrder){

             # Validation
             objType <- class(H0)
             if(objType[1] != "ccc_class"){
               warning("This test requires a valid instance of an estimated ccc model as the null (H0)")
               return(0)
             }

             # Get the common variables:
             g <- matrix(1,H0@Tobs,H0@N)
             h <- matrix(1,H0@Tobs,H0@N)
             beta <- matrix(1,1,H0@N)
             for (n in 1:H0@N) {
               g[,n] <- H0$mtvgarch[[n]]$tv@g
               h[,n] <- H0$mtvgarch[[n]]$garch@h
               beta[1,n] <- H0$mtvgarch[[n]]$garch$Estimated$pars["beta",1]
             }
             w <- e/sqrt(g)
             z <- w/sqrt(h)

             # Get x_garch
             x_garch <- .x_garch(w,H0,h,beta)

             # Get x_tv
             x_tv <- .x_tv(z,H0,g,h,beta)

             # Get x_tau, dlldrho_A
             rtn <- .x_tau(z,H0,H1,testOrder)
             x_tau <- rtn$x_tau
             dlldrho_A <- rtn$dlldrho_A

             # Get im_garch
             im_garch <- .im_garch(H0,x_garch)

             # Get im_garch_cor
             im_garch_cor <- .im_garch_cor_parsim(H0,x_garch,x_tau)

             # Get im_tv,
             im_tv <- .im_tv(H0,x_tv)

             # Get im_tv_cor
             im_tv_cor <- .im_tv_cor_parsim(H0,x_tv,x_tau)

             # Get im_tv_garch
             Pt <- H0$Estimated$P
             im_tv_garch <- .im_tv_garch(H0,Pt,x_tv,x_garch)

             # Get im_cor
             im_cor <- .im_cor_parsim(H0,x_tau,testOrder)

             # Get LM using all InfoMatrix blocks
             IM_list <- list()
             IM_list$IM_tv <- im_tv
             IM_list$IM_tv_cor <- im_tv_cor
             IM_list$IM_garch <- im_garch
             IM_list$IM_garch_cor <- im_garch_cor
             IM_list$IM_tv_garch <- im_tv_garch
             IM_list$IM_cor <- im_cor

             LM <- .LM(H0,IM_list,dlldrho_A,testOrder)

             return(LM)

           }
)

##===  test.CCCvSTCC1 ===####
setGeneric(name="test.CCCvSTCC1",
           valueClass = "numeric",
           signature = c("e","H0","H1","testOrder"),
           def = function(e,H0,H1,testOrder){

             # Validation
             objType <- class(H0)
             if(objType[1] != "ccc_class"){
               warning("This test requires a valid instance of an estimated ccc model as the null (H0)")
               return(0)
             }

             # Get the common variables:
             g <- matrix(1,H0@Tobs,H0@N)
             h <- matrix(1,H0@Tobs,H0@N)
             beta <- matrix(1,1,H0@N)
             for (n in 1:H0@N) {
               g[,n] <- H0$mtvgarch[[n]]$tv@g
               h[,n] <- H0$mtvgarch[[n]]$garch@h
               beta[1,n] <- H0$mtvgarch[[n]]$garch$Estimated$pars["beta",1]
             }
             w <- e/sqrt(g)
             z <- w/sqrt(h)

             # Get x_garch
             x_garch <- .x_garch(w,H0,h,beta)

             # Get x_tv
             x_tv <- .x_tv(z,H0,g,h,beta)

             # Get v_rho, dlldrho_A
             rtn <- .v_rho(z,H0,H1,testOrder)
             v_rho <- rtn$v_rho
             dlldrho_A <- rtn$dlldrho_A

             # Get im_garch
             im_garch <- .im_garch(H0,x_garch)

             # Get im_garch_cor
             im_garch_cor <- .im_garch_cor(H0,x_garch,v_rho)

             # Get im_tv,
             im_tv <- .im_tv(H0,x_tv)

             # Get im_tv_cor
             im_tv_cor <- .im_tv_cor(H0,x_tv,v_rho)

             # Get im_tv_garch
             Pt <- H0$Estimated$P
             im_tv_garch <- .im_tv_garch(H0,Pt,x_tv,x_garch)

             # Get im_cor
             im_cor <- .im_cor(H0,Pt,v_rho)

             # Get LM using all InfoMatrix blocks
             IM_list <- list()
             IM_list$IM_tv <- im_tv
             IM_list$IM_tv_cor <- im_tv_cor
             IM_list$IM_garch <- im_garch
             IM_list$IM_garch_cor <- im_garch_cor
             IM_list$IM_tv_garch <- im_tv_garch
             IM_list$IM_cor <- im_cor

             LM <- .LM(H0,IM_list,dlldrho_A,testOrder)

             return(LM)

           }
)

##===  test.TVCC1vTVCC2 ===####
setGeneric(name="test.TVCC1vTVCC2",
           valueClass = "numeric",
           signature = c("e","H0","testOrder"),
           def = function(e,H0,testOrder){
             # Purpose:
             # Test a time-series to identify evidence of a second transition...

             # Validation
             objType <- class(H0)
             if(objType[1] != "stcc1_class"){
               warning("This test requires a valid instance of an estimated stcc model as the null (H0)")
               return(0)
             }

             # Get the common variables:
             g <- matrix(1,H0@Tobs,H0@N)
             h <- matrix(1,H0@Tobs,H0@N)
             beta <- matrix(1,1,H0@N)
             for (n in 1:H0@N) {
               g[,n] <- H0$mtvgarch[[n]]$tv@g
               h[,n] <- H0$mtvgarch[[n]]$garch@h
               beta[1,n] <- H0$mtvgarch[[n]]$garch$Estimated$pars["beta",1]
             }
             w <- e/sqrt(g)
             z <- w/sqrt(h)

             # Get x_tv - T x Total Nr.TvPars
             x_tv <- .x_tv(z,H0,g,h,beta)

             # Get x_garch - T x Total Nr.GarchPars
             x_garch <- .x_garch(w,H0,h,beta)

             # Get v_cor - T x (2 + TestOrder)
             v_cor_A <- .v_cor_H1(H0,testOrder) # T x testOrder
             v_cor <- cbind(.v_cor_H0(H0), v_cor_A)

             # T x 2 derivative of G w.r.t. tr pars (speed and loc)
             dGdtr <- .dG_dtr(H0,trNum = 1)

             I <- diag(H0@N)
             U <- .get_U(H0@N)

             # Note: this has been written with STCC (single transition) as H0, against one additional transition
             # TO DO: generalise to have STCC (with k transitions) as H0, against one more additional transition

             # Number of transitions in the Null model plus 1:
             trH0 <- 1 + 1

             im_cor_dim   <- (trH0 + testOrder) * H0@N * (H0@N-1)/2

             im_tv        <- matrix(0,NCOL(x_tv),NCOL(x_tv))
             im_tv_garch  <- matrix(0,NCOL(x_tv),NCOL(x_garch))
             im_tv_tr     <- matrix(0,NCOL(x_tv),H0@nr.trPars)
             im_tv_cor    <- matrix(0,NCOL(x_tv),im_cor_dim)

             im_garch     <- matrix(0,NCOL(x_garch),NCOL(x_garch))
             im_garch_tr  <- matrix(0,NCOL(x_garch),H0@nr.trPars)
             im_garch_cor <- matrix(0,NCOL(x_garch),im_cor_dim)

             im_tr        <- matrix(0,H0@nr.trPars,H0@nr.trPars)
             im_tr_cor    <- matrix(0,H0@nr.trPars,im_cor_dim)
             im_cor       <- matrix(0,im_cor_dim,im_cor_dim)

             dlldrho_A    <- matrix(0,testOrder * H0@N * (H0@N-1)/2,1)

             for(t in 1:H0@Tobs){
               Pt <- .unVecL(H0$Estimated$Pt[t,])
               Ptinv <- solve(Pt)
               P1 <- H0$Estimated$P1
               P2 <- H0$Estimated$P2

               # im_tv
               I.Pt.Ptinv <- .I.P.Pinv_scale(H0,Pt,"tv","tv")
               im_tv <- im_tv + (t(x_tv[t,,drop=FALSE]) %*% x_tv[t,,drop=FALSE]) * I.Pt.Ptinv

               # im_garch
               I.Pt.Ptinv <- .I.P.Pinv_scale(H0,Pt,"garch","garch")
               im_garch <- im_garch + (t(x_garch[t,,drop=FALSE]) %*% x_garch[t,,drop=FALSE]) * I.Pt.Ptinv

               # im_tv_garch
               I.Pt.Ptinv <- .I.P.Pinv_scale(H0,Pt,"tv","garch")
               im_tv_garch <- im_tv_garch + (t(x_tv[t,,drop=FALSE]) %*% x_garch[t,,drop=FALSE]) * I.Pt.Ptinv

               # im_tr
               Pinv.K <- .Pinv.K(H0,Pt)

               x_tr <- -0.5 * t(dGdtr[t,,drop=FALSE]) %*% matrix(.vec(P2-P1),1,H0@N^2)  # nr.trPars x N^2
               im_tr <- im_tr + 0.25 * x_tr %*% Pinv.K %*% t(x_tr)

               # im_tr_cor
               im_tr_cor <- im_tr_cor + x_tr %*% Pinv.K %*% (v_cor[t,,drop=FALSE] %x% U)

               # im_cor
               U.Pinv.K.U <- .U.Pinv.K.U(H0,Pt)
               im_cor <- im_cor + (t(v_cor[t,,drop=FALSE]) %*% v_cor[t,,drop=FALSE]) %x% U.Pinv.K.U

               # im_tv_cor AND im_garch_cor AND im_tv_tr AND im_garch_tr
               idxrange_tv <- 0
               idxrange_garch <- 0
               rowblock_tv_cor <- matrix(NA,nrow=0,ncol=(NCOL(v_cor)*NCOL(U)))
               rowblock_garch_cor <- matrix(NA,nrow=0,ncol=(NCOL(v_cor)*NCOL(U)))
               rowblock_tv_tr <- matrix(NA,nrow=0,ncol=NROW(x_tr))
               rowblock_garch_tr <- matrix(NA,nrow=0,ncol=NROW(x_tr))
               for (n in 1:H0@N){
                 endblock <- (Ptinv[n,,drop=FALSE]%x%I[n,,drop=FALSE]+I[n,,drop=FALSE]%x%Ptinv[n,,drop=FALSE])
                 endblock_cor <- endblock %*% (v_cor[t,,drop=FALSE]%x%U)
                 endblock_tr  <- endblock %*% t(x_tr)
                 idxrange_tv <- (max(idxrange_tv)+1):(max(idxrange_tv)+H0$mtvgarch[[n]]$tv@nr.pars)
                 idxrange_garch <- (max(idxrange_garch)+1):(max(idxrange_garch)+H0$mtvgarch[[n]]$garch@nr.pars)
                 rowblock_tv_cor <- rbind(rowblock_tv_cor,t(x_tv[t,idxrange_tv,drop=FALSE])%*%endblock_cor)
                 rowblock_garch_cor <- rbind(rowblock_garch_cor,t(x_garch[t,idxrange_garch,drop=FALSE])%*%endblock_cor)
                 rowblock_tv_tr <- rbind(rowblock_tv_tr,t(x_tv[t,idxrange_tv,drop=FALSE])%*%endblock_tr)
                 rowblock_garch_tr <- rbind(rowblock_garch_tr,t(x_garch[t,idxrange_garch,drop=FALSE])%*%endblock_tr)

               }
               im_tv_cor <- im_tv_cor + rowblock_tv_cor
               im_garch_cor <- im_garch_cor + rowblock_garch_cor
               im_tv_tr <- im_tv_tr + rowblock_tv_tr
               im_garch_tr <- im_garch_tr + rowblock_garch_tr

               # dlldrho_A: testOrder*N(N-1)/2 x 1
               vP_PxP_ZxZ <- .vec(Ptinv) - ((Ptinv %x% Ptinv) %*% t(z[t,,drop=FALSE] %x% z[t,,drop=FALSE]))
               dlldrho_A <- dlldrho_A + ((t(v_cor_A[t,,drop=FALSE]) %x% t(U) ) %*% vP_PxP_ZxZ)


             }  # End of for(t in 1:Tobs)

             im_tv <- im_tv/H0@Tobs
             im_tv_garch <- im_tv_garch/H0@Tobs
             im_tv_tr <- im_tv_tr/H0@Tobs
             im_tv_cor <- -0.5*im_tv_cor/H0@Tobs
             im_garch <- im_garch/H0@Tobs
             im_garch_tr <- im_garch_tr/H0@Tobs
             im_garch_cor <- -0.5*im_garch_cor/H0@Tobs
             im_tr <- im_tr/H0@Tobs
             im_tr_cor <- -0.5*im_tr_cor/H0@Tobs
             im_cor <- 0.25*im_cor/H0@Tobs
             dlldrho_A <- -0.5* dlldrho_A


             # Get LM using all InfoMatrix blocks
             IM_list <- list()
             IM_list$IM_tv <- im_tv
             IM_list$IM_tv_garch <- im_tv_garch
             IM_list$IM_tv_tr <- im_tv_tr
             IM_list$IM_tv_cor <- im_tv_cor
             IM_list$IM_garch <- im_garch
             IM_list$IM_garch_tr <- im_garch_tr
             IM_list$IM_garch_cor <- im_garch_cor
             IM_list$IM_tr <- im_tr
             IM_list$IM_tr_cor <- im_tr_cor
             IM_list$IM_cor <- im_cor

             LM <- .LM_v2(H0,IM_list,dlldrho_A,testOrder)

             return(LM)

           }
)




## ===== Test Sub Functions =====####

.get_U <- function(N){
  # Construct the U matrix:Dimensions = N^2 x N*(N-1)/2
  U <- NULL
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      block <- matrix(0,N,N)
      block[i,j] <- block[j,i] <- 1
      Ucol <- as.vector(block)
      U <- cbind(U,Ucol)
    }
  }
  return(U)
}

.get_K <- function(N){
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
  return(K)
}

##=== .I.P.Pinv_scale( tv or garch ) ===####
.I.P.Pinv_scale <- function(H0,Pt,type1,type2){

  # H0 is the Null-Hypothesis - must be an estimated CCC object
  # type: "tv" or "garch" - case sensitive

  N <- H0@N
  I <- diag(N,N) # NxN Identity matrix
  I.P.Pinv <- I + Pt * solve(Pt)
  I.P.Pinv_scale <- NULL

  for (i in 1:N) {
    # i = row index
    if(H0$mtvgarch[[i]][[type1]]@nr.pars > 0) {
      I.P.Pinv_scale_row <- NULL

      for (j in 1:N) {
        # j = col index
        if (H0$mtvgarch[[j]][[type2]]@nr.pars > 0){
          #TODO: Confirm this works when N contains different Garch Types
          scaleFactor <- matrix(1,H0$mtvgarch[[i]][[type1]]@nr.pars, H0$mtvgarch[[j]][[type2]]@nr.pars)
          I.P.Pinv_scale_row <- cbind(I.P.Pinv_scale_row, (I.P.Pinv[i,j] %x% scaleFactor))
        }
      } # End: for (j in 1:N)

      I.P.Pinv_scale <- rbind(I.P.Pinv_scale, I.P.Pinv_scale_row)
    }
  } # End: for (i in 1:N)
  return(I.P.Pinv_scale)

}

##===  .x_garch ===####
setGeneric(name=".x_garch",
           valueClass = "matrix",
           signature = c("w","H0","h","beta"),
           def = function(w,H0,h,beta){

             # partial derivatives of h1,h2,...,hN w.r.t garch_pars
             v_garch <- NULL
             h_scale <- NULL
             beta_scale <- NULL
             x_garch <- matrix(NA,0,0)  # Initialise return matrix

             # Loop over every series:
             for (n in 1:H0@N){

               if(H0$mtvgarch[[n]]$garch$type==garchtype$noGarch ) {
                 # Do nothing, move on to next series

               } else if(H0$mtvgarch[[n]]$garch$type==garchtype$general) {
                 # General GARCH(1,1) case
                 v_garch <- cbind(v_garch,c(0,rep(1,(H0@Tobs-1))),c(0,w[1:(H0@Tobs-1),n]^2),c(0,h[1:(H0@Tobs-1),n])) # T x Num_garch_pars, each row = "1~w(i,t-1)^2~h(i,t-1)", i=1,...,N
                 beta_scale <- beta_scale <- c(beta_scale,(beta[1,n] %x% c(1,1,1)))    # vector, length 3
                 scaleFactor <- matrix(1,nrow=1,ncol=3)
                 h_scale <- cbind(h_scale,h[,n,drop=FALSE] %x% scaleFactor)    # T x Num_garch_pars
               } else if (H0$mtvgarch[[n]]$garch$type == garchtype$gjr) {
                 # GJR GARCH(1,1) case
                 v_garch <- cbind(v_garch,c(0,rep(1,(H0@Tobs-1))),c(0,w[1:(H0@Tobs-1),n]^2),c(0,h[1:(H0@Tobs-1),n]),c(0,(min(w[1:(H0@Tobs-1),n],0))^2) ) # T x Num_garch_pars, each row = "1~w(i,t-1)^2~h(i,t-1)~min[w(i,t-1),0]^2", i=1,...,N
                 beta_scale <- beta_scale <- c(beta_scale,(beta[1,n] %x% c(1,1,1,1)))    # vector, length 3
                 scaleFactor <- matrix(1,nrow=1,ncol=4)
                 h_scale <- cbind(h_scale,h[,n,drop=FALSE] %x% scaleFactor)    # T x Num_garch_pars
                 #warning("Annastiina - please confirm this matrix is correct for GJR garch!")
               }

             } # End: For..loop

             if(is.null(v_garch) && is.null(beta_scale)){
               # All series have noGarch
               x_garch <- matrix(nrow=H0@Tobs,ncol=0)
             } else {
               dhdt <- .ar1.Filter(v_garch,beta_scale) # T x Num_garch_pars, each row = "(dh(i,t)/dtheta(i))' ", i=1,...,N
               x_garch <- -0.5*dhdt/h_scale # T x Num_garch_pars, each row = "x_it'", i=1,...,N
             }

             return(x_garch)

           }
)
##===  .x_tv ===####
setGeneric(name=".x_tv",
           valueClass = "matrix",
           signature = c("w","H0","g","h","beta"),
           def = function(w,H0,g,h,beta){

             v_tv <- NULL
             h_scale <- NULL
             g_scale <- NULL
             beta_scale <- NULL
             dgdt <- NULL
             dhdt <- NULL
             x_tv <- matrix(NA,0,0)  # Initialise return matrix

             # partial derivatives of g1,g2,...,gN w.r.t tv_pars
             for (n in 1:H0@N){

               if(H0$mtvgarch[[n]]$tv@nr.pars == 0 && H0$mtvgarch[[n]]$garch@nr.pars == 0){
                 break  #Goto next 'n' in for..loop
               }

               ## All code below is assured that there are tv parameters, so always calc dgdt

               dgdt_n <- .dg_dt(H0$mtvgarch[[n]]$tv)  # T x #of TVpars in TV[n], includes d0's derivative if it is a free param
               scaleFactor <- matrix(1,nrow=1,ncol=NCOL(dgdt_n))
               g_scale <- cbind( g_scale,g[,n,drop=FALSE] %x% scaleFactor ) # cbinds g(nt) as many times as g(n) has tv parameters
               dgdt <- cbind(dgdt,dgdt_n)

               if (H0$mtvgarch[[n]]$garch$type==garchtype$noGarch){
                 # Has TV, but noGarch => no beta_ or h_scale to calculate


               } else {
                 # all other GARCH(1,1) cases
                 v_tv_n <- ((-H0$mtvgarch[[n]]$garch$Estimated$pars["alpha",1] * c(0,1/g[1:(H0@Tobs-1),n])*c(0,w[1:(H0@Tobs-1),n]^2)) %x% scaleFactor) * dgdt_n
                 v_tv <- cbind(v_tv,v_tv_n)

                 beta_scale_n <- as.vector(beta[,n,drop=FALSE] %x% scaleFactor)
                 beta_scale <- c(beta_scale,beta_scale_n)

                 h_scale_n <- h[,n,drop=FALSE] %x% scaleFactor
                 h_scale <- cbind(h_scale,h_scale_n)

               }

             } # End: for loop

             if(is.null(g_scale) && is.null(h_scale)){
               # All series have noGarch AND noTV
               x_tv <- matrix(1,nrow=H0@Tobs,ncol=0)
             } else if(is.null(h_scale)){
               # All series have noGarch, but some do have TV
               x_tv <- -0.5*dgdt/g_scale
             } else {
               dhdt <- .ar1.Filter(v_tv,beta_scale) # T x Num_tv_pars, each row = dh(i,t).dtvpar(i), i=1...N
               x_tv <- -0.5*dhdt/h_scale -0.5*dgdt/g_scale  # T x Num_tv_pars, each row = x_it, i=1...N
             }

             return(x_tv)

           }
)

##===  .x_tau ===####
setGeneric(name=".x_tau",
           valueClass = "list",
           signature = c("z","H0","H1","testOrder"),
           def = function(z,H0,H1,testOrder){
             st <- H1$st
             if (testOrder==1) x_tau <- (-0.5)*cbind(st)
             if (testOrder==2) x_tau <- (-0.5)*cbind(st,st^2)
             if (testOrder==3) x_tau <- (-0.5)*cbind(st,st^2,st^3)

             N <- H0@N
             P <- H0$Estimated$P
             tmp <- eigen(P)
             Q <- tmp$vectors
             L.inv <- diag(tmp$values^(-1)) # matrix

             One_N_1.1 <- matrix(1,nrow=(N-1),ncol=1)
             I <- diag(N,N) # NxN Identity matrix
             I.P.Pinv <- I + P*solve(P)
             I.P.Pinv_scale <- NULL
             scaleFactor <- NULL

             # score for rho
             qxq <- matrix(0,nrow=N,ncol=(N^2))  # each row = qi' %x% qi', 1 x N^2, N rows
             for (n in 1:N){
               qxq[n,] <- t(Q[,n,drop=FALSE] %x% Q[,n,drop=FALSE])
             }
             mHelp0 <- matrix(0,nrow=H0@Tobs,ncol=(N-1))
             for (t in 1:H0@Tobs){
               vec_zz_t <- .vec(t(z[t,,drop=FALSE]) %*% (z[t,,drop=FALSE]))  # N^2 x 1
               mHelp0[t,] <-  L.inv[1:(N-1),1:(N-1)] %*% (One_N_1.1 - L.inv[1:(N-1),1:(N-1)] %*% qxq[1:(N-1),,drop=FALSE] %*% vec_zz_t) - One_N_1.1 %*% L.inv[N,N] %*% (1-L.inv[N,N] * (qxq[N,,drop=FALSE] %*% vec_zz_t))  # 1 x N-1
             }

             dlldrho_A <- t(mHelp0) %*% x_tau   # (N-1)xT %*% (T x testorder) = N-1 x testorder, SUM OVER TIME
             dlldrho_A <- .vec(dlldrho_A)  # testorder*(N-1) x 1, SUM OVER TIME
             x_tau <- cbind(rep(-0.5,H0@Tobs),x_tau) # T x 2 or T x 3, now add column of ones at the front (Note:xtau includes -0.5)

             rtn <- list()
             rtn$x_tau <- x_tau
             rtn$dlldrho_A <- dlldrho_A
             return(rtn)

           }
)

##===  .v_rho ===####
setGeneric(name=".v_rho",
           valueClass = "list",
           signature = c("z","H0","H1","testOrder"),
           def = function(z,H0,H1,testOrder){

             st <- H1@st
             if (testOrder==1) v_rho <- cbind(st)
             if (testOrder==2) v_rho <- cbind(st,st^2)
             if (testOrder==3) v_rho <- cbind(st,st^2,st^3)

             N <- H0@N
             P <- H0$Estimated$P
             Pinv <- solve(P)

             # U matrix: N^2 x N*(N-1)/2
             U <- .get_U(H0@N)

             # score for rho
             zKRONz <- matrix(0,nrow=N^2,ncol=H0@Tobs)
             for (t in 1:H0@Tobs) {
               zKRONz[,t] <- t(z[t,,drop=FALSE] %x% z[t,,drop=FALSE]) # (N^2 x T), each col = "z_t kron z_t"
             }
             scaleFactor <- matrix(1,nrow=1,ncol=H0@Tobs)
             dlldrho_A <- -0.5*t(U) %*% ( .vec(Pinv) %x% scaleFactor-(Pinv%x%Pinv)%*%(zKRONz) ) %*% v_rho # N*(N-1)/2 x testorder
             dlldrho_A <- .vec(dlldrho_A) # testorder*N*(N-1)/2 x 1, SUM OVER TIME
             v_rho <- cbind(1,v_rho) # T x 2 or T x 3, now add column of ones at the front

             rtn <- list()
             rtn$v_rho <- v_rho
             rtn$dlldrho_A <- dlldrho_A
             return(rtn)

           }
)

##===  .im_garch ===####
setGeneric(name=".im_garch",
           valueClass = "matrix",
           signature = c("H0","x_garch"),
           def = function(H0,x_garch){

             # IM_garch (Num_garch_pars x Num_garch_pars), SUM OVER TIME

             Pt <- H0$Estimated$P
             I.P.Pinv_scale <- .I.P.Pinv_scale(H0,Pt,"garch","garch")

             IM_garch <- matrix(NA,0,0)

             if(!is.null(I.P.Pinv_scale)) {
               IM_garch <- ((t(x_garch) %*% x_garch) * I.P.Pinv_scale) / H0@Tobs
             }
             return(IM_garch)

           }
)

##===  .im_garch_cor_parsim(...,x_tau) ===####
setGeneric(name=".im_garch_cor_parsim",
           valueClass = "matrix",
           signature = c("H0","x_garch","x_tau"),
           def = function(H0,x_garch,x_tau){

             # IM_garch (Num_garch_pars x Num_garch_pars), SUM OVER TIME
             N <- H0@N
             P <- H0$Estimated$P
             tmp <- eigen(P)
             Q <- tmp$vectors
             L.inv <- diag(tmp$values^(-1)) # matrix
             One_1.N_1 <- matrix(1,nrow=1,ncol=(N-1))

             # IM_garch_cor, will be 3N x (testorder+1)*(N-1)
             mHelp1 <- matrix(0,nrow=0,ncol=(N - 1) )

             for (n in 1:H0@N) {
               scaleFactor <- matrix(1, nrow=H0$mtvgarch[[n]]$garch@nr.pars, ncol=1)
               mHelp1 <- rbind(mHelp1, scaleFactor %x% (2*Q[n,1:(N-1),drop=FALSE]^2 %*% L.inv[1:(N-1),1:(N-1)] - 2*Q[n,N]^2 * L.inv[N,N] * One_1.N_1))
             }
             mHelp2 <- t(x_garch) %*% x_tau
             mHelp3 <- matrix(0,nrow=NROW(mHelp2),ncol=(N-1)*NCOL(mHelp2))
             for (i in 1:NROW(mHelp2)){
               mHelp3[i,] <- mHelp2[i,,drop=FALSE] %x% mHelp1[i,,drop=FALSE]
             }
             IM_garch_cor <- mHelp3/H0@Tobs # (Num_Garch_Pars x (testorder+1)*N), SUM OVER TIME

             return(IM_garch_cor)

           }
)

##===  .im_garch_cor(...,v_rho) ===####
setGeneric(name=".im_garch_cor",
           valueClass = "matrix",
           signature = c("H0","x_garch","v_rho"),
           def = function(H0,x_garch,v_rho){

             N <- H0@N
             P <- H0$Estimated$P
             Pinv <- solve(P)
             I <- diag(nrow = N,ncol = N) # NxN Identity matrix
             # U matrix: N^2 x N*(N-1)/2
             U <- .get_U(H0@N)

             # IM_garch_cor, 3N x (testorder+1)*N*(N-1)/2
             mHelp1 <- matrix(0,nrow = 0,ncol=N*N)  # N x N^2
             for (n in 1:N) mHelp1 <- rbind(mHelp1, (Pinv[n,] %x% I[n,] + I[n,] %x% Pinv[n,])) # N x N^2

             mHelp1 <- -0.5*(mHelp1 %*% U)
             mHelp1_scale <- matrix(0,nrow = 0,ncol=N*(N-1)/2)  # (Num_Garch_Pars x N*(N-1)/2), required for rbind() below to work
             for (n in 1:N) {
               scaleFactor <- matrix(1, nrow=H0$mtvgarch[[n]]$garch@nr.pars, ncol=1)
               mHelp1_scale <- rbind(mHelp1_scale, (mHelp1[n, ,drop=FALSE] %x% scaleFactor))
             }

             mHelp2 <- t(t(v_rho) %*% x_garch)/H0@Tobs # Num_garch_pars x 2 (or x3 if TestOrder=2), SUM OVER TIME
             IM_garch_cor <- matrix(NA,nrow=NROW(mHelp2),ncol=0)  # (Num_Garch_Pars x (testorder+1)*N*(N-1)/2)
             for (i in 1:NCOL(mHelp2)){
               IM_garch_cor <- cbind(IM_garch_cor,(mHelp2[,i,drop=FALSE] %x% t(rep(1,(N*(N-1))/2))) * mHelp1_scale)  # (Num_Garch_Pars x (testorder+1)*N*(N-1)/2),  SUM OVER TIME
             }
             return(IM_garch_cor)

           }
)






##===  .im_tv ===####
setGeneric(name=".im_tv",
           valueClass = "matrix",
           signature = c("H0","x_tv"),
           def = function(H0,x_tv){

             Pt <- H0$Estimated$P
             I.P.Pinv_scale <- .I.P.Pinv_scale(H0,Pt,"tv","tv")

             IM_tv <- matrix(NA,0,0)

             if(!is.null(I.P.Pinv_scale)) {
               IM_tv <- ((t(x_tv) %*% x_tv) * I.P.Pinv_scale) / H0@Tobs  # sum(num_tv_pars) x sum(num_tv_pars)
             }

             return(IM_tv)
           }
)

##===  .im_tv_cor_parsim(...,x_tau) ===####
setGeneric(name=".im_tv_cor_parsim",
           valueClass = "matrix",
           signature = c("H0","x_tv","x_tau"),
           def = function(H0,x_tv,x_tau){

             N <- H0@N
             P <- H0$Estimated$P
             tmp <- eigen(P)
             Q <- tmp$vectors
             L.inv <- diag(tmp$values^(-1)) # matrix
             One_1.N_1 <- matrix(1,nrow=1,ncol=(N-1))

             # IM_tv_cor, (total #of tv pars in model) x (testorder+1)*(N-1)
             mHelp1 <- matrix(0,nrow=0,ncol=(N-1) )
             for (n in 1:N) {
               if (H0$mtvgarch[[n]]$tv@nr.pars > 0){
                 scaleFactor <- matrix(1,nrow=H0$mtvgarch[[n]]$tv@nr.pars,ncol=1)
                 mHelp1 <- rbind(mHelp1, scaleFactor %x% (2*Q[n,1:(N-1),drop=FALSE]^2 %*% L.inv[1:(N-1),1:(N-1)] - 2*Q[n,N]^2 * L.inv[N,N] * One_1.N_1) )
               }
             } # End: for (n in 1:N)
             mHelp2 <- t(x_tv) %*% x_tau
             mHelp3 <- matrix(0,nrow=NROW(mHelp2),ncol=(N-1)*NCOL(mHelp2))
             if (NROW(mHelp3)>0){
               for (i in 1:NROW(mHelp3)){
                 mHelp3[i,] <- mHelp2[i,,drop=FALSE] %x% mHelp1[i,,drop=FALSE]
               }
             }
             IM_tv_cor <- mHelp3 / H0@Tobs  # (Num_tv_Pars x (testorder+1)*N),  SUM OVER TIME

             return(IM_tv_cor)

           }
)

##===  .im_tv_cor(...,v_rho) ===####
setGeneric(name=".im_tv_cor",
           valueClass = "matrix",
           signature = c("H0","x_tv","v_rho"),
           def = function(H0,x_tv,v_rho){

             N <- H0@N
             Pinv <- solve(H0$Estimated$P)
             I <- diag(nrow = N,ncol = N) # NxN Identity matrix
             # U matrix: N^2 x N*(N-1)/2
             U <- .get_U(H0@N)

             # IM_tv_cor, (#of tv pars in n,n=1...N)N x (testorder+1)*N*(N-1)/2
             mHelp1 <- matrix(0,nrow=0,ncol=N*N )
             for (n in 1:N) {
               if (H0$mtvgarch[[n]]$tv@nr.pars > 0){
                 mHelp1 <- rbind(mHelp1, (Pinv[n,,drop=FALSE] %x% I[n,,drop=FALSE] + I[n,,drop=FALSE] %x% Pinv[n,,drop=FALSE])) # N x N^2
               }
             } # End: for (n in 1:N)
             mHelp1 <- -0.5*(mHelp1 %*% U)  # N x N*(N-1)/2
             mHelp1_scale <- matrix(0,nrow=0,ncol=(N*(N-1)/2))  # Total Num_tv_Pars x N*(N-1)/2
             for (n in 1:N) {
               if (H0$mtvgarch[[n]]$tv@nr.pars > 0){
                 scaleFactor <- matrix(1,nrow=H0$mtvgarch[[n]]$tv@nr.pars,ncol=1)
                 mHelp1_scale <- rbind(mHelp1_scale, mHelp1[n,,drop=FALSE] %x% scaleFactor)
               }
             }
             mHelp2 <- t(t(v_rho)%*%x_tv) / H0@Tobs # Num_tv_pars x 2 (or x3 if TestOrder=2), SUM OVER TIME

             IM_tv_cor <- matrix(NA,nrow = NROW(mHelp2),ncol = 0)
             for (i in 1:NCOL(mHelp2)){
               scaleFactor <- matrix(1,nrow=1,ncol=N*(N-1)/2)
               IM_tv_cor <- cbind(IM_tv_cor,(mHelp2[,i,drop=FALSE]%x%scaleFactor) * mHelp1_scale)  # (Num_tv_Pars x (testorder+1)*N*(N-1)/2),  SUM OVER TIME
             }

             return(IM_tv_cor)

           }
)


##===  .im_tv_garch ===####
setGeneric(name=".im_tv_garch",
           valueClass = "matrix",
           signature = c("H0","Pt","x_tv","x_garch"),
           def = function(H0,Pt,x_tv,x_garch){

             N <- H0@N

             I <- diag(nrow = N,ncol = N) # NxN Identity matrix
             I.P.Pinv <- I + Pt * solve(Pt)
             I.P.Pinv_scale <- .I.P.Pinv_scale(H0,Pt,"tv","garch")

             IM_tv_garch <- (t(x_tv) %*% x_garch) * (I.P.Pinv_scale/H0@Tobs)  #  (total # of tvpars)x(total # of garchpars)

             return(IM_tv_garch)

           }
)
##===  .im_cor_parsim(...,x_tau) ===####
setGeneric(name=".im_cor_parsim",
           valueClass = "matrix",
           signature = c("H0","x_tau","testOrder"),
           def = function(H0,x_tau,testOrder){

             # IM_cor (testorder+1)*(N-1) x (testorder+1)*(N-1), SUM OVER TIME
             N <- H0@N
             P <- H0$Estimated$P
             tmp <- eigen(P)
             L2.inv <- diag(tmp$values^(-2)) # matrix
             One_N_1.N_1 <- matrix(1,nrow=(N-1),ncol=(N-1))

             if (testOrder==1){
               mHelp7 <- matrix(c(1,1/2,1/2,1/3),nrow=2,ncol=2)
             }
             if (testOrder==2){
               mHelp7 <- matrix(c(1,1/2,1/3,1/2,1/3,1/4,1/3,1/4,1/5),nrow=3,ncol=3)
             }
             if (testOrder==3){
               mHelp7 <- matrix(c(1,1/2,1/3,1/4,1/2,1/3,1/4,1/5,1/3,1/4,1/5,1/6,1/4,1/5,1/6,1/7),nrow=4,ncol=4)
             }

             IM_cor <- ( t(x_tau)%*%x_tau ) %x% (2*L2.inv[1:(N-1),1:(N-1)] + 2*L2.inv[N,N]*One_N_1.N_1) / H0@Tobs # (testorder+1)*N x (testorder+1)*N

             return(IM_cor)

           }
)

##===  .U.Pinv.K.U(H0,Pt) ===####
setGeneric(name=".U.Pinv.K.U",
           valueClass = "matrix",
           signature = c("H0","Pt"),
           def = function(H0,Pt){

             N <- H0@N
             Pinv <- solve(Pt)
             I <- diag(nrow = N,ncol = N) # NxN Identity matrix

             # U matrix: N^2 x N*(N-1)/2
             U <- .get_U(N)

             # K matrix: N^2 x N^2
             K <- .get_K(N)

             # N*(N-1)/2 x N*(N-1)/2
             # TODO: Fix Error = (Pinv %x% I) %*% K : non-conformable arguments
             return( t(U) %*% (Pinv %x% Pinv + (Pinv %x% I) %*% K %*% (Pinv %x% I)) %*% U )

           }
)

##===  .Pinv.K(...,v_rho) ===####
setGeneric(name=".Pinv.K",
           valueClass = "matrix",
           signature = c("H0","Pt"),
           def = function(H0,Pt){

             N <- H0@N
             Pinv <- solve(Pt)
             I <- diag(nrow = N,ncol = N) # NxN Identity matrix

             # K matrix: N^2 x N^2
             K <- .get_K(N)

             # N*(N-1)/2 x N*(N-1)/2
             return( Pinv %x% Pinv + (Pinv %x% I) %*% K %*% (Pinv %x% I) )

           }
)

##===  .im_cor(...,v_rho) ===####
setGeneric(name=".im_cor",
           valueClass = "matrix",
           signature = c("H0","Pt","v_rho"),
           def = function(H0,Pt,v_rho){

             mHelp1 <- .U.Pinv.K.U(H0,Pt)
             mHelp2 <- t(v_rho) %*% v_rho/H0@Tobs # 2x2 or 3x3, SUM OVER TIME
             IM_cor <- 0.25*(mHelp2 %x% mHelp1) # (testorder+1)*N*(N-1)/2 x (testorder+1)*N*(N-1)/2

             return(IM_cor)

           }
)
##===  .LM ===####
setGeneric(name=".LM",
           valueClass = "numeric",
           signature = c("H0","IM_list","dlldrho_A","testOrder"),
           def = function(H0,IM_list,dlldrho_A,testOrder){

             IM_tv <- IM_list$IM_tv
             IM_tv_cor <- IM_list$IM_tv_cor
             IM_garch <- IM_list$IM_garch
             IM_garch_cor <- IM_list$IM_garch_cor
             IM_tv_garch <- IM_list$IM_tv_garch
             IM_cor <- IM_list$IM_cor

             #IM <- rbind(IM_TV, IM_GARCH, IM_COR)
             IM <- rbind(cbind(IM_tv,IM_tv_garch,IM_tv_cor), cbind(t(IM_tv_garch),IM_garch,IM_garch_cor), cbind(t(IM_tv_cor),t(IM_garch_cor),IM_cor)) # the whole IM, not really needed...
             IM_inv <- solve(IM)
             ## The above is not necessarily very efficient - could use block inversion methods - TO DO LATER

             if (!is.matrix(IM_inv)) stop("LM Test: Can't invert the information matrix IM_inv")

             ##--- Block corresponding to the corr.parameters that are set to zero under null ---##
             block_start <- 1 + NCOL(IM_inv) - NROW(dlldrho_A)
             block_end <- NCOL(IM_inv)
             IM_inv_SE <- IM_inv[(block_start:block_end),(block_start:block_end)]

             ##--- Return LM ---##
             LM <- (1/H0@Tobs)*t(dlldrho_A) %*% IM_inv_SE %*% dlldrho_A
             return( as.numeric(LM[1,1]) )

           }
)

##===  .LM_v2 ===####
setGeneric(name=".LM_v2",
           valueClass = "numeric",
           signature = c("H0","IM_list","dlldrho_A","testOrder"),
           def = function(H0,IM_list,dlldrho_A,testOrder){


             IM_tv <- IM_list$IM_tv
             IM_tv_garch <- IM_list$IM_tv_garch
             IM_tv_tr <- IM_list$IM_tv_tr
             IM_tv_cor <- IM_list$IM_tv_cor
             IM_garch <- IM_list$IM_garch
             IM_garch_tr <- IM_list$IM_garch_tr
             IM_garch_cor <- IM_list$IM_garch_cor
             IM_tr <- IM_list$IM_tr
             IM_tr_cor <- IM_list$IM_tr_cor
             IM_cor <- IM_list$IM_cor

             #IM <- rbind(IM_TV, IM_GARCH, IM_COR)
             IM <- rbind(cbind(IM_tv,IM_tv_garch,IM_tv_tr,IM_tv_cor), cbind(t(IM_tv_garch),IM_garch,IM_garch_tr,IM_garch_cor),cbind(t(IM_tv_tr),t(IM_garch_tr),IM_tr,IM_tr_cor), cbind(t(IM_tv_cor),t(IM_garch_cor),t(IM_tr_cor),IM_cor)) # the whole IM, not really needed...
             IM_inv <- solve(IM)
             ## The above is not necessarily very efficient - could use block inversion methods - TO DO LATER

             if (!is.matrix(IM_inv)) stop("LM Test: Can't invert the information matrix IM_inv")

             ##--- Block corresponding to the corr.parameters that are set to zero under null ---##
             block_start <- 1 + NCOL(IM_inv) - NROW(dlldrho_A)
             block_end <- NCOL(IM_inv)
             IM_inv_SE <- IM_inv[(block_start:block_end),(block_start:block_end)]

             ##--- Return LM ---##
             LM <- (1/H0@Tobs)*t(dlldrho_A) %*% IM_inv_SE %*% dlldrho_A
             return( as.numeric(LM[1,1]) )

           }
)


##===  .v_cor_H0 ===####
setGeneric(name=".v_cor_H0",
           valueClass = "matrix",
           signature = c("H0"),
           def = function(H0){
             ret <- matrix(nrow = H0@Tobs,ncol = 2)
             ret[,1] <- 1 - calc.Gt(H0)
             ret[,2] <- calc.Gt(H0)
             return(ret)

           }
)
##===  .v_cor_H1 ===####
setGeneric(name=".v_cor_H1",
           valueClass = "matrix",
           signature = c("H0","testOrder"),
           def = function(H0,testOrder){
             ret <- matrix(nrow = H0@Tobs,ncol = testOrder)
             # Test requires the transition variable to have zero mean (by design)
             H1_st <- H0@st - mean(H0@st)
             for(n in 1:testOrder){
               ret[,n] <- H1_st^n
             }
             return(ret)
           }
)
