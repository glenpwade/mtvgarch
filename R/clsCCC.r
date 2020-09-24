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
             this$P <- cor(z)

             return(this)

           }
)

## --- estimateCCC --- ####
setGeneric(name="estimateCCC",
           valueClass = "ccc_class",
           signature = c("z","cccObj","estimationCtrl"),
           def = function(z,cccObj,estimationCtrl){
             this <- cccObj

             this$Estimated <- list()
             #this$Estimated$P <- cor(z)
             #this$Estimated$error <- FALSE

             optimpars <- this$P[lower.tri(this$P)]

             this$Estimated <- .loglik.ccc(optimpars,z,this)

             return(this)
           }
)



##====  TESTS  ====####

##===  test.CCCParsim ===####

setGeneric(name="test.CCCParsim",
           valueClass = "list",
           signature = c("e","H0","H1","testOrder"),
           def = function(e,H0,H1,testOrder){

             # Validation
             objType <- class(H0)
             if(objType[1] != "ccc_class"){
               warning("This test requires a valid instance of an estimated ccc model as the null (H0)")
               return(list())
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

             # Get im_garch, im_garch_cor
             rtn <- .im_garch(H0,x_garch,x_tau)
             im_garch <- rtn$im_garch
             im_garch_cor <- rtn$im_garch_cor

             # Get im_tv, im_tv_cor
             rtn <- .im_tv(H0,x_tv,x_tau)
             im_tv <- rtn$im_tv
             im_tv_cor <- rtn$im_tv_cor

             # Get im_tv_garch
             im_tv_garch <- .im_tv_garch(H0,x_tv,x_garch)

             # Get im_cor
             im_cor <- .im_cor(H0,x_tau,testOrder)

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


## ===== Test Sub Functions =====####



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
                 v_garch <- cbind(v_garch,c(0,rep(1,(H0@Tobs-1))),c(0,w[1:(H0@Tobs-1),n]^2),c(0,h[1:(H0@Tobs-1),n])) # T x Num_garch_pars, each row = "1~w(i,t-1)^2~h(i,t-1)", i=1,...,N
                 beta_scale <- beta_scale <- c(beta_scale,(beta[1,n] %x% c(1,1,1)))    # vector, length 3
                 scaleFactor <- matrix(1,nrow=1,ncol=3)
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

##===  .im_garch ===####
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
             One_1.N_1 <- matrix(1,nrow=1,ncol=(N-1))
             I <- diag(N,N) # NxN Identity matrix
             I.P.Pinv <- I + P*solve(P)
             I.P.Pinv_scale <- NULL
             scaleFactor <- NULL

             for (i in 1:N) {
               # i = row index
               I.P.Pinv_scale_row <- NULL
               if (H0$mtvgarch[[i]]$garch$type == garchtype$noGarch){
                 # do nothing
               } else {
                 for (j in 1:N) {
                   if(H0$mtvgarch[[j]]$garch$type == garchtype$noGarch){
                     # do nothing
                   } else {
                     # j = col index
                     #if(H0$mtvgarch[[i]]$garch$type==garchtype$general && H0$mtvgarch[[j]]$garch$type==garchtype$general) scaleFactor <- One_33 else stop("Garch type is not supported")
                     scaleFactor <- One_33
                     I.P.Pinv_scale_row <- cbind(I.P.Pinv_scale_row, (I.P.Pinv[i,j] %x% scaleFactor))
                   }
                 } # End: for "j" loop
               }
               I.P.Pinv_scale <- rbind(I.P.Pinv_scale,I.P.Pinv_scale_row)
             } # End: for "i" loop


             if(!is.null(I.P.Pinv_scale)){
               IM_garch <- ((t(x_garch) %*% x_garch) * I.P.Pinv_scale) / H0@Tobs
               # IM_garch_cor, will be 3N x (testorder+1)*(N-1)
               mHelp1 <- matrix(0,nrow=0,ncol=(N - 1) )

               for (n in 1:H0@N) {
                 ## TODO: What about GJR??
                 if (H0$mtvgarch[[n]]$garch$type != garchtype$noGarch) {
                   scaleFactor <- matrix(1,nrow=3,ncol=1)
                   mHelp1 <- rbind(mHelp1, scaleFactor %x% (2*Q[n,1:(N-1),drop=FALSE]^2 %*% L.inv[1:(N-1),1:(N-1)] - 2*Q[n,N]^2 * L.inv[N,N] * One_1.N_1))
                 } else mHelp1 <- matrix(1,nrow = (3*Ho@N), ncol=(N-1))
               }
               mHelp2 <- t(x_garch) %*% x_tau
               mHelp3 <- matrix(0,nrow=NROW(mHelp2),ncol=(N-1)*NCOL(mHelp2))
               for (i in 1:NROW(mHelp2)){
                 mHelp3[i,] <- mHelp2[i,,drop=FALSE] %x% mHelp1[i,,drop=FALSE]
               }
               IM_garch_cor <- mHelp3/H0@Tobs # (Num_Garch_Pars x (testorder+1)*N), SUM OVER TIME

             } # End: if(!is.null(I.P.Pinv_scale))  => Only process when we have 'some' Garch

             rtn <- list()
             rtn$im_garch <- IM_garch
             rtn$im_garch_cor <- IM_garch_cor
             return(rtn)

           }
)

##===  .im_tv ===####
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

             I <- diag(N,N) # NxN Identity matrix
             I.P.Pinv <- I + P * solve(P)
             I.P.Pinv_scale <- NULL
             scaleFactor <- NULL

             for (i in 1:N) {
               # i = row index
               if(H0$mtvgarch[[i]]$tv@nr.pars > 0) {
                I.P.Pinv_scale_row <- NULL

                for (j in 1:N) {
                   # j = col index
                   if (H0$mtvgarch[[j]]$tv@nr.pars > 0){
                     scaleFactor <- matrix(1,H0$mtvgarch[[i]]$tv@nr.pars,H0$mtvgarch[[j]]$tv@nr.pars)
                     I.P.Pinv_scale_row <- cbind(I.P.Pinv_scale_row, (I.P.Pinv[i,j] %x% scaleFactor))
                   } # End: if (jj==0)
                 } # End: for (j in 1:N)
                 I.P.Pinv_scale <- rbind(I.P.Pinv_scale, I.P.Pinv_scale_row)

               } # End: if (ii==0)
             } # End: for (i in 1:N)

             IM_tv <- ((t(x_tv) %*% x_tv) * I.P.Pinv_scale) / H0@Tobs  # sum(num_tv_pars) x sum(num_tv_pars)


             # IM_tv_cor, (total #of tv pars in model) x (testorder+1)*(N-1)
             mHelp4 <- matrix(0,nrow=0,ncol=(N-1) )
             for (n in 1:N) {
               if (H0$mtvgarch[[n]]$tv@nr.pars > 0){
                 scaleFactor <- matrix(1,nrow=H0$mtvgarch[[n]]$tv@nr.pars,ncol=1)
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
             rtn$im_tv <- IM_tv
             rtn$im_tv_cor <- IM_tv_cor
             return(rtn)


           }
)
##===  .im_tv_garch ===####
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
               if (H0$mtvgarch[[i]]$tv@nr.pars > 0){
                 I.P.Pinv_scale_row <- NULL
                 for (j in 1:N){
                   # j = col index
                   if (H0$mtvgarch[[j]]$garch@nr.pars > 0){
                     #TODO: Hack below forces garch@nr.pars==3 (Need to work out how to handle GJR Garch)
                     #scaleFactor <- matrix(1,H0$mtvgarch[[i]]$tv@nr.pars,H0$mtvgarch[[j]]$garch@nr.pars)
                     scaleFactor <- matrix(1,H0$mtvgarch[[i]]$tv@nr.pars, 3)
                     I.P.Pinv_scale_row <- cbind(I.P.Pinv_scale_row, (I.P.Pinv[i,j] %x% scaleFactor))
                   } # End if (jj==0)
                 } # End for (j in 1:N)
                 I.P.Pinv_scale <- rbind(I.P.Pinv_scale,I.P.Pinv_scale_row)  # (total of tvpars) x (total of garchpars)
               } # End if (ii==0)
             }
             #IM_tv_garch <- NULL
             IM_tv_garch <- (t(x_tv) %*% x_garch) * (I.P.Pinv_scale/H0@Tobs)  #  (total # of tvpars)x(total # of garchpars)

             return(IM_tv_garch)

           }
)
##===  .im_cor ===####
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
##===  .LM ===####
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

             #IM <- rbind(IM_TV, IM_GARCH, IM_COR)
             IM <- rbind(cbind(IM_tv,IM_tv_garch,IM_tv_cor), cbind(t(IM_tv_garch),IM_garch,IM_garch_cor), cbind(t(IM_tv_cor),t(IM_garch_cor),IM_cor)) # the whole IM, not really needed...
             IM_inv <- solve(IM)
             ## The above is not necessarily very efficient - could use block inversion methods - TO DO LATER

             if (!is.matrix(IM_inv)) stop("LM Test: Can't invert the information matrix IM_inv")

             ##--- Block corresponding to the corr.parameters that are set to zero under null ---##
             block_start <- NCOL(IM_inv) - (testOrder*(N-1) + 1)
             block_end <- NCOL(IM_inv)
             IM_inv_SE <- IM_inv[(block_start:block_end),(block_start:block_end)]

             ##--- Return LM ---##
             return( (1/H0@Tobs)*t(dlldrho_A) %*% IM_inv_SE %*% dlldrho_A )

           }
)

