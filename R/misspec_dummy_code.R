missTest1(one tvgarch obj at a time){
#
g <- # T x 1
h <- # T x 1

# derivatives:
# (1/ht)*(dh/dtvpars)
dg_dtv <- # T x nr.tv.pars
dh_dtv <- # T x nr.tv.pars
dh_dga <- # T x nr.garch.pars
df_dli <- cbind(matrix(1,Tobs,1),st,st^2,st^3) # T x (testorder+1)

# matrices
u <- g^(-1) * dg_dtv     # (1/g)*(dg/dtvpars); T x nr.tv.pars
x1 <- (h^(-1)) * dh_dtv  # (1/ht)*(dh/dtvpars); T x nr.tv.pars
x2 <- (h^(-1)) * dh_dga  # (1/ht)*(dh/dgarchpars); T x nr.garch.pars
#x3 <- (h^(-1)) * dh_dli # (1/ht)*(dh/dlinpars); T x nr.lin.pars  DON'T THINK WE NEED THIS AT ALL
v3 <- (g^(-1)) * df_dli  # (1/gt)*(df/dlinpars); T x nr.lin.pars (=testorder+1)
r1 <- cbind(u+x1,x2)     # T x (tot.tv+garch.pars)
r2 <- v3                 # T x (testorder+1)

# 1 estimate tvgjrgarch, get std.residuals z=e/sqrt(gh), get SSR0 = sum(z^2-1)
z2_1 <- (e/sqrt(h*g))^2-1
SSR0 <- t(z2_1)%*%z2_1
# 2: regress z2_1 on r1~r2, get SSR
X <- cbind(r1,r2) # T x ncol(r1)+nocl(r2)
Y <- z2_1# T x 1
b = solve(t(X)%*%X)%*%t(X)%*%Y  # vector len=ncol(X)
resid = Y-X%*%t(b)   # T x 1
SSR1 <- t(resid)%*%resid # 1x1
# LM stat
LM <- Tobs * (SSR0-SSR1)/SSR0
p <- pchisq(LM,df=NCOL(v3),lower.tail=FALSE)

# ROBUST
# 1 as above
# 2 regress r2 on r1, get residual vectors
resid <- matrix(0,Tobs,NCOL(r2))
X = r1
for (i in 1:NCOL(r2)){
  Y = r2[,i]
  b = solve(t(X)%*%X)%*%t(X)%*%Y  # vector len=ncol(X)
  resid[,i] <- Y-X%*%t(b)   # T x 1
}
# regress 1 on (z2_1)resid, get SSR
Y = matrix(1,Tobs,1)
X <- z2_1*resid
b = solve(t(X)%*%X)%*%t(X)%*%Y  # vector len=ncol(X)
resid = Y-X%*%t(b)   # T x 1
SSR <- t(resid)%*%resid # 1x1
# LM Robust
LMrob <- Tobs-SSR
prob <- (LMrob,df=NCOL(v3),lower.tail=FALSE)
}
