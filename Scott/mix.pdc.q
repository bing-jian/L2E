# June 15, 2004    few tweeks

# Copyright David W. Scott     June 10, 2004      general mixture K components via L2E
# Copyright David W. Scott     February 5, 2004   
# June 10, 2001

# Mixture of Multivariate Partial Density Component Estimation (mpdc)

# Purpose:  fit model   sum k=1,K w_k N(mu_k,sig_k)  by L2E

# Inputs:
#          X - n x d input data matrix (d=1 vector OK)
#          K - number of components desired (K=1 default)
#       grps - optional way of inputting initial guesses  (data labels 1,2,...,K)
#                  vector of length n    labels = 0 are ignored (useful if w.sum=F)
#      w.sum - constrain weights to sum to 1  (T/F)
#         mu - input guess for means  d x K (matrix)
#        sig - d x d x K input guess for covariance matrices (array)
#          w - input guess for K weights (always length K, even with w.sum=T)
#        nit - max number of iterations for nlmin
#        nev - max number of function evaluations for nlmin
#         pl - print level for nlmin (0=none 1=some 2=lots)

# Output:
#          list containing estimated parameters (m=mean sig=sig w=w)

# Note:  set parameter  br  to positive value to invoke browser

br <- 0

mix.pdc <- function(X,K=1,grps,w.sum=T,mu,sig,w,nit=100,nev=200,pl=1) {
  if(!is.matrix(X)) { X <- cbind(X) }; n <- nrow(X); d <- ncol(X)
  assign("X",X,f=1); assign("n",n,f=1); assign("d",d,f=1); assign("K",K,f=1)
  assign("onen",rep(1,n),f=1); assign("oned",rep(1,d),f=1)

  if(!missing(grps)) { if(length(grps)!=n) stop("grps length wrong")
    if(max(grps)!=K) print("warning -- max of grps does not match K")
    nk <- rep(0,K); mu <- matrix(0,d,K); sig <- array(0,c(d,d,K))
    for(k in 1:K) { ii <- seq(n)[grps==k]; nk[k] <- length(ii)
      if( length(ii)<2 ) stop(paste("grps insufficient for class ",as.character(k)))
      mu[,k] <- apply(X[ii,,drop=F],2,"mean"); sig[,,k] <- var(X[ii,]) }
    w <- nk/n }

  if(d==1) { mu <- matrix(c(mu),1,K); sig <- array(c(sig),c(1,1,K)) }
  if(any(dim(mu)!=c(d,K))) stop("input mu matrix wrong dimension")
  if(any(dim(sig)!=c(d,d,K))) stop("input sig array wrong dimension")
  if(length(w)!=K) stop("input w vector wrong length")

 pdc <- function(x) {		  # first extract parameters from  x
    mu <- matrix(x[1:(d*K)],d,K)  # assumes fortran order
    ns <- d*(d+1)/2; sigi.u <- matrix(x[(d*K+1):(d*K+ns*K)],ns,K) # will make an array
    tw <- x[-(1:(d*K+ns*K))]; if(length(tw)==K) { w <- exp(tw) } else
      { tw <- c(exp(tw),1); w <- tw/sum(tw) }	# if length=K then no sum constraint

    tot1 <- 0; tot2 <- 0; cc <- 2^(d/2); sig <- array(0,c(d,d,K))
    deti <- rep(0,K)   # determinant inverse sq root
    for(k in 1:K) { muk <- mu[,k]; U <- matrix(0,d,d)
        U[row(U)<=col(U)] <- sigi.u[,k]; deti[k] <- exp(sum(diag(U)))
        dU <- exp(diag(U)); diag(U) <- dU; sig[,,k] <- solve( t(U)%*%U )
       tot1 <- tot1 + w[k]^2*deti[k]/cc
       tot2 <- tot2 + w[k] * deti[k] *
             sum( exp( (-.5* ( (X-outer(onen,muk))%*%t(U) )**2) %*%oned ) )
 if(br>3) { print(paste("in loop k = ",as.ch(k))); browser() }
      if(k>1) { for(m in 1:(k-1)) { mum <- mu[,m]; sigi.km <- solve( sig[,,k]+sig[,,m] )
        U0 <- chol( (sigi.km+t(sigi.km))/2 );  deti.km <- exp(sum(log(diag(U0))))
        dd <- U0%*%(muk-mum); tot2 <- tot2 + 2*w[k]*w[m]*deti.km*exp(-.5*sum(dd^2))
 if(br>3) { print(paste("in loop m = ",as.ch(m))); browser() }  
  } } }
    tot <- tot1 - 2*tot2/n
 if(br>2) {print("exciting pdc"); browser()}
    tot
 }

  x0 <- c(mu)	# assumes fortran column order
  for(k in 1:K) { sig0 <- sig[,,k]
    if(d==1) { xU0 <- log(1/sqrt(sig0)) } else { tmp <- solve(sig0); sigi <- (tmp+t(tmp))/2
      U0 <- chol(sigi); diag(U0) <- log(diag(U0)); xU0 <- U0[row(U0)<=col(U0)]}
    x0 <- c(x0,xU0) }
  if(w.sum) { if(K==1) {w0 <- NULL} else { w0 <- log(w[-K]/w[K]) } } else { w0 <- log(w) }
    x0 <- c(x0,w0)	# w0 of length K-1 if sum constraint on (NULL if K=1)

if(br>1) {print("calling nlmin");browser()}

  ans <- nlmin(pdc,x0,print.level=pl,max.fcal=nev,max.iter=nit)

  xx <- ans$x; ns <- d*(d+1)/2
  mu <- matrix(xx[1:(d*K)],d,K); xx <- xx[-(1:(d*K))]
  for(k in 1:K) {
    U.ans <- matrix(0,d,d); U.ans[row(U.ans)<=col(U.ans)] <- xx[1:ns]
    diag(U.ans) <- exp(diag(U.ans)); sigi <- t(U.ans) %*% U.ans;
    sig[,,k] <- solve(sigi); xx <- xx[-(1:ns)] }
  w <- xx; if(length(w)==K) { w <- exp(w) } else { w <- exp(c(w,0)); w <- w/sum(w) }
if(br>0) {print("wrapping up");browser()}

  list(m=mu,sig=sig,w=w)
}
