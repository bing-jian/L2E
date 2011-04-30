## Copyright David W. Scott     February 5, 2004   
## June 10, 2001

## Multivariate Partial Density Component Estimation (mpdc)

## Purpose:  fit model   w N(mu,sig)

## Inputs:
##          x - n x d input data matrix (d=1 vector OK)
##         mu - input guess for mean (default is sample mean)
##        sig - d x d input guess for covariance matrix (default is sample)
##          w - input guess for weight (default w=1)
##	sig0 - (do not use)

## Output:
##          list containing estimated  (m=mean sig=sig w=w)

mpdc <- function(x, mu, sig, w, sig0) {
    if (is.matrix(x)) {
        n <- nrow(x)
        d <- ncol(x)
        X <- x
    } else {
        n <- length(x)
        d <- 1
        X <- matrix(x, n, 1)
    }
    onen <- rep(1, n)
    oned <- rep(1, d)
    
    if (missing(mu)) {
        mu <- if (d > 1) colMeans(x) else mean(x)
    } else {
        if (length(mu) != d)  stop("mu wrong length")
    }
    if (missing(w)) {
        w <- 1
    } else {
        if (w <= 0) stop("w must be positive")
    }
    if(missing(sig)) {
        sig <- var(x)
    } else {
        if (length(sig) != d * d) stop("sig wrong size")
    } 
    if (missing(sig0)) {
        sig0 <- diag(d)
    } else {
        if (length(sig0) != d*d) stop("sig0 wrong size")
    }
    dU0 <- diag(chol(sig0))
    ##    assign("dU0",dU0,f=1)	# for determinate of sig0
    
    pdc <- function(x) {
        mu <- x[1:d]
        U <- matrix(0, d, d)
        k <- d * (d + 3) / 2
        U[upper.tri(U, diag=TRUE)] <- x[(d+1):k]
        dU <- exp(diag(U))
        diag(U) <- dU
        w <- exp(x[k+1])
        t1 <- exp( sum( log(dU) + log(dU0) ) )  # cholesky of sig0 and sig-inv
        t2 <- sum(exp( (-.5 * ( (X - outer(onen, mu)) %*% t(U) )^2) %*% oned ))
        w^2 * t1 - 2 * w / n * 2^(d/2) * prod(dU) * t2 + 1
    }
    
    if (d == 1) {
        xU0 <- log(1 / sig0)
    } else {
        tmp <- solve(sig)
        sigi <- (tmp + t(tmp))/2 # make sure sigi is symmetric
        U0 <- chol(sigi)
        diag(U0) <- log(diag(U0))
        xU0 <- U0[upper.tri(U0, diag=TRUE)]
    }
    x0 <- c(mu, xU0, log(w))
    ans <- nlm(pdc, x0, print.level=2, iterlim=100)
    m.ans <- ans$estimate[1:d]
    k <- d * (d + 3) / 2
    w.ans <- exp(ans$estimate[k+1])
    U.ans <- matrix(0, d, d)
    U.ans[row(U.ans)<=col(U.ans)] <- ans$estimate[(d+1):k]
    diag(U.ans) <- exp(diag(U.ans))
    sigi <- crossprod(U.ans)
    sig.ans <- solve(sigi)
    
    list(m=m.ans, sig=sig.ans, w=w.ans)
}

