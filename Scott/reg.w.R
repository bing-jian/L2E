#	R version June 16, 2004

# Copyright David W. Scott     February 5, 2004   

# Partial L2E multivariate linear regression

# Purpose:  fit model   y = a + t(b)*x + e  by L2E

#	the error distribution can be either  N(0,sig^2)  or  w N(0,sig^2)

# Inputs:
#          X - n x p input data matrix (p=1 vector OK)
#          y - n responses
#        xin - input guesses for (a,b,sig) length 1+p+1
#		defaults = c(mean(y),0,..,0,sd(y))
#      w.opt - False (default e ~ N(0,sig^2) or  True (e ~ w N(0,sig^2))
#        win - fixed value for w   or   initial guess for w (see w.opt)
#         pl - print level for nlmin (1=default  0=none)

# Output:      list containing estimated coefficients
#         $a - intercept
#         $b - vector of estimated beta's
#       $sig - estimated standard deviation of residuals
#         $w - estimated weight (note:  w can be greater than 1)
#       $res - residual vector


reg.w.mat <- function(X,y,xin,w.opt=F,win=1,pl=1) {

    if(!is.matrix(X)) {X <- cbind(X)}; p <- ncol(X)
    if(length(y)!=nrow(X)) stop("X and y differ on n")

    reg.w.mat.crit <- function(x) {
	a <- x[1]; b <- x[2:(p+1)]; sig <- exp(x[p+2])
	if(w.opt) { w <- exp(x[p+3]) } else { w <- win }
if(br>1) {br <<- br-1} else {browser()}
	ei <- y - a - X%*%b
	w^2/(2*sqrt(pi)*sig) - 2*w*mean(dnorm(ei,0,sig))
    }

    if(missing(xin)) { xin <- c(mean(y),rep(0,p),sqrt(var(y))) }
    if(length(xin)!=p+2) stop("xin error")
    x0 <- c(xin[1:(p+1)],log(xin[p+2]))
    if(w.opt) { x0 <- c(x0,log(win)) }

    ans <- nlm(reg.w.mat.crit,x0,iterlim=100,print.level=pl)

    a <- ans$est[1]; b <- ans$est[2:(p+1)]; sig <- exp(ans$est[p+2])
    if(w.opt) {w <- exp(ans$est[p+3])} else {w <- win}

    invisible(list(a=a,b=b,sig=sig,w=w,res=c(y-a-X%*%b)))
}
