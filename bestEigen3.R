require(Matrix)
bestEigen <- function(x=x,a1=NULL,a2=NULL,tol=1.0e-6,maxiter=100) {
#	set.seed(12345)
	n <- nrow(x)
	s1 <- colSums(x)/n
	s2 <- colSums(x^2)/n
	z <- s2-s1*s1
#
	bad <- which(z==0)
	good <- setdiff(1:n,bad)
	x <- x[good,good]
	for (i in 1:3) junk <- gc()
	k <- nrow(x)
	s1 <- colSums(x)/k
	s2 <- colSums(x^2)/k
	sd <- sqrt(s2-s1*s1)
	x <- x %*% Diagonal(x=1/sd)
	C <- s1/sd
	if (is.null(a1)) a1 <- rnorm(nrow(x))
	if (is.null(a2)) a2 <- rnorm(nrow(x))
	xx <- t(x)
	for (i in 1:3) junk <- gc()

	niter <- 0

	while (niter < maxiter) {
		for (i in 1:10) {
			a1 <- x %*% a1 - sum(C*a1)
			a1 <- xx %*% a1 - sum(a1)*C
			a1 <- a1/sqrt(sum(a1*a1))
		}
		ev1 <- a1
		a1 <- x %*% a1 - sum(C*a1)
		a1 <- xx %*% a1 - sum(a1)*C
		lam1 <- sqrt(sum(a1*a1))
		er1 <- sqrt(sum((a1 - lam1*ev1)^2))
		a1 <- a1/sqrt(sum(a1*a1))

		a2 <- a2 - a1*sum(a2*a1)

		for (i in 1:10) {
			a2 <- x %*% a2 - sum(C*a2)
			a2 <- xx %*% a2 - sum(a2)*C
			a2 <- a2 - a1*sum(a2*a1)
			a2 <- a2/sqrt(sum(a2*a2))
		}

		ev2 <- a2
		a2 <- x %*% a2 - sum(C*a2)
		a2 <- xx %*% a2 - sum(a2)*C
		lam2 <- sqrt(sum(a2*a2))
		a2 <- a2/sqrt(sum(a2*a2))
		
		if (lam1 < lam2) {
			temp <- a1
			a1 <- a2
			a2 <- temp
			temp <- lam1
			lam1 <- lam2
			lam2 <- temp
			for (junk in 1:3) gc()
		}

		niter <- niter + 10
		if (er1/(lam1-lam2) < tol) break
	}

	v <- rep(NA,n)
	v[good] <- ev1
	return(list(v=v,lam1=lam1,er=er1,lam2=lam2,niter=niter,a1=a1,a2=a2))
}
