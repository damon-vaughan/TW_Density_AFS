# ***************************************************************************************************
# **                         SPLUS Code for Furnival's Index of Fit Function                       **
# **                                   Written by Sean M. Garber                                   **
# **                                      December 16, 2001                                        **
# **                                    Revised May 24, 2004                                       **
# **                        (furnival.nls revised by John Moore - August 1, 2009)                  **
# ***************************************************************************************************


furnival <- function(object, ...) UseMethod("furnival")

# This function computes George Furnival's Index of Fit
# A modified likelihood criterion allowing concurrent 
# comparisons of models of different MSE's, weights
# and transformations
# Version 5.24.2004


# Furnival method for the object class "lm"
# corresponding to ordinary least squares linear models
furnival.lm <- function(object, value=1)
{
	SSR		<- sum((residuals(object, type="pearson"))^2)
	n		<- length(object$resid)
	k		<- length(object$coef)
	yobs	<- object$fitted + object$resid
	yhat	<- object$fitted
	degf	<- object$df.residual
	sig		<- sqrt(SSR/degf)		
	FI		<- furncalc(value, sig)
	r2		<- list()
	r2		<- rsquare(yobs, yhat, n, k)
	
	furnout			<- list(FI = FI, rmse = sig, rsg = r2$gen, rsa = r2$adj)
	class(furnout)	<- "furnival.lm"
	furnout
}

# Furnival method for the object class "nls"
# corresponding to ordinary least squares nonlinear models
# changed function to be compatible with nls syntax in R
furnival.nls <- function(object, value=1)
{
	SSR		<- sum(resid(object)^2)
	n		<- length(resid(object))
	k		<- length(summary(object)$parameters[,1])
	yobs	<- (fitted(object) + resid(object)) * value
	yhat	<- fitted(object) * value
	degf	<- n - k
	sig		<- sqrt(SSR/degf)		
	FI		<- furncalc(value, sig)
	r2		<- list()
	r2		<- rsquare(yobs, yhat, n, k)
	
	furnout			<- list(FI = FI, rmse = sig, rsg = r2$gen, rsa = r2$adj)
	class(furnout)	<- "furnival.nls"
	furnout
}

# Furnival method for the object class "gls" and "gnls"
# corresponding to generalized least squares linear and 
# nonlinear models
furnival.gls <- function(object)
{
	n		<- object$dims$N
	k		<- object$dims$p
	yobs	<- object$fitted + object$resid
	yhat	<- object$fitted
	sig		<- object$sigma
	LL		<- object$logLik
	aic		<- AIC(object)
	bic		<- BIC(object)
	r2		<- list()
	r2		<- rsquare(yobs, yhat, n, k)
	
	furnout			<- list(LL = LL, aic = aic, bic = bic, rmse = sig, rsg = r2$gen, rsa = r2$adj)
	class(furnout)	<- "furnival.gls"
	furnout
}

# Furnival method for the object class "lme" and "nlme"
# corresponding to generalized least squares linear and 
# nonlinear mixed effects models
furnival.lme <- function(object)
{	
	nobs					<- object$dims$N
	nparms					<- length(object$coefficients$fixed)
	nlevels				<- object$dims$Q + 1
	levnames				<- labels(object$fitted[1,])
	nranef					<- vector("integer", nlevels)
	ngrps					<- vector("integer", nlevels)
	r2						<- list()
	coefd 					<- matrix(nrow = 2, ncol = nlevels)
	for(i in 1:nlevels)
	{
			nranef[i]		<- object$dims$qvec[nlevels + 1 - i]
			ngrps[i] 		<- object$dims$ngrps[nlevels + 1 - i]
	}
	yobs					<- object$fitted + object$resid	
	for(i in 1:nlevels) 
	{
		yhat				<- object$fitted[,i]
		r2					<- rsquare(yobs, yhat, nobs, nparms, ngrps[i])
		coefd[1, i] 		<- r2$generalized
		coefd[2, i] 		<- r2$adjusted
	}	
	sig						<- object$sigma
	LL						<- object$logLik
	aic						<- AIC(object)
	bic						<- BIC(object)
		
	furnout				<- list(ngrps = ngrps, nlevels = nlevels, levnames = levnames, LL = LL, aic = aic, bic = bic, rmse = sig,
									 coefd = coefd)
	class(furnout) 		<- "furnival.lme"
	furnout
}

furnival.lmer <- function(object)
{	
	nobs<- nrow(model.frame(object))
	nparms<- length(object@fixef)
	nlevels<-  ncol(object@flist)
  levnames<- labels(object@flist)[[2]]
	nranef<- vector("integer", nlevels)
	ngrps	<- vector("integer", nlevels)
	r2<- list()
	coefd<- matrix(nrow = 2, ncol = nlevels)
	for(i in 1:nlevels)
	{
			nranef[i]		<- object$dims$qvec[nlevels + 1 - i]
			ngrps[i] 		<- object$dims$ngrps[nlevels + 1 - i]
	}
	yobs					<- object@y
	for(i in 1:nlevels) 
	{
		yhat				<- fitted(object)
		r2					<- rsquare(yobs, yhat, nobs, nparms, ngrps[i])
		coefd[1, i] 		<- r2$generalized
		coefd[2, i] 		<- r2$adjusted
	}	
	sig						<- object$sigma
	LL						<- object$logLik
	aic						<- AIC(object)
	bic						<- BIC(object)
		
	furnout				<- list(ngrps = ngrps, nlevels = nlevels, levnames = levnames, LL = LL, aic = aic, bic = bic, rmse = sig,
									 coefd = coefd)
	class(furnout) 		<- "furnival.lmer"
	furnout
}


furncalc <- function(invyprime, rmse)
# furncalc(invyprime, rmse) goes through the 4 steps in
# calculating George Furnival's (1961) index of fit
# The subfunction accepts two arguments:
#	invyprime = the inverse of the first derivative of the response function
#		log(Y) = Y, log10(Y) = Ylog(10)
#		weight of W^-2 = W
#	rmse = root mean square error of the model
{
	step1 <- log10(invyprime)
	step2 <- mean(step1)
	step3 <- 10^step2
	step4 <- step3*rmse
	FI    <- step4
	FI
}

rsquare <- function(Y, Ypred, n, k, r=0)
# rsquare(Y, Ypred) calculates Kvalseth's (1985)
# generalized and adjsuted coefficients of determination
# Called by the Furnival() function, but can be
# called separately.
# This subfunction accepts four arguments:
#	Y = observed response data
#	Ypred = fitted values
#	n = number of data points
#	k = number of parameters
{
	rsq <- list()
	SSE <- sum((Y-Ypred)^2)
	SST <- sum((Y-mean(Y))^2)
	a   <- (n-1)/(n-k-r)
	grs <- 1 - (SSE/SST)
	ars <- 1 - a*(SSE/SST)
			
	rsq$generalized <- grs
	rsq$adjusted <- ars
	rsq
}

print.furnival.lm <- function(results, dig=5)
# Print function method for linear ordinary least squares models
# lm objects
{
	cat("\nFit summary stats:\n")
	cat(paste("\nFurnival's Index of Fit:", round(results$FI, digits=dig),
				"\nResidual standard error:", round(results$rmse, digits=dig),
				"\nGeneralized R-square:   ", round(results$rsg, digits=dig),
				"\nAdjusted R-square:      ", round(results$rsa, digits=dig),"\n\n"))
	invisible(results)
}

print.furnival.nls <- function(results, dig=5)
# Print function method for nonlinear ordinary least squares models
# nls objects
{
	cat("\nFit summary stats:\n")
	cat(paste("\nFurnival's Index of Fit:", round(results$FI, digits=dig),
				"\nResidual standard error:", round(results$rmse, digits=dig),
				"\nGeneralized R-square:   ", round(results$rsg, digits=dig),
				"\nAdjusted R-square:      ", round(results$rsa, digits=dig),"\n\n"))
	invisible(results)
}

print.furnival.gls <- function(results, dig=5)
# Print function method for linear and nonlinear generalized
# least squares models gls and glns objects
{
	cat("\nFit summary stats:\n")
	cat(paste("\nResidual standard error:       ", round(results$rmse, digits=dig),
				"\nLog-likelihood:                ", round(results$LL, digits=dig),
				"\nAkaike's information criterion:", round(results$aic, digits=dig),
				"\nBayesian information criterion:", round(results$bic, digits=dig),"\n\n"))
	cat("\nCoefficients of determination:\n")
	cat(paste("\nGeneralized R-square:   ", round(results$rsg, digits=dig),
				"\nAdjusted R-square:      ", round(results$rsa, digits=dig),"\n\n"))
	invisible(results)
}

print.furnival.lme <- function(results, dig=5)
# Print function method for linear and nonlinear mixed effects models
# lme and nlme objects
{
	coefd				<- results$coefd
	cnames				<- results$levnames
	nlevels			<- results$nlevels
	coefd 				<- array(coefd, c(2, nlevels))
	dimnames(coefd) 	<- list(c("Generalized", "Adjusted"), cnames)	
	cat("\nFit summary stats:\n")
	cat(paste("\nResidual standard error:       ", round(results$rmse, digits=dig),
				"\nLog-likelihood:                ", round(results$LL, digits=dig),
				"\nAkaike's information criterion:", round(results$aic, digits=dig),
				"\nBayesian information criterion:", round(results$bic, digits=dig),"\n\n"))
	cat("\nCoefficients of determination:\n")
	print(format(round(coefd, digits = dig)), quote = F)
	cat("\n")
	invisible(results)
}
