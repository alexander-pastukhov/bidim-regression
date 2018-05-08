#' Calculates the bidimensional regression between two 2D configurations
#'
#' @description Calculates the bidimensional regression between two 2D configurations using both Euclidean and Affine transformations following the approach by Tobler (1965).
#' This function assumes strict data format and returns all coefficients and statistics in a single structure. Same functionality is now re-implemented in a R-friendly style, see \code{\link{lm2}} function.
#'
#' @param coord table that must contain two columns for dependent variables (named \code{depV1} and \code{depV2}) and two columns for independent variables (named \code{indepV1} and \code{indepV2}).
#'
#' @return an S3 class \code{BiDimRegression} containing all essential measures of the bidimensional regression
#' * \code{euclidean.r, affine.r} - the regression coefficient, defined analogously to Pearson's r.
#' * \code{euclidean.rsqr, affine.rsqr} - the squared regression coefficient.
#' * \code{euclidean.diABSqr, affine.diABSqr} - the squared distortion index for dependent variables; following Waterman and Gordon's (1984) extension of the bidimensional regression, it provides a measure of comparison of distortions, but the range of values is 0 to 1 following Friedman and Kohler (2003).
#' * \code{euclidean.dMaxABSqr, affine.dMaxABSqr} - the maximal squared distortion index for dependent variables.
#' * \code{euclidean.diXYSqr, affine.diXYSqr} - the distortion index for independent variables.
#' * \code{euclidean.dMaxXYSqr, affine.dMaxXYSqr} - the maximal squared distortion index for independent variables.
#' * \code{euclidean.scaleFactorX, affine.scaleFactorX} - the scaling factor of the first dimension (1.0 means no scaling; values below 1.0 indicate a contraction, values above 1.0 indicate an expansion).
#' * \code{euclidean.scaleFactorY, affine.scaleFactorY} - the scaling factor of the second dimension.
#' * \code{euclidean.angleDEG, affine.angleDEG} - the rotation angle \bold{in degrees}.
#' * \code{euclidean.shear, affine.shear} - shearing of the transformed configuration, always zero for the Euclidean transformation.
#' * \code{euclidean.ttestDF, affine.ttestDF} - degrees of freedom (DF) for the t-tests regarding the model parameters (alphas and betas).
#' * \code{euclidean.alpha1.*, euclidean.alpha2.*, affine.alpha1.*, affine.alpha2.*} - intercept vectors, information includes \code{.coeff} for coefficient, \code{.SE} for standard error, \code{tValue} for _t_-statistics, and \code{pValue} for significance.
#' * \code{euclidean.beta1.*, euclidean.beta2.*, affine.beta1.*, affine.beta2.*, affine.beta3.*, affine.beta4.*} - slope vectors, information includes \code{.coeff} for coefficient, \code{.SE} for standard error, \code{tValue} for _t_-statistics, and \code{pValue} for significance.
#' * \code{euclidean.fValue, affine.fValue} -  F-statistics, following the advice of Nakaya (1997).
#' * \code{euclidean.df1, affine.df1} -  degrees of freedom of the nominator used for the F-statistics propagated by Nakaya (1997); df1 = p-2, with p is the number of elements  needed to calculate the referring model: p=4 for the Euclidean and p=6 for the affine geometry {Nakaya, 1997, Table 1}.
#' * \code{euclidean.df2, affine.df2} -  degrees of freedom of the denominator used for the F-statistics propagated by Nakaya (1997); df2 = 2n-p, with p is the number of elements needed to calculate the referring model (see df1) and n is the number of coordinate pairs.
#' * \code{euclidean.pValue, affine.pValue} -  the significance level based on the preceding F-statistics.
#' * \code{euclidean.dAICso, affine.dAICso} -  the AIC difference between the regarding bidimensional regression model and the bidimensional null model (S0) according to Nakaya (1997), formula 56.
#' * \code{eucVSaff.*} - statistical comparison between Euclidean and Affine models, include \code{.fValue} for F-statistics, \code{.df1} and \code{.df2} for the degrees of freedom, \code{.pValue} for the significance level, and \code{.dAIC} for AIC difference between two models.
#' @md
#'
#' @export
#'
#' @examples
#' resultingMeasures <- BiDimRegression(NakayaData)
#' print(resultingMeasures)
#'
#' @seealso \code{\link{lm2}}
#' @importFrom stats pf lm
BiDimRegression <-
function (coord)
{
	# set standard variables
	n <- dim(coord)[1]   # number of coordinates
	vecZero <- c(rep(0, n))
	vecOne <- c(rep(1, n))
	A <- coord$depV1
	B <- coord$depV2

	X <- coord$indepV1
	Y <- coord$indepV2

	# calculating means
	Am <- mean(A)
	Bm <- mean(B)
	Xm <- mean(X)
	Ym <- mean(Y)

	# calculating (co)variances
	X2 <- sum(X^2)
	Y2 <- sum(Y^2)
	sumX2Y2 <- sum(X^2+Y^2)
	A2 <- sum(A^2)
	B2 <- sum(B^2)
	sumA2B2 <- sum(A^2+B^2)
	varA <- (sum((A-Am)*(A-Am)))/(n)
	varB <- (sum((B-Bm)*(B-Bm)))/(n)
	varX <- (sum((X-Xm)*(X-Xm)))/(n)
	varY <- (sum((Y-Ym)*(Y-Ym)))/(n)
	covABXY <- sum(((A-Am)+(B-Bm))*((X-Xm)+(Y-Ym)))/n

	# ----------- Calculating the Euclidean regression model
	euc_par <- data.frame(
		ps = 4L,
		name = "Euclidean"
		)
	euc_dataMatrix <- matrix(c(vecOne, vecZero, vecZero, vecOne, X, Y, -Y, X), ncol=4)
	euc_target <- matrix(c(A, B),ncol=1)
	euc_data <- data.frame(
			y = (euc_target),
			x1 = (euc_dataMatrix[,1]),
			x2 = (euc_dataMatrix[,2]),
			x3 = (euc_dataMatrix[,3]),
			x4 = (euc_dataMatrix[,4])
		)
	euc_regression <- lm(y ~ 0 + x1 + x2 + x3 +x4, data=euc_data)
	euc_alpha <- data.frame(
			coeff = c(NA, NA),
			SE = c(NA, NA),
			tValue = c(NA, NA),
			pValue = c(NA, NA)
		)
	euc_beta <- data.frame(
			coeff = c(NA, NA),
			SE = c(NA, NA),
			tValue = c(NA, NA),
			pValue = c(NA, NA)
		)

	# retrieving alphas
	for(iAlpha in 1:2) {
		for(iPar in 1:4) {
			euc_alpha[iAlpha, iPar] <- summary(euc_regression)$coeff[iAlpha,iPar]
		}
	}
	# retrieving betas
	for(iBeta in 1:2) {
		for(iPar in 1:4) {
			euc_beta[iBeta , iPar] <- summary(euc_regression)$coeff[iBeta+2,iPar]
		}
	}

	# calculating the scale factors (=sigma)
	euc_scaleFactorX <- sqrt(euc_beta$coeff[1]*euc_beta$coeff[1] + euc_beta$coeff[2]*euc_beta$coeff[2])
	euc_scaleFactorY <- euc_scaleFactorX  # no shear (by definition) for the Euclidean solution
                        # ==> same scaling factors for dimension 1 and dimension 2

	# calculating the rotation (angle) (=theta)
	euc_angleRAD <- atan(euc_beta$coeff[2]/euc_beta$coeff[1])
	euc_angleDEG <- euc_angleRAD*180/pi

	# calculating the shear (gamma)
	euc_shear = 0L # per definition shear must be ZERO within an Euclidean geometry

	# calculating the predicted values
	euc_Apred <- euc_alpha$coeff[1]+euc_beta$coeff[1]*X-euc_beta$coeff[2]*Y
	euc_Bpred <- euc_alpha$coeff[2]+euc_beta$coeff[2]*X+euc_beta$coeff[1]*Y


	euc_Xpred <- (euc_alpha$coeff[2]*euc_beta$coeff[2]-B*euc_beta$coeff[2]-
		euc_alpha$coeff[1]*euc_beta$coeff[1]+A*euc_beta$coeff[1])/
		(euc_beta$coeff[1]*euc_beta$coeff[1]-euc_beta$coeff[2]*euc_beta$coeff[2])

	euc_Ypred <- (euc_beta$coeff[1]*B-euc_alpha$coeff[2]*euc_beta$coeff[1]+
		euc_alpha$coeff[1]*euc_beta$coeff[2]-A*euc_beta$coeff[2])/
		(euc_beta$coeff[1]*euc_beta$coeff[1]+euc_beta$coeff[2]*euc_beta$coeff[2])

	# calculating the bidimensional correlation coefficient
	euc_r <- sqrt(sum(((euc_Apred-Am)*(euc_Apred-Am)) + ((euc_Bpred-Bm)*(euc_Bpred-Bm)))/
	sum(((A-Am)*(A-Am)) + ((B-Bm)*(B-Bm))))


	euc_rsqr <- euc_r*euc_r

	# conducting the inference statistics following Nakaya (1997)
	euc_F <- ((2*n - euc_par$ps)/2)*(euc_rsqr/(1-euc_rsqr))
	euc_df1 <- 2L
	euc_df2 <- 2*n - euc_par$ps 	# set the degrees of freedom: df1/df2
	euc_p <- pf(euc_F, euc_df1, euc_df2, lower.tail = FALSE, log.p = FALSE)



	# ------- Calculating the distortion index following Waterman and Gordon (1984),
	# adjusted by Friedman and Kohler (2003)
	# --- first: calculating distortion index for original configuration
	euc_dDistanceXY <- sqrt(sum((X-euc_Xpred)*(X-euc_Xpred))+sum((Y-euc_Ypred)*(Y-euc_Ypred)))
	euc_dDistanceXYSqr <- euc_dDistanceXY*euc_dDistanceXY
	euc_dMaxXYSqr <- sum((X-Xm)*(X-Xm)+((Y-Ym)*(Y-Ym)))
	euc_dMaxXY <- sqrt(euc_dMaxXYSqr)
	# --- second: calculating distortion index for target configuration
	euc_dDistanceAB <- sqrt(sum((A-euc_Apred)*(A-euc_Apred))+sum((B-euc_Bpred)*(B-euc_Bpred)))
	euc_dDistanceABSqr <- euc_dDistanceAB*euc_dDistanceAB
	euc_dMaxABSqr <- sum((A-Am)*(A-Am)+((B-Bm)*(B-Bm)))  # referring to target configuration
	euc_dMaxAB <- sqrt(euc_dMaxABSqr)
	euc_diABSqr <- euc_dDistanceABSqr/euc_dMaxABSqr
	euc_diAB <- sqrt(euc_diABSqr)

	euc_diXYSqr <- euc_dDistanceABSqr/euc_dMaxXYSqr
	euc_diXY <- sqrt(euc_diXYSqr)

	# ------- Calculation of DAIC (Difference AIC = Akaike Information Criterion)
	#      DAICso: AIC difference DAICso between a bidimensional regression model and
	#              the bidimensional null model
	#              --> if DAICso < 0, the bidimensional regression model is
	#                  better than the bidimensional null model.
	#      [calculation according to Nakaya (1997), formula 56]
	euc_dAICso <- 2L*n*log(1-euc_rsqr)+2L*(euc_par$ps-2L)


	# +++++++ end of Euclidean regression model


	# ----------- Calculating the Affine regression model
      	aff_par <- data.frame(
		ps = 6L,
		name = "Affine"
		)
	aff_dataMatrix <- matrix(c(
		vecOne, vecZero,
		vecZero, vecOne,
		X, vecZero,
		Y, vecZero,
		vecZero, X,
		vecZero, Y), ncol=6)
	aff_target <- matrix(c(A, B),ncol=1)
      	aff_data <- data.frame(
			y = (aff_target),
			x1 = (aff_dataMatrix[,1]),
			x2 = (aff_dataMatrix[,2]),
			x3 = (aff_dataMatrix[,3]),
			x4 = (aff_dataMatrix[,4]),
			x5 = (aff_dataMatrix[,5]),
			x6 = (aff_dataMatrix[,6])
		)
	aff_regression <- lm(y ~ 0 + x1 + x2 + x3 + x4 + x5 + x6, data=aff_data)
	aff_alpha <- data.frame(
			coeff = c(NA, NA),
			SE = c(NA, NA),
			tValue = c(NA, NA),
			pValue = c(NA, NA)
		)
	aff_beta <- data.frame(
			coeff = c(NA, NA),
			SE = c(NA, NA),
			tValue = c(NA, NA),
			pValue = c(NA, NA)
		)

	# retrieving alphas
	for(iAlpha in 1:2) {
		for(iPar in 1:4) {
			aff_alpha[iAlpha, iPar] <- summary(aff_regression)$coeff[iAlpha,iPar]
		}
	}
	# retrieving betas
	for(iBeta in 1:4) {
		for(iPar in 1:4) {
			aff_beta[iBeta , iPar] <- summary(aff_regression)$coeff[iBeta+2,iPar]
		}
	}

	# calculating the rotation (angle) (=theta)
	aff_angleRAD <- atan(aff_beta$coeff[3]/aff_beta$coeff[1])
	aff_angleDEG <- aff_angleRAD*180/pi
	if (aff_beta$coeff[1] < 0)
		{
    		aff_angleDEG <- aff_angleDEG+180
		}

	# calculating the shear (gamma)
	aff_shear <- (((aff_beta$coeff[4]/aff_beta$coeff[2])*sin(aff_angleRAD))+cos(aff_angleRAD))/
				 (((aff_beta$coeff[4]/aff_beta$coeff[2])*cos(aff_angleRAD))-sin(aff_angleRAD))

	# calculating the scale factors (=sigma)
	aff_scaleFactorX <- sqrt(aff_beta$coeff[1]*aff_beta$coeff[1]+aff_beta$coeff[3]*aff_beta$coeff[3])
	if (is.nan(aff_shear))
    		{
		aff_shear <- (aff_beta$coeff[1]-cos(aff_angleRAD)*aff_scaleFactorX)/aff_beta$coeff[3]
		}
	if (is.nan(aff_shear))
    		{
    		aff_shear <- (sin(aff_angleRAD)*aff_scaleFactorX+aff_beta$coeff[2])/aff_beta$coeff[4]
		}
  	aff_scaleFactorY <- aff_beta$coeff[2]/(aff_shear*cos(aff_angleRAD)-sin(aff_angleRAD))
	if (is.nan(aff_scaleFactorY))
    		{
	    	aff_scaleFactorY <- aff_scaleFactorX
		}

	# calculating the predicted values
	aff_Apred <- aff_alpha$coeff[1]+aff_beta$coeff[1]*X+aff_beta$coeff[2]*Y
	aff_Bpred <- aff_alpha$coeff[2]+aff_beta$coeff[3]*X+aff_beta$coeff[4]*Y

	aff_Xpred <- -(aff_Bpred*aff_beta$coeff[2]-aff_Apred*aff_beta$coeff[4]-
		aff_alpha$coeff[2]*aff_beta$coeff[2]+aff_alpha$coeff[1]*aff_beta$coeff[4])/
		(aff_beta$coeff[2]*aff_beta$coeff[3]+aff_beta$coeff[1]*aff_beta$coeff[4])
	aff_Ypred <- -(aff_Apred*aff_beta$coeff[3]-aff_Bpred*aff_beta$coeff[1]-
		aff_alpha$coeff[1]*aff_beta$coeff[3]+aff_alpha$coeff[2]*aff_beta$coeff[1])/
		(aff_beta$coeff[2]*aff_beta$coeff[3]+aff_beta$coeff[1]*aff_beta$coeff[4])

	# calculating the bidimensional correlation coefficient
	aff_r <- sqrt(sum(((aff_Apred-Am)*(aff_Apred-Am))+((aff_Bpred-Bm)*(aff_Bpred-Bm)))/sum(((A-Am)*(A-Am))+((B-Bm)*(B-Bm))))
	aff_rsqr <- aff_r*aff_r

	# conducting the inference statistics according to
	# Nakaya (1997), formula 50
	aff_F <- ((2*n-aff_par$ps)/4) * (aff_rsqr/(1-aff_rsqr))
	aff_df1 <- 4L
	aff_df2 <- 2*n-aff_par$ps;   # set the degrees of freedom: df1/df2
	aff_p <-pf(aff_F, aff_df1, aff_df2, lower.tail = FALSE, log.p = FALSE)

	# ------- Calculating the distortion index following Waterman and Gordon (1984),
	# adjusted by Friedman and Kohler (2003)
	# --- first: calculating distortion index for original configuration
	aff_dDistanceXY <- sqrt(sum((X-aff_Xpred)*(X-aff_Xpred))+sum((Y-aff_Ypred)*(Y-aff_Ypred)))
	euc_dDistanceXY <- sqrt(sum((X-euc_Xpred)*(X-euc_Xpred))+sum((Y-euc_Ypred)*(Y-euc_Ypred)))

	aff_dDistanceXYSqr <- aff_dDistanceXY*aff_dDistanceXY
	euc_dDistanceXYSqr <- euc_dDistanceXY*euc_dDistanceXY

	aff_dMaxXYSqr <- sum((X-Xm)*(X-Xm)+((Y-Ym)*(Y-Ym)))
	euc_dMaxXYSqr <- sum((X-Xm)*(X-Xm)+((Y-Ym)*(Y-Ym)))
	aff_dMaxXY <- sqrt(aff_dMaxXYSqr)
	euc_dMaxXY <- sqrt(euc_dMaxXYSqr)

	# --- second: calculating distortion index for target configuration
	aff_dDistanceAB <- sqrt(sum((A-aff_Apred)*(A-aff_Apred))+sum((B-aff_Bpred)*(B-aff_Bpred)))
	aff_dDistanceABSqr <- aff_dDistanceAB*aff_dDistanceAB
	aff_dMaxABSqr <- sum((A-Am)*(A-Am)+((B-Bm)*(B-Bm)))  # referring to target configuration
	aff_dMaxAB <- sqrt(aff_dMaxABSqr)
	aff_diABSqr <- aff_dDistanceABSqr/aff_dMaxABSqr
	aff_diAB <- sqrt(aff_diABSqr)

	aff_diXYSqr <- aff_dDistanceABSqr/aff_dMaxXYSqr
	aff_diXY <- sqrt(aff_diXYSqr)


	# ------- Calculation of DAIC (Difference AIC = Akaike Information Criterion)
	#      DAICso: AIC difference DAICso between a bidimensional regression model and
	#              the bidimensional null model
	#              --> if DAICso < 0, the bidimensional regression model is
	#                  better than the bidimensional null model.
	#      [calculation according to Nakaya (1997), formula 56]
	aff_dAICso <- 2*n*log(1-aff_rsqr)+2*(aff_par$ps-2)
	euc_dAICso <- 2*n*log(1-euc_rsqr)+2*(euc_par$ps-2)


	#+++++++++++ end of affine solution


	# ---- Calculation of DAICs (Difference AIC = Akaike Information Criterion)
	#      between the different fitted bidimensional regression models
	#      [see Nakaya (1997), table 4]
	# -- comparative test between Euclidean and Affine regression model
	dAICea <- 2*n*log((1-aff_rsqr)/(1-euc_rsqr))+2*(aff_par$ps-euc_par$ps)

	f_ea <- ((2*n-aff_par$ps)/(aff_par$ps-euc_par$ps))*
		((aff_rsqr-euc_rsqr)/(1-aff_rsqr))

	df1_ea <- as.integer(aff_par$ps-euc_par$ps)

	df2_ea <- as.integer(2*n-aff_par$ps)
	p_ea <- pf(f_ea, df1_ea, df2_ea, lower.tail = FALSE, log.p = FALSE)

	# df of parameter t-tests
	tTestDF <- aff_regression$df[1]

	# return all the results to a data.frame
	res_euc <- data.frame(r=euc_r, rsqr=euc_rsqr, diABSqr=euc_diABSqr, dMaxABSqr=euc_dMaxABSqr, diXYSqr=euc_diXYSqr, dMaxXYSqr=euc_dMaxXYSqr, 	scaleFactorX=euc_scaleFactorX, scaleFactorY=euc_scaleFactorY, angleDEG=euc_angleDEG, shear=euc_shear, ttestDF=tTestDF, alpha1=euc_alpha[1,], 	alpha2=euc_alpha[2,], beta1=euc_beta[1,], beta2=euc_beta[2,], beta3=euc_beta[3,], beta4=euc_beta[4,], fValue=euc_F, df1=euc_df1, df2=euc_df2, 	pValue=euc_p, dAICso=euc_dAICso)
	res_aff <- data.frame(r=aff_r, rsqr=aff_rsqr, diABSqr=aff_diABSqr, dMaxABSqr=aff_dMaxABSqr, diXYSqr=aff_diXYSqr, dMaxXYSqr=aff_dMaxXYSqr, 	scaleFactorX=aff_scaleFactorX, scaleFactorY=aff_scaleFactorY, angleDEG=aff_angleDEG, shear=aff_shear, ttestDF=tTestDF, alpha1=aff_alpha[1,], 	alpha2=aff_alpha[2,], beta1=aff_beta[1,], beta2=aff_beta[2,], beta3=aff_beta[3,], beta4=aff_beta[4,], fValue=aff_F, df1=aff_df1, df2=aff_df2, 	pValue=aff_p, dAICso=aff_dAICso)
	euclideanVSaffine <- data.frame(dAIC=dAICea, fValue=f_ea, pValue=p_ea, df1=df1_ea, df2=df2_ea)

	results_sum <- data.frame(euclidean=res_euc, affine=res_aff, eucVSaff=euclideanVSaffine)
	class(results_sum) <- "BiDimRegression"

	invisible(results_sum)   # returns the measures of the bidimensional regression
}

#' @export
#' @importFrom stats printCoefmat
print.BiDimRegression <-
function(x, ...)
{
  # version of the routine
  verRoutine <- "1.0.6"

  # ---- OUTPUT the data ----
  for(iOutputLoop in 0:1) {
    if (iOutputLoop == 1) {
      # declaration of the output file
      sink(file="BiDimRegrOutput.txt", type="output", append=TRUE)
    }

    # print the results of the BiDimensional Regression Analysis
    cat(sprintf('\n-------------------\n'))
    cat(sprintf('BiDimRegression %s \n', verRoutine))
    cat(sprintf('Date-Time: %s\n', date()))
    cat(sprintf('----------------------------------------------------------------------\n'))

    euc_alpha=matrix(c(x$euclidean.alpha1.coeff, x$euclidean.alpha1.SE, x$euclidean.alpha1.tValue,
                       x$euclidean.ttestDF, x$euclidean.alpha1.pValue,
                       x$euclidean.alpha2.coeff, x$euclidean.alpha2.SE, x$euclidean.alpha2.tValue,
                       x$euclidean.ttestDF, x$euclidean.alpha2.pValue),
                     nrow=2, byrow = TRUE)
    euc_beta=matrix(c(x$euclidean.beta1.coeff, x$euclidean.beta1.SE,
                      x$euclidean.beta1.tValue, x$euclidean.ttestDF, x$euclidean.beta1.pValue,
                      x$euclidean.beta2.coeff, x$euclidean.beta2.SE, x$euclidean.beta2.tValue,
                      x$euclidean.ttestDF, x$euclidean.beta2.pValue),
                    nrow=2, byrow = TRUE)
    aff_alpha=matrix(c(x$affine.alpha1.coeff, x$affine.alpha1.SE, x$affine.alpha1.tValue,
                       x$affine.ttestDF, x$affine.alpha1.pValue,
                       x$affine.alpha2.coeff, x$affine.alpha2.SE, x$affine.alpha2.tValue,
                       x$affine.ttestDF, x$affine.alpha2.pValue),
                     nrow=2, byrow = TRUE)
    aff_beta=matrix(c(x$affine.beta1.coeff, x$affine.beta1.SE, x$affine.beta1.tValue,
                      x$affine.ttestDF, x$affine.beta1.pValue,
                      x$affine.beta2.coeff, x$affine.beta2.SE, x$affine.beta2.tValue, x$affine.ttestDF, x$affine.beta2.pValue,
                      x$affine.beta3.coeff, x$affine.beta3.SE, x$affine.beta3.tValue, x$affine.ttestDF, x$affine.beta3.pValue,
                      x$affine.beta4.coeff, x$affine.beta4.SE, x$affine.beta4.tValue, x$affine.ttestDF, x$affine.beta4.pValue),
                    nrow=4, byrow = TRUE)

    # -- Overall analysis
    overallStats=matrix(c(
      x$euclidean.r, x$euclidean.rsqr, x$euclidean.fValue, x$euclidean.df1, x$euclidean.df2, x$euclidean.pValue,
      x$affine.r, x$affine.rsqr, x$affine.fValue, x$affine.df1, x$affine.df2, x$affine.pValue),
      nrow=2, byrow = TRUE)
    colnames(overallStats) <- c("r", "r-sqr", "F-value", "df1", "df2", "p-value")
    rownames(overallStats) <- c("Euclidean", "Affine")

    cat(sprintf('\n--- overall statistics ---\n'))
    printCoefmat(overallStats, P.values=TRUE, has.Pvalue=TRUE, digits=3, tst.ind=4, signif.stars=TRUE)

    # show warning for unidentified model
    if (x$euclidean.df2 < x$euclidean.df1)	{
      cat(sprintf('WARNING: Euclidean model is not defined'))
    }
    if (x$affine.df2 < x$affine.df1) {
      cat(sprintf('\t\t\t\t\tWARNING: Affine model is not defined\n'))
    }

    cat(sprintf('----------------------------------------------------------------------\n'))
    cat(sprintf('\n--- parameters ---\n'))
    cat(sprintf('\n- Euclidean -\n'))

    # -- Parameter analyses
    colNames <- c("Parameter", "Std.Err", "t-value", "df", "p-value")
    rowNamesAlphas <- c("alpha1", "alpha2")
    rowNames2Betas <- c("beta1 ", "beta2 ")
    rowNames4Betas <- c("beta1 ", "beta2 ", "beta3 ", "beta4 ")
    colnames(euc_alpha) <- colNames
    colnames(euc_beta) <- colNames
    colnames(aff_alpha) <- colNames
    colnames(aff_beta) <- colNames
    rownames(euc_alpha) <- rowNamesAlphas
    rownames(euc_beta) <- rowNames2Betas
    rownames(aff_alpha) <- rowNamesAlphas
    rownames(aff_beta) <- rowNames4Betas

    printCoefmat(euc_alpha, P.values = TRUE, has.Pvalue = TRUE, digits=3, tst.ind=4)
    cat(sprintf('\n'))
    printCoefmat(euc_beta, P.values = TRUE, has.Pvalue = TRUE, digits=3, tst.ind=4)
    cat(sprintf('\n'))

    cat(sprintf('\n- Affine -\n'))
    printCoefmat(aff_alpha, P.values = TRUE, has.Pvalue = TRUE, digits=3, tst.ind=4)
    cat(sprintf('\n'))
    printCoefmat(aff_beta, P.values = TRUE, has.Pvalue = TRUE, digits=3, tst.ind=4)
    cat(sprintf('\n'))

    cat(sprintf('----------------------------------------------------------------------\n'))
    cat(sprintf('\n--- details ---\n'))
    cat(sprintf('\n- Euclidean -\t\t\t- Affine -\n'))
    cat(sprintf('scaleX\t= scaleY = %4.3f\tscaleX\t= %4.3f, scaleY = %4.3f\n',
                x$euclidean.scaleFactorX, x$affine.scaleFactorX, x$affine.scaleFactorY))
    cat(sprintf('shear\t= %4.3f\t\t\tshear\t= %4.3f\n',
                x$euclidean.shear, x$affine.shear))
    cat(sprintf('angle\t= %4.3f DEG\t\tangle\t= %4.3f DEG\n',
                x$euclidean.angleDEG, x$affine.angleDEG))

    cat(sprintf('---\t\t\t\t---\n'))
    cat(sprintf('DAIC (agst.0)\t= %4.2f\tDAIC (agst.0)\t= %4.2f\n',
                x$euclidean.dAICso, x$affine.dAICso))

    cat(sprintf('---\t\t\t\t---\n'))
    cat(sprintf('dMaxABSqr\t= %4.3f\tdMaxABSqr\t= %4.3f\n',
                x$euclidean.dMaxABSqr, x$affine.dMaxABSqr))
    cat(sprintf('diABSqr  \t= %4.3f\t\tdiABSqr  \t= %4.3f\n',
                x$euclidean.diABSqr, x$affine.diABSqr))
    cat(sprintf('dMaxXYSqr\t= %4.3f\tdMaxXYSqr\t= %4.3f\n',
                x$euclidean.dMaxXYSqr, x$affine.dMaxXYSqr))
    cat(sprintf('diXYSqr  \t= %4.3f\t\tdiXYSqr  \t= %4.3f\n',
                x$euclidean.diXYSqr, x$affine.diXYSqr))
    cat(sprintf('----------------------------------------------------------------------\n'))
    cat(sprintf('\n--- comparative statistics of fitted bidimensional regxsion models\n\n'))

    comparativeStats=matrix(c(
      x$eucVSaff.fValue, as.integer(x$eucVSaff.df1), as.integer(x$eucVSaff.df2), x$eucVSaff.pValue),
      nrow=1, byrow = TRUE)
    colnames(comparativeStats) <- c("F-value", "df1", "df2", "p-value")
    rownames(comparativeStats) <- c("Euclidean vs. Affine")
    printCoefmat(comparativeStats, P.values=TRUE, has.Pvalue=TRUE, digits=3, cs.ind=1, signif.stars=TRUE)

    cat(sprintf('\n'))

    if (x$eucVSaff.df2 < x$eucVSaff.df1) {
      cat(sprintf('WARNING: model is not defined\n'))
    }
    if (x$eucVSaff.pValue <= .05) {
      if (x$eucVSaff.dAIC<0) {
        superiorSolution <- '(significantly better: Affine solution)'
      } else {
        superiorSolution <- '(significantly better: Euclidean solution)'
      }
    } else {
      superiorSolution = '(not significantly different solutions)'
    }

    cat(sprintf('\nDAICea = %4.3f %s\n\n', x$eucVSaff.dAIC, superiorSolution))
    cat(sprintf('**********************************************************************\n\n\n'))
  }

  sink()   # close the output result file
}


#' @export
summary.BiDimRegression <-
function(object, ...)
{
  summary.BiDimRegression <- object
  class(summary.BiDimRegression) <- "summary.BiDimRegression"
  return(summary.BiDimRegression)
}

#' @export
print.summary.BiDimRegression <-
function(x, ...)
{
  # -- short version of the statistics
  # Overall analysis
  overallStats=matrix(c(
    x$euclidean.r, x$euclidean.rsqr, x$euclidean.fValue, x$euclidean.df1, x$euclidean.df2, x$euclidean.pValue,
    x$affine.r, x$affine.rsqr, x$affine.fValue, x$affine.df1, x$affine.df2, x$affine.pValue),
    nrow=2, byrow = TRUE)
  colnames(overallStats) <- c("r", "r-sqr", "F-value", "df1", "df2", "p-value")
  rownames(overallStats) <- c("Euclidean", "Affine")

  cat(sprintf('\n--- summary statistics from bidimensional regressions ---\n'))
  printCoefmat(overallStats, P.values=TRUE, has.Pvalue=TRUE, digits=3, tst.ind=4, signif.stars=TRUE)

}
