library(Formula)


#' Fitting Bidimensional Regression Models
#'
#' lm2 is used to fit bidimensional linear regression models using
#' Euclidean and Affine transformations following the approach by Tobler (1965).
#'
#' @usage
#' lm2(formula, data, transformation)
#'
#' @param formula a symbolic description of the model to be fitted in the format \code{A + B ~ C + D}, where
#' \code{A} and \code{B} are dependent and \code{C} and \code{D} are independent variables
#' @param data a data frame containing variables for the model.
#' @param transformation the transformation to be used, either \code{'euclidean'}, \code{'affine'}, or \code{'projective'}.
#'
#' @return lm2 returns an object of class "lm2".
#' An object of class "lm" is a list containing at least the following components:
#' \item{\code{transformation}}{string with the transformation type (\code{euclidean}, \code{affine}, or \code{projective})}
#' \item{\code{npredictors}}{number of predictors used in the model: 4 for euclidean, 6 for affine, 8 for projective.}
#' \item{\code{df_model, df_residual}}{degrees of freedom for the model and for the residuals}
#' \item{\code{transformation_matrix}}{\code{3x3} transformation matrix}
#' \item{\code{coeff}}{transformation coefficients, with \code{a} denoting the intercept terms.}
#' \item{\code{transformed_coeff}}{\code{scale}, \code{angle}, and \code{sheer} coefficients, depends on transformation.}
#' \item{\code{fitted_values}}{data frame containing fitted values for the original data set}
#' \item{\code{residuals}}{data frame containing residuals  for the original fit}
#' \item{\code{r.squared, adj.r.squared}}{R-squared and adjusted R-squared.}
#' \item{\code{F, p.value}}{F-statistics and the corresponding p-value, given the \code{df_model} and \code{df_residual} degrees of freedom.}
#' \item{\code{dAIC}}{Akaike Information Criterion (AIC) difference between the regression model and the null model. A negative values indicates that the regression model is better. See \cite{Nakaya (1997)}.}
#' \item{\code{distortion_index}}{Distortion index following \cite{Waterman and Gordon (1984)}, as adjusted by \cite{Friedman and Kohler (2003)}}
#' \item{\code{lm}}{an underlying \link[=lm]{linear model} for \code{Euclidean} and \code{affine} transformations.}
#' \item{\code{formula}}{formula, describing input and output columns}
#' \item{\code{data}}{data used to fit the model}
#' \item{\code{Call}}{function call information, incorporates the \code{formula}, \code{transformation}, and \code{data}.}
#' @export
#' @seealso \code{\link{anova.lm2}} \code{\link{BiDimRegression}}
#'
#' @examples
#' lm2euc <- lm2(depV1 + depV2 ~ indepV1 + indepV2, NakayaData, 'euclidean')
#' lm2aff <- lm2(depV1 + depV2 ~ indepV1 + indepV2, NakayaData, 'affine')
#' lm2prj <- lm2(depV1 + depV2 ~ indepV1 + indepV2, NakayaData, 'projective')
#' anova(lm2euc, lm2aff, lm2prj)
#' predict(lm2euc)
#' summary(lm2euc)
lm2 <- function(formula, data, transformation) { UseMethod("lm2") }

#' @export
#' @importFrom methods is
lm2.formula <-  function(formula, data, transformation){

  # Check arguments ---------------------------------------------------------
  # Are they present?
  if(missing(formula))
  {
    stop("'formula' is missing or incorrect")
  }
  if (missing(data)){
    stop('argument "data" is missing')
  }
  if (missing(transformation)){
    stop('argument "transformation" is missing')
  }

  # Valid type and values?
  if (!is.data.frame(data))
  {
    stop('argument "data" must be a data frame')
  }
  if (!is.character(transformation) || !(tolower(transformation) %in% c('euclidean', 'affine', 'projective'))){
    stop("unknown transformation, please use either 'euclidean', 'affine', or 'projective'")
  }


  # Extract variables from dataframe ----------------------------------------
  model_formula <- Formula::Formula(formula)
  DV <- Formula::model.part(model_formula, data = data, lhs = 1)
  if (any(!sapply(DV, is.numeric))){
    stop('Non-numeric dependent varaible')
  }
  IV <- Formula::model.part(model_formula, data = data, rhs = 1)
  if (any(!sapply(IV, is.numeric))){
    stop('Non-numeric independent varaible')
  }

  # Fit the model -----------------------------------------------------------
  lm2model <- lm2fit(cbind(DV, IV), tolower(transformation))




  # Adding information about the call ---------------------------------------
  lm2model$Call <- match.call(expand.dots = FALSE)
  lm2model$formula <- formula
  m <- match(c("formula", "data", "transformation"), names(lm2model$Call), 0L)
  lm2model$formula <- lm2model$Call[2]
  lm2model$data <- lm2model$Call[3]

  return(lm2model)
}


#' Fits the specified model and computes stats
#'
#' Calls a specific transformation model function and then computes statistics
#' that is common across all transformations.
#' This function should not be called directly, please use \code{\link{lm2}}.
#'
#' @param data the preprocessed data frame from \code{\link{lm2}} function,
#' so that the first two columns are the dependent variables and the other
#' two are independent variables
#' @param transformation the transformation to be used, either \code{'euclidean'} or \code{'affine'}.
#'
#' @return returns an object of class "lm2", see \code{\link{lm2}}
#' for the description.
#'
#' @keywords internal
#' @importFrom stats pf
lm2fit <- function(data, transformation){
  lm2model <- switch(transformation,
                     euclidean = lm2euclidean(data),
                     affine = lm2affine(data),
                     projective= lm2projective(data),
                     stop("unknown transformation, please use either 'euclidean' or 'affine'"))
  class(lm2model) <- 'lm2'


  # common stats for the bidimensional regresion
  var_mean <- colMeans(data)
  n <- nrow(data)
  lm2model$r.squared <- 1- sum((lm2model$fitted_values[, 1]-data[, 1])^2 +(lm2model$fitted_values[, 2]-data[, 2])^2)/
    sum((data[, 1]-var_mean[[1]])^2 +(data[, 2]-var_mean[[2]])^2)
  lm2model$adj.r.squared <- 1-( ( (n-1)/(n-lm2model$npredictors-1)) * ( (n-2)/(n-lm2model$npredictors-2)) * ((n+1)/n))*(1-lm2model$r.squared)
  lm2model$dAIC<- 2*n*log(1-lm2model$r.squared)+2*lm2model$df_model
  lm2model$F <- (lm2model$df_residual/lm2model$df_model)*(lm2model$r.squared/(1-lm2model$r.squared))
  lm2model$p.value<- pf(lm2model$F, lm2model$df_model, lm2model$df_residual, lower.tail= FALSE, log.p= FALSE)

  ## ------- the distortion index following Waterman and Gordon (1984), adjusted by Friedman and Kohler (2003)
  di<- data.frame(D.sqr= c(NA,NA), Dmax.sqr= c(NA,NA), DI.sqr= c(NA,NA), row.names = c('Dependent', 'Independent'))
  di$D.sqr[1]<- sum((data[, 1]-lm2model$fitted_values[, 1])^2)+
    sum((data[, 2]-lm2model$fitted_values[, 2])^2)
  di$D.sqr[2]<- sum((data[, 3]-lm2model$fitted.I[, 3])^2)+
    sum((data[, 4]-lm2model$fitted.I[, 4])^2)

  di$Dmax.sqr[1] <- sum((data[, 1]-var_mean[[1]])^2 +(data[, 2]-var_mean[[2]])^2)
  di$Dmax.sqr[2] <- sum((data[, 3]-var_mean[[3]])^2 +(data[, 4]-var_mean[[4]])^2)
  di$DI.sqr <- di$D.sqr/di$Dmax.sqr
  lm2model$distortion_index <- di

  return(lm2model)
}


# Euclidean ---------------------------------------------------------------


#' Computes model for the euclidean transformation
#'
#' @param data the preprocessed data frame from \code{\link{lm2}} function,
#' so that the first two columns are the dependent variables and the other
#' two are independent variables
#'
#' @return object with transformation specific data to be supplemented with further stats
#' @keywords internal
#' @importFrom stats lm predict setNames
lm2euclidean <- function(data){

  lm2model <- list(transformation= 'euclidean',
                   npredictors= 4,
                   df_model= 2L,
                   df_residual= 2*nrow(data)-4L)

  # arranging the data frame for the lm function
  cZeros <- c(rep(0, nrow(data)))
  cOnes <- c(rep(1, nrow(data)))
  lm_data <- data.frame(
    y= c(data[, 1], data[, 2]),
    a1= c(cOnes, cZeros),
    a2= c(cZeros, cOnes),
    b1 = c( data[, 3], data[, 4]),
    b2 = c(-data[, 4], data[, 3]))

  # using lm to fit the model
  lm2model$lm <- stats::lm(y ~ 0 + a1 + a2 + b1 +b2, data= lm_data)

  # coefficients and the transformation matrix
  lm2model$coeff <- summary(lm2model$lm)$coeff[, 1]
  lm2model$transformation_matrix <- matrix(c(lm2model$coeff['b1'], -lm2model$coeff['b2'], lm2model$coeff['a1'],
                                           lm2model$coeff['b2'],  lm2model$coeff['b1'], lm2model$coeff['a2'],
                                           0,0,1), nrow=3)
  # calculating the transformed coefficients
  lm2model$transformed_coeff <- c(
    sqrt(lm2model$coeff[['b1']]^2 + lm2model$coeff[['b2']]^2),
    sqrt(lm2model$coeff[['b1']]^2 + lm2model$coeff[['b2']]^2),
    atan2(lm2model$coeff[['b2']], lm2model$coeff[['b1']])
  )
  names(lm2model$transformed_coeff) <- c('scale1', 'scale2', 'angle')

  # getting the predicted values for dependent variables
  lm2model$fitted_values <- setNames(data.frame(matrix(predict(lm2model$lm), ncol=2)), colnames(data)[1:2])


  # getting the residuals
  lm2model$residuals <- setNames(data.frame(matrix(lm2model$lm$residuals, ncol=2)), colnames(data)[1:2])

  return(lm2model)
}



# Affine ------------------------------------------------------------------


#' Computes model for the affine transformation
#'
#' @param data the preprocessed data frame from \code{\link{lm2}} function,
#' so that the first two columns are the dependent variables and the other
#' two are independent variables
#'
#' @return object with transformation specific data to be supplemented with further stats
#' @keywords internal
#' @importFrom stats lm predict setNames
lm2affine <- function(data){
  lm2model <- list(transformation= 'affine',
                   npredictors= 6,
                   df_model= 4L,
                   df_residual= 2*nrow(data)-6L)

  # re-arraging data for affine regression model
  cZeros <- c(rep(0, nrow(data)))
  cOnes <- c(rep(1, nrow(data)))
  lm_data <- data.frame(
    y= c(data[, 1], data[, 2]),
    a1= c(cOnes, cZeros),
    a2= c(cZeros, cOnes),
    b1= c(data[, 3], cZeros),
    b2= c(data[, 4], cZeros),
    b3= c(cZeros, data[, 3]),
    b4= c(cZeros, data[, 4]))

  # using lm to fit the model
  lm2model$lm <- stats::lm(y ~ 0 + a1 + a2 + b1 +b2 +b3 +b4, data= lm_data)

  # coefficients and the transformation matrix
  lm2model$coeff <- summary(lm2model$lm)$coeff[, 1]
  lm2model$transformation_matrix <- matrix(c(lm2model$coeff['b1'],  lm2model$coeff['b2'], lm2model$coeff['a1'],
                                             lm2model$coeff['b3'],  lm2model$coeff['b4'], lm2model$coeff['a2'],
                                             0,0,1), nrow=3)

  # Calculating the transformed coefficients
  aff_angle <- atan2(lm2model$coeff[['b3']], lm2model$coeff[['b1']])
  aff_shear <- ((lm2model$coeff[['b4']]/lm2model$coeff[['b2']])*sin(aff_angle)+cos(aff_angle))/
    ((lm2model$coeff[['b4']]/lm2model$coeff[['b2']])*cos(aff_angle)-sin(aff_angle))

  aff_scale1 <- sqrt(lm2model$coeff[['b1']]^2+lm2model$coeff[['b3']]^2)
  if (is.nan(aff_shear))
  {
    aff_shear <- (lm2model$coeff[['b1']]-cos(aff_angle)*aff_scale1)/lm2model$coeff[['b3']]
  }
  if (is.nan(aff_shear))
  {
    aff_shear <- (sin(aff_angle)*aff_scale1+lm2model$coeff[['b2']])/lm2model$coeff[['b4']]
  }

  aff_scale2 <- lm2model$coeff[['b2']]/(aff_shear*cos(aff_angle)-sin(aff_angle))
  if (is.nan(aff_scale2))
  {
    aff_scale2 <- aff_scale1
  }
  lm2model$transformed_coefficients <- c(
    aff_scale1, aff_scale2, aff_shear, aff_angle
  )
  names(lm2model$transformed_coefficients) <- c('scale1', 'scale2', 'shear', 'angle')


  # getting the predicted values for dependent variables
  lm2model$fitted_values <- setNames(data.frame(matrix(predict(lm2model$lm), ncol=2)), colnames(data)[1:2])


  # getting the residuals
  lm2model$residuals <- setNames(data.frame(matrix(lm2model$lm$residuals, ncol=2)), colnames(data)[1:2])

  return(lm2model)
}


# Projective --------------------------------------------------------------

#' Computes model for the projective transformation
#'
#' @param data the preprocessed data frame from \code{\link{lm2}} function,
#' so that the first two columns are the dependent variables and the other
#' two are independent variables
#'
#' @return object with transformation specific data to be supplemented with further stats
#' @keywords internal
#' @importFrom stats setNames
lm2projective <- function(data){
  ## Preparing a placeholder for the class
  lm2model <- list(transformation= 'projective',
                   npredictors= 8,
                   df_model= 6L,
                   df_residual= 2*nrow(data)-8L)

  # preparing  matrix for the projective  model
  A <- matrix(0, nrow= nrow(data) * 2, ncol = 9)
  cZeros <- rep(0, nrow(data))
  cOnes <- rep(1, nrow(data))

  # laying out parameters
  A[, 1] <- c(data[, 3], cZeros)
  A[, 2] <- c(data[, 4], cZeros)
  A[, 3] <- c(cOnes,    cZeros)

  A[, 4] <- c(cZeros, data[, 3])
  A[, 5] <- c(cZeros, data[, 4])
  A[, 6] <- c(cZeros, cOnes)

  A[, 7] <- c( -data[, 1]*data[, 3], -data[, 2]*data[, 3])
  A[, 8] <- c( -data[, 1]*data[, 4], -data[, 2]*data[, 4])
  A[, 9] <- c(  data[, 1],            data[, 2])

  # using singular value decomposition
  V <- svd(A)$v

  # copying over results
  lm2model$coeff <- -V[1:8, 9]/V[9,9]
  names(lm2model$coeff) <- c('a1', 'a2', 'a3', 'b1', 'b1', 'b3', 'c1', 'c2')
  lm2model$transformation_matrix <- matrix(c(lm2model$coeff, 1), nrow=3)

  # computing the predicted values for dependent variables
  IV <- data[, 3:4]
  IV$z <- 1
  fitted_DV <- data.matrix(IV) %*% lm2model$transformation_matrix
  lm2model$fitted_values <- data.frame(fitted_DV[, 1:2]/fitted_DV[, 3])
  colnames(lm2model$fitted_values) <- colnames(data)[1:2]

  # getting the residuals
  lm2model$residuals <- setNames(data[, 1:2]-lm2model$fitted_values, colnames(data)[1:2])

  return(lm2model)
}


# Printing and summary ----------------------------------------------------

#' @export
print.lm2 <- function(x, ...){
  cat(sprintf('Call:\n'))
  cat(deparse(x$Call))
  cat('\n\n')
  cat('Coefficients:\n')
  coeff <- data.frame(as.list(x$coeff))
  rownames(coeff) <- ''
  printCoefmat(coeff)

  # transformed coefficients, if applicable
  if ("transformed_coeff" %in% names(x)){
    cat('\nTransformed coefficients:\n')
    transformed_coeff <- data.frame(as.list(x$transformed_coeff))
    rownames(transformed_coeff) <- ''
    printCoefmat(transformed_coeff)
  }

  # correlation strength
  cat('\nMultiple R-squared:', x$r.squared, '\tAdjusted R-squared:', x$adj.r.squared)
}


#' Makes a lightweight summary lm2 object
#'
#' Drops heavy bits, like the data frame with predicted values or the lm object.
#' However, the print tells more! :)
#'
#' @param object an object of class "lm2", see \code{\link{lm2}}
#'
#' @export
#' @keywords internal
summary.lm2 <- function(object, ...){
  # copying most of the object
  object_summary<- object
  class(object_summary) <- "summary.lm2"

  # dropping heavy bits
  if ('lm' %in% names(object_summary)){
    object_summary$coeff <- summary(object_summary$lm)$coeff
  }
  else{
    object_summary$coeff <- data.frame(as.list(object$coeff))
    rownames(object_summary$coeff) <- ''
  }
  object_summary$fitted_values <- NULL
  object_summary$lm <- NULL

  return(object_summary)
}

#' @export
print.summary.lm2 <- function(x, ...){
  cat(sprintf('Call:\n'))
  cat(deparse(x$Call))
  cat('\n\n')
  cat('Coefficients:\n')
  printCoefmat(x$coeff)

  # transformed coefficients
  if ('transformed_coeff' %in% names(x)){
    cat('\nTransformed coefficients:\n')
    transformed_coeff <- data.frame(as.list(x$transformed_coeff))
    rownames(transformed_coeff) <- ''
    printCoefmat(transformed_coeff)
  }

  # distortion index
  cat('\nDistortion index:\n')
  di <- t(x$distortion_index)
  rownames(di)<- c('Distortion distance, squared',
                   'Maximal distortion distance, squared',
                   'Distortion index, squared')
  printCoefmat(di)

  # statistics
  cat('\nMultiple R-squared:', x$r.squared, '\tAdjusted R-squared:', x$adj.r.squared)
  cat('\nF-statistic:', x$F, 'on', x$df_model, 'and', x$df_residual, 'DF, p-value:', format.pval(x$p.value))
  cat('\nDifference in AIC to the null model:', x$dAIC)
  if (x$dAIC<2){
    cat('*')
  }
}


# Predicting  -------------------------------------------------------------

#' Predict method for Bidimensional Regression Model Fits
#'
#' Predicted values based on the bidimensional regressional model object.
#'
#' @param object an object of class "lm2"
#' @param newdata An optional two column data frame with independent variables.
#' If omitted, the fitted values are used.
#' @param ... optional arguments
#'
#' @return a two column data frame with predicted values for dependent variables.
#' @export
#'
#' @seealso \code{\link{lm2}}
#' @examples
#' lm2euc <- lm2(depV1+depV2~indepV1+indepV2, NakayaData, transformation = 'Euclidean')
#' predict(lm2euc, NakayaData[, 3:4])
predict.lm2 <-  function(object, newdata, ...) {
  # returning predictions for original independent variable values
  if (missing(newdata)){
    return(object$fitted_values)
  }

  # otherwise, checking dimensionality
  if (ncol(newdata)!=2) {
    stop('New data must be a two column matrix/data.frame.')
  }

  newdata$z <- 1
  newly_fitted <- data.matrix(newdata) %*% object$transformation_matrix
  newly_fitted <- newly_fitted[, 1:2]/newly_fitted[, 3]
  # colnames(newly_fitted) <- colnames(newdata)[1:2]
  return(newly_fitted)
}


# Comparing models --------------------------------------------------------


#' Anova for lm2 objects
#'
#' Anova for lm2 objects, returns a table with pairwise comparisons
#' between models or, if only one model was supplied, with the null model.
#'
#'
#' @param object an object of class "lm2"
#' @param ... further objects of class "lm2"
#'
#' @return an anova data frame
#' @export
#'
#' @seealso \code{\link{lm2}}
#' @examples
#' lm2euc <- lm2(depV1+depV2~indepV1+indepV2, NakayaData, transformation = 'Euclidean')
#' lm2aff <- lm2(depV1+depV2~indepV1+indepV2, NakayaData, transformation = 'Affine')
#' anova(lm2euc, lm2aff)
#' @importFrom stats pf
anova.lm2 <- function(object, ...)
{
  # checkings whether dots are lm2 objects
  dots_are_lm2 <- as.logical(vapply(list(...), is, NA, "lm2"))

  # merging models into a list
  all_models <- c(list(object), list(...)[dots_are_lm2])

  # checking that they all were fitted to the same data
  same_df <- vapply(all_models,  function(lm2object, df_name) { lm2object$data == df_name }, NA, all_models[[1]]$data)
  if (any(!same_df)){
    warning('Not all models are based on the same data as the first one, ignoring them')
  }
  all_models <- all_models[same_df]

  # checking for duplicate transforms
  transforms <- sapply(all_models, function(lm2object){lm2object$transformation})
  retain <- rep(TRUE, length(all_models))
  for(iModel in 2:length(transforms)){
    if (transforms[iModel] %in% transforms[1:(iModel-1)]){
      retain[iModel] <- FALSE
    }
  }
  if (any(!retain)){
    warning('DUplicate models, ignoring them')
  }
  all_models <- all_models[retain]

  if (length(all_models)==1) {
    # it all boiled down to a single model, for which we actually already computed statistics relative to the null model
    # thus, we just copy the numbers into a table
    anova_tbl <- data.frame(dAIC = object$dAIC,
                            df1 = as.integer(object$df_model),
                            df2 = as.integer(object$df_residual),
                            F= object$F,
                            p.value= object$p.value)
    row.names(anova_tbl) <- c(paste(object$transformation, 'null', sep= ' vs. '))
  }
  else {
    # we have more than one! First, let's order them based on complexity
    predictorsN <- sapply(all_models, function(lm2object){lm2object$npredictors})
    all_models <- all_models[order(predictorsN)]


    # Let's do pairwise comparisons.
    comparisonsN <- (length(all_models) * (length(all_models)-1))/2
    anova_tbl <- data.frame(dAIC = rep(NA, comparisonsN),
                            df1 = rep(NA, comparisonsN),
                            df2 = rep(NA, comparisonsN),
                            F= rep(NA, comparisonsN),
                            p.value= rep(NA, comparisonsN))
    iRow <- 1
    pairs_labels <- rep('', comparisonsN)
    for(iModel1 in 1:(length(all_models)-1)){
      for(iModel2 in (iModel1+1):length(all_models)){
        pairs_labels[iRow] <- paste(all_models[[iModel1]]$transformation, all_models[[iModel2]]$transformation, sep=' vs. ')

        anova_tbl$df1[iRow]<- all_models[[iModel2]]$df_model-all_models[[iModel1]]$df_model
        anova_tbl$df2[iRow]<- all_models[[iModel2]]$df_residual
        anova_tbl$F[iRow] <- (anova_tbl$df2[iRow]/anova_tbl$df1[iRow])*((all_models[[iModel2]]$r.squared-all_models[[iModel1]]$r.squared)/(1-all_models[[iModel2]]$r.squared))
        anova_tbl$p.value[iRow]<- pf(anova_tbl$F[iRow], anova_tbl$df1[iRow], anova_tbl$df2[iRow], lower.tail = FALSE, log.p = FALSE)
        anova_tbl$dAIC[iRow] <- 2*nrow(all_models[[iModel2]]$fitted_values)*log((1-all_models[[iModel2]]$r.squared)/
                                (1-all_models[[iModel1]]$r.squared))+2*(all_models[[iModel2]]$npredictors-all_models[[iModel1]]$npredictors)


        iRow <- iRow + 1
      }
      row.names(anova_tbl) <- pairs_labels
    }
  }

  # packaging
  anova_object <- list(anova_table = anova_tbl)
  class(anova_object) <- 'anova.lm2'
  return(anova_object)
}


#' @export
print.anova.lm2 <- function(x, ...){
  cat('Bidimensional regression:\n')
  printCoefmat(x$anova_table, cs.ind = c(1,4), P.values= TRUE, has.Pvalue=TRUE, na.print = '')
}
