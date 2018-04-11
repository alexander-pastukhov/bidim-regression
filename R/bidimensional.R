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
#' \code{A} and \code{B} are dependent and \code{C} and \code{D} are indepdent variables
#' @param data a data frame containing variables for the model.
#' @param transformation the transformation to be used, either \code{'euclidean'} or \code{'affine'}.
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
#' \item{\code{lm}}{an underlying \link[=lm]{linear model}}
#' \item{\code{formula}}{formula, describing input and output columns}
#' \item{\code{data}}{data used to fit the model}
#' \item{\code{Call}}{function call information, incorporates the \code{formula}, \code{transformation}, and \code{data}.}

#' @export
#'
#' @examples
#' nakayaAffine <- lm2(depV1 + depV2 ~ indepV1 + indepV2, NakayaData, 'Affine')
lm2 <- function(formula, data, transformation) { UseMethod("lm2") }

#' @export
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
  if (!is.character(transformation) || !(tolower(transformation) %in% c('euclidean', 'affine'))){
    stop("unknown transformation, please use either 'euclidean' or 'affine'")
  }


  # Extract variables from dataframe ----------------------------------------
  model_formula <- Formula::Formula(formula)
  DV <- Formula::model.part(model_formula, data = data, lhs = 1)
  IV <- Formula::model.part(model_formula, data = data, rhs = 1)


  # Fit the model -----------------------------------------------------------
  lm2model <- lm2fit(cbind(DV, IV), tolower(transformation))



  # Common information ------------------------------------------------------

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

  # adding information about the call
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
#' two are indepdent variables
#' @param transformation the transformation to be used, either \code{'euclidean'} or \code{'affine'}.
#'
#' @return returns an object of class "lm2", see \code{\link{lm2}}
#' for the description.
#'
#' @keywords internal
lm2fit <- function(data, transformation){
  lm2model <- switch(transformation,
                     euclidean = lm2euclidean(data),
                     affine = lm2affine(data),
                     stop("unknown transformation, please use either 'euclidean' or 'affine'"))
  class(lm2model) <- 'lm2'




  return(lm2model)
}


# Euclidean ---------------------------------------------------------------


#' Computes model for the euclidean transformation
#'
#' @param data the preprocessed data frame from \code{\link{lm2}} function,
#' so that the first two columns are the dependent variables and the other
#' two are indepdent variables
#'
#' @return object with transformation specific data to be supplemented with further stats
#' @keywords internal
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
#' two are indepdent variables
#'
#' @return object with transformation specific data to be supplemented with further stats
#' @keywords internal
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


# Printing and summary ----------------------------------------------------

#' @export
print.lm2 <- function(object){
  cat(sprintf('Call:\n'))
  cat(deparse(object$Call))
  cat('\n\n')
  cat('Coefficients:\n')
  coeff <- data.frame(as.list(object$coeff))
  rownames(coeff) <- ''
  printCoefmat(coeff)

  # transformed coefficients, if applicable
  if ("transformed_coeff" %in% names(object)){
    cat('\nTransformed coefficients:\n')
    transformed_coeff <- data.frame(as.list(object$transformed_coeff))
    rownames(transformed_coeff) <- ''
    printCoefmat(transformed_coeff)
  }

  # correlation strength
  cat('\nMultiple R-squared:', object$r.squared, '\tAdjusted R-squared:', object$adj.r.squared)
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
summary.lm2 <- function(object){
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
print.summary.lm2 <- function(object){
  cat(sprintf('Call:\n'))
  cat(deparse(object$Call))
  cat('\n\n')
  cat('Coefficients:\n')
  printCoefmat(object$coeff)

  # transformed coefficients
  if ('transformed_coeff' %in% names(object)){
    cat('\nTransformed coefficients:\n')
    transformed_coeff <- data.frame(as.list(object$transformed_coeff))
    rownames(transformed_coeff) <- ''
    printCoefmat(transformed_coeff)
  }

  # distortion index
  cat('\nDistortion index:\n')
  di <- t(object$distortion_index)
  rownames(di)<- c('Distortion distance, squared',
                   'Maximal distortion distance, squared',
                   'Distortion index, squared')
  printCoefmat(di)

  # statistics
  cat('\nMultiple R-squared:', object$r.squared, '\tAdjusted R-squared:', object$adj.r.squared)
  cat('\nF-statistic:', object$F, 'on', object$df_model, 'and', object$df_residual, 'DF, p-value:', format.pval(object$p.value))
  cat('\nDifference in AIC to the null model:', object$dAIC)
  if (object$dAIC<2){
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
#'
#' @return a two column data frame with predicted values for dependent variables.
#' @export
#'
#' @seealso \code{\link{lm2}}
#' @examples
#' lm2euc <- lm2(depV1+depV2~indepV1+indepV2, NakayaData, transformation = 'Euclidean')
#' predict(lm2euc, NakayaData[, 3:4])
predict.lm2 <-  function(object, newdata) {
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
  colnames(newly_fitted) <- colnames(newdata)[1:2]
  return(newly_fitted)
}


# Comparing models --------------------------------------------------------


#' Anova for lm2 objects
#'
#' Anova for lm2 objects.
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
anova.lm2 <- function(object, ...)
{
  # checkings whether dots are lm2 objects
  dots_are_lm2 <- as.logical(vapply(list(...), is, NA, "lm2"))

  if (!any(dots_are_lm2)) {
    # single model, for which we actually already computed statistics relative to the null model
    # thus, we just copy the numbers into a table
    anova_tbl <- data.frame(dAIC = object$dAIC,
                            df1 = as.integer(object$df_model),
                            df2 = as.integer(object$df_residual),
                            F= object$F,
                            p.value= object$p.value)
    row.names(anova_tbl) <- c(paste(object$transformation, 'null', sep= ' x '))

    anova_object <- list(anova_table = anova_tbl)
    class(anova_object) <- 'anova.lm2'
    return(anova_object)
  }
}

#   # sorting model based on number of predictors
#   npredictors <- c()
#   for(current.model in models){
#     npredictors <- c(npredictors, current.model$npredictors)
#   }
#   models <- models[order(npredictors)]
#
#   # comparing each model to the previous one
#   anova.table <- data.frame(transformation= rep(NA, length(models)),
#                             df1= NA,
#                             df2= NA,
#                             dAIC= NA,
#                             F= NA,
#                             p.value= NA)
#   anova.table$transformation[1] <- tolower(models[[1]]$transformation)
#
#   for(i.Model in 2:length(models)){
#     # no point comparing same transformation
#     if (models[[i.Model-1]]$transformation==models[[i.Model]]$transformation){
#       next;
#     };
#
#     anova.table$transformation[i.Model] <- tolower(models[[i.Model]]$transformation)
#     anova.table$df1[i.Model]<- models[[i.Model]]$df1-models[[i.Model-1]]$df1
#     anova.table$df2[i.Model]<- models[[i.Model]]$df2
#     anova.table$F[i.Model] <- (anova.table$df2[i.Model]/anova.table$df1[i.Model])*((models[[i.Model]]$r.squared-models[[i.Model-1]]$r.squared)/(1-models[[i.Model]]$r.squared))
#     anova.table$p.value[i.Model]<- pf(anova.table$F[i.Model], anova.table$df1[i.Model], anova.table$df2[i.Model], lower.tail = FALSE, log.p = FALSE)
#     anova.table$dAIC[i.Model] <- 2*nrow(models[[i.Model]]$fitted.values)*log((1-models[[i.Model]]$r.squared)/
#                                                                                (1-models[[i.Model-1]]$r.squared))+2*(models[[i.Model]]$npredictors-models[[i.Model-1]]$npredictors)
#   }
#
#   anova.table <- subset(anova.table, !is.na(transformation))
#   row.names(anova.table)<-anova.table$transformation
#   anova.table <- subset(anova.table, select = -transformation)
#
#   object <- list(anova.table= anova.table, dimsN= models[[1]]$dimsN)
#
#   class(object) <- 'anova.rmNDim'
#   return(object)
# }
#
# rm2printAnova <- function(object){
#   cat('Bidimensional regression:\n')
#   printCoefmat(object$anova.table, P.values= TRUE, has.Pvalue=TRUE, na.print = '')
# }
#


#' @export
print.anova.lm2 <- function(object){
  cat('Bidimensional regression:\n')
  printCoefmat(object$anova_table, cs.ind = c(1,4), P.values= TRUE, has.Pvalue=TRUE, na.print = '')
}
