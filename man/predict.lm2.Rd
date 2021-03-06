% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bidimensional.R
\name{predict.lm2}
\alias{predict.lm2}
\title{Predict method for Bidimensional Regression Model Fits}
\usage{
\method{predict}{lm2}(object, newdata, ...)
}
\arguments{
\item{object}{an object of class "lm2"}

\item{newdata}{An optional two column data frame with independent variables.
If omitted, the fitted values are used.}

\item{...}{optional arguments}
}
\value{
a two column data frame with predicted values for dependent variables.
}
\description{
Predicted values based on the bidimensional regressional model object.
}
\examples{
lm2euc <- lm2(depV1+depV2~indepV1+indepV2, NakayaData, transformation = 'Euclidean')
predict(lm2euc, NakayaData[, 3:4])
}
\seealso{
\code{\link{lm2}}
}
