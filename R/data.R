#' Eye gaze calibration data
#'
#' A dataset containing a monocular eye gaze recording with calibration sequence.
#' Courtesy of Bamberger Baby Institut: \href{https://www.uni-bamberg.de/entwicklungspsychologie/transfer/unsichtbar-fuer-bambi/babyforschung-bambi/}{BamBI}.
#'
#' @format A data frame with 365 rows and 6 variables:
#' \describe{
#'   \item{time}{sample timestamp, in milliseconds}
#'   \item{x, y}{recorded gaze, in internal eye tracker units}
#'   \item{target_x, target_y}{location of the calibration target on the screen, in pixels}
#'   \item{target}{index of the target within the sequence}
#'
#'   ...
#' }
"EyegazeData"
