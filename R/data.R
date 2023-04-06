#' Time-to-event dummy data with time-updating variables
#'
#' The dummy data is generated to mimic the following observational study:
#' Patients eligible for treatments A or B are recruited and followed up for a maximum of 24 months.
#' Age and sex are measured at baseline; levels of a biomarker are measured at each monthly visit.
#' Treatment is more likely to be initiated when biomarker level > 1200.
#' Higher biomarker level associated with higher risk of death (not currently the case in dummydata, need to update).
#'
#'
#' @format ## `dummydata`
#' A data frame with  10,753 rows and 8 columns:
#' \describe{
#'   \item{id}{Patient id}
#'   \item{time}{Time since study baseline in months}
#'   \item{sex}{Patient's sex (measured at baseline)}
#'   \item{age_grp}{Patient's age group (measured at baseline)}
#'   \item{biomarker}{Patient's biomarker level (measured at "time")}
#'   \item{treated}{Whether the patient is treated}
#'   \item{treatment}{Which treatment}
#'   \item{death}{Had the patient died by "time"}
#'   ...
#' }
"dummydata"
