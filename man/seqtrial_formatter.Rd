% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seqtrial_matchit.R
\name{seqtrial_formatter}
\alias{seqtrial_formatter}
\title{Formatting a dataset for the sequential trials approach}
\usage{
seqtrial_formatter(
  id_var,
  time_var,
  treated_var,
  treatment_var = NULL,
  treatment_vals = NULL,
  matching_vars,
  other_vars = NULL,
  data
)
}
\arguments{
\item{id_var}{The id variable given as a character vector.}

\item{time_var}{The time variable given as a character vector. Must be a numeric variable.}

\item{treated_var}{The variable indicating the time of treatment initiation. Must be coercible to 0 / 1.}

\item{treatment_var}{The treatment variable given as a character vector. Must be NA before treatment initiated. Must be NULL if comparison is treated / untreated.}

\item{matching_vars}{The matching variables given as a character vector. Must be a factor variable.}

\item{other_vars}{Other variables that may be used as covariates in the model or to assess the balance after matching. Can be NULL.}

\item{data}{A data frame in which these variables exist. All variables must be in this data frame.}

\item{treatment_val}{A value to define the treated group. Must be one of the unique values in treated_var. Must be NULL if comparison is treated / untreated.}
}
\value{
A data frame with class "seqtrial_fmt"
}
\description{
Formatting a dataset for the sequential trials approach
}
\examples{
data("dummydata")

# for treated vs untreated comparison
data_seqtrial_tu <- seqtrial_formatter(
id_var ="id",
time_var = "time",
treated_var = "treated",
matching_vars = c("age_grp", "biomarker", "sex"),
other_vars = NULL,
data = dummydata
)

# for treatment A vs treatment B comparison
data_seqtrial_ab <- seqtrial_formatter(
id_var ="id",
time_var = "time",
treated_var = "treated",
treatment_var = "treatment",
treatment_vals = c("A", "B"),
matching_vars = c("age_grp", "biomarker", "sex"),
other_vars = NULL,
data = dummydata
)
}
