#' Formatting a dataset for the sequential trials approach
#'
#' @param id_var The id variable given as a character vector.
#' @param time_var The time variable given as a character vector. Must be a numeric variable.
#' @param treated_var The variable indicating the time of treatment initiation. Must be coercible to 0 / 1. Takes the value 0 until the time of treatment initiation, and 1 on and after the time of treatment initiation.
#' @param treatment_var The treatment variable given as a character vector. Must be NA before treatment initiated. Must be NULL if comparison is treated / untreated.
#' @param treatment_vals Values to define the treatment groups in a "treatment A vs treatment B" comparison. Must match the nonmissing values in treated_var. Must be NULL if treatment_var is NULL.
#' @param outcome_vars Outcome variables given as a character vector. Must specify at least one outcome variable. Variables must be coercible to 0 / 1, taking the value 0 before the outcome event and 1 at the time of the event. Only the first occurrence of the event will be used, so it does not matter whether the value returns to 0 after the event, or remains at 1.
#' @param censor_vars Censoring variables given as a character vector. Variables must be coercible to 0 / 1, taking the value 0 before the censoring event and 1 at the time of the event. Only the first occurrence of the event will be used, so it does not matter whether the value returns to 0 after the event, or remains at 1.
#' @param keep_vars The variables to keep in the dataset given as a character vector. Any variables used in matching, assessing balance post-matching or used as covariates in the model should be specified here. At least one variable must be specified.
#' @param data A data frame in which these variables exist. All variables must be in this data frame.
#'
#' @return A data frame with class "seqtrial_fmt"
#' @export
#'
#' @examples
#' data("dummydata")
#'
#' # for treated vs untreated comparison
#' data_seqtrial_tu <- seqtrial_formatter(
#' id_var ="id",
#' time_var = "time",
#' treated_var = "treated",
#' outcome_vars = death,
#' censor_vars = NULL,
#' keep_vars = c("age_grp", "biomarker_grp", "sex"),
#' data = dummydata
#' )
#'
#' # for treatment A vs treatment B comparison
#' data_seqtrial_ab <- seqtrial_formatter(
#' id_var ="id",
#' time_var = "time",
#' treated_var = "treated",
#' treatment_var = "treatment",
#' treatment_vals = c("A", "B"),
#' keep_vars = c("age_grp", "biomarker_grp", "sex"),
#' data = dummydata
#' )
seqtrial_formatter <- function(id_var,
                               time_var,
                               treated_var,
                               treatment_var = NULL,
                               treatment_vals = NULL,
                               outcome_vars = NULL,
                               censor_vars = NULL, # TODO add a censor var to dummydata so this can be tested
                               keep_vars,
                               data) {

  # TODO handle censoring and outcome variables

  # create symbols for programming
  id_sym <- rlang::sym(id_var)
  time_sym <- rlang::sym(time_var)
  treated_sym <- rlang::sym(treated_var)
  treatment_sym <- NULL
  if (!is.null(treatment_var)) treatment_sym <- rlang::sym(treatment_var)

  stopifnot(
    "`data` must have class \"data.frame\"" =
      "data.frame" %in% class(data)
    )

  all_vars <- c(id_var, time_var, treated_var, treatment_var, outcome_vars, censor_vars, keep_vars)

  stopifnot(
    "All supplied vars must be in `data`" =
      all(all_vars %in% names(data))
  )

  # TODO (extend to handle dates for time_var)
  stopifnot(
    "`time_var` must be numeric" =
      is.numeric(data[[time_var]])
    )

  stopifnot(
    "`treated_var` must have two unique values that can be coerced to ingeter values 0 and 1" =
      all(as.integer(unique(data[[treated_var]])) == c(0,1))
    )

  # treated_var can be 0 or 1 throughout, or switch from 0 to 1, but not from 1 to 0
  treatment_stopped <- data %>%
    dplyr::group_by(!!id_sym) %>%
    dplyr::mutate(treated_diff = c(NA_integer_, diff(as.integer(treated)))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(treated_diff < 0) %>%
    nrow()
  stopifnot(
    "Treatment may be initiated once, then must be continued for the remainder of followup" =
      treatment_stopped == 0
  )


  stopifnot(
    "If one of `treatment_var` and `treatment_vals` specified, both must be specified." =
      is.null(treatment_var) == is.null(treatment_vals)
  )

  if (!is.null(treatment_var)) {

    treatment_var_levs <- unique((data[[treatment_var]]))
    stopifnot(
      "`treatment_var` must have two unique non-missing values." =
        length(treatment_var_levs[!is.na(treatment_var_levs)]) >= 2
    )
    stopifnot(
      "`treatment_vals` must have length 2." =
        length(treatment_vals) == 2
    )
    stopifnot(
      "`treatment_vals` must be coerced to values in `treatment_var`." =
        all(treatment_vals %in% treatment_var_levs)
    )

    treatment_when_untreated <- data %>%
      dplyr::filter(
        (as.integer(!!treated_sym) == 0L) & !is.na(!!treatment_sym)
        ) %>%
      nrow()
      stopifnot(
        "`treatment_var` must be missing when `treated_var` is 0." =
          treatment_when_untreated == 0
      )

      multiple_treatments <- data %>%
        dplyr::filter(as.integer(!!treated_sym) == 0L) %>%
        dplyr::group_by(!!id_sym) %>%
        dplyr::distinct(!!treatment_sym) %>%
        dplyr::count() %>%
        dplyr::ungroup() %>%
        dplyr::filter(n>1)
      stopifnot(
        "Once treatment is intiated, `treatment_var` can only take one value per individual." =
          multiple_treatments == 0
      )

  }

  # check for missing values
  check_missing <- data %>%
    dplyr::select(dplyr::all_of(c(id_var, time_var, treated_var, outcome_vars, censor_vars, keep_vars))) %>%
    purrr::map_lgl(\(x) any(is.na(x)))
  if (any(check_missing)) {
    stop(
      paste0(
        "Missing values in ",
        paste(names(check_missing)[check_missing], collapse = " and "),
        "."
      )
    )
  }

  names_map <- list(
    "id" = id_var,
    "time" = time_var,
    "treated" = treated_var,
    "treatment" = treatment_var
  )

  # only keep the specified variables and rename
  data <- data %>%
    dplyr::select(dplyr::all_of(all_vars)) %>%
    dplyr::rename(
      "id" = !!id_sym,
      "time" = !!time_sym,
      "treated" = !!treated_sym,
      "treatment" = !!treatment_sym
    ) %>%
    dplyr::mutate(dplyr::across(treated, as.integer))

  attributes(data)$names_map <- unlist(names_map) # for reverting to original names in seqtrial_unformatter
  attributes(data)$comparison <- ifelse(is.null(treatment_var), "treated vs untreated", "treatment A vs treatment B")
  attributes(data)$treatment_vals <- treatment_vals

  if (attributes(data)$comparison == "treatment A vs treatment B") {
    data <- data %>%
      dplyr::filter(treated == 1)
  }

  data <- data %>%
    # use tibble::new_tibble() otherwise dplyr functions applied to data fail
    tibble::new_tibble(class = "seqtrial_fmt")

  return(data)

}

seqtrial_fmt_checks <- function(data) {

  stopifnot(
    "`data` must have class \"seqtrial_fmt\"" =
      "seqtrial_fmt" %in% class(data)
  )

  vars <- c("id", "time", "treated")
  if ("treatment_vals" %in% names(attributes(data))) {
    vars <- c(vars, "treatment")
  }

  error_message <- paste0("`data` must contain the following columns: ", paste0(vars, collapse = ", "))
  stopifnot(
    error_message =
      all(vars %in% names(data))
  )

}

# get all the times at which people initiate treatment
treatment_initiation_times <- function(data) {

  seqtrial_fmt_checks(data)

  vars <- c("id", "time", "treated")
  if ("treatment_vals" %in% names(attributes(data))) {
    vars <- c(vars, "treatment")
  }

  data_treatment_initiation_times <- data %>%
    dplyr::filter(treated == 1) %>%
    dplyr::arrange(id, time) %>%
    dplyr::group_by(id) %>%
    # slice is quicker than taking the min
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(all_of(vars)) %>%
    dplyr::arrange(time)

  return(data_treatment_initiation_times)

}

#' Perform matching according to the sequential trial approach for target trial emulation
#'
#' @param seqtrial_data A data frame with class "seqtrial_fmt".
#' @param formula a two-sided [`formula`] object containing the treatment and
#' covariates to be used in creating the distance measure used in the matching.
#' This formula will be supplied to the functions that estimate the distance
#' measure. The formula should be specified as `A ~ X1 + X2 + ...` where
#' `A` represents the treatment variable and `X1` and `X2` are
#' covariates.
#' @param \dots additional arguments passed to [`matchit`]
#'
#' @return A data frame with class "seqtrial_matched"
#' @export
#'
#' @examples
#' data("dummydata")
#'
#' # for treated vs untreated comparison
#' data_seqtrial_tu <- seqtrial_formatter(
#' id_var ="id",
#' time_var = "time",
#' treated_var = "treated",
#' outcome_vars = "death",
#' keep_vars = c("age_grp", "biomarker_grp", "sex"),
#' data = dummydata
#' )
#' # exact 1:1 matching on age_grp, biomarker_grp and sex
#' data_seqtrial_tu_matched <- seqtrial_matchit(
#' data_seqtrial_tu,
#' formula = treated ~ 1,
#' exact = c("age_grp", "biomarker_grp", "sex"),
#' replace=FALSE,
#' ratio=1,
#' m.order = "random",
#' estimand = "ATT",
#' )
#'
#' # for treatment A vs treatment B comparison
#' data_seqtrial_ab <- seqtrial_formatter(
#' id_var ="id",
#' time_var = "time",
#' treated_var = "treated",
#' treatment_var = "treatment",
#' treatment_vals = c("A", "B"),
#' outcome_vars = "death",
#' keep_vars = c("age_grp", "biomarker_grp", "sex"),
#' data = dummydata
#' )
#' # exact 1:1 matching on age_grp, biomarker_grp and sex
#' data_seqtrial_ab_matched <- seqtrial_matchit(
#' data_seqtrial_ab,
#' formula = treatment ~ 1,
#' exact = c("age_grp", "biomarker_grp", "sex"),
#' replace=FALSE,
#' ratio=1,
#' m.order = "random",
#' estimand = "ATT",
#' )
seqtrial_matchit <- function(seqtrial_data, formula, ...) {

  # TODO allow grouping dates for trials (e.g. data are daily, but trials are weekly)

  seqtrial_fmt_checks(seqtrial_data)

  # get LHS of formula as character
  tf <- terms(formula)
  tf.vars <- attributes(tf)$variables
  lhs.char <- as.character(tf.vars[[2]])
  # check that LHS is compatible with the comparison
  if (attributes(seqtrial_data)$comparison == "treated vs untreated") {
    stopifnot(
      "LHS of formula must be \"treated\" when comparison is \"treated vs untreated\""  =
        lhs.char == "treated"
    )
  }
  if (attributes(seqtrial_data)$comparison == "treatment A vs treatment B") {
  stopifnot(
    "LHS of formula must be \"treatment\" when comparison is \"treatment A vs treatment B\""  =
      lhs.char == "treatment"
  )
  }

  # get all the times at which people initiate treatment
  data_treatment_initiation_times <- treatment_initiation_times(seqtrial_data)

  # create empty objects for output
  # data_successful_matches <- list()
  obj_matchit <- list()
  previously_matched_ids <- vector(mode = mode(seqtrial_data$id), length=0)

  # create function to catch errors
  safely_matchit <- purrr::safely(MatchIt::matchit, quiet = TRUE)

  for (trial_start in unique(data_treatment_initiation_times$time)) {

    treated_ids_i <- data_treatment_initiation_times %>%
      dplyr::filter(time == trial_start) %>%
      dplyr::pull(id)

    data_i <- seqtrial_data %>%
      dplyr::filter(time == trial_start)

    if (attributes(seqtrial_data)$comparison == "treated vs untreated") {

      # data for match candidates
      match_candidates_i <- dplyr::bind_rows(
        # treated: people who initiated treatment on trial_start
        # (people can be treated if previously matched as a control)
        data_i %>%
          dplyr::filter(id %in% treated_ids_i),
        # untreated: people who remained untreated on trial_start and have not been
        # matched as untreated in a previous trial
        data_i %>%
          dplyr::filter(!(treated | id %in% previously_matched_ids))
      )

    }

    if (attributes(seqtrial_data)$comparison == "treatment A vs treatment B") {

      match_candidates_i <- data_i %>%
        dplyr::mutate(across(
          treatment,
          ~as.integer(.x == attributes(seqtrial_data)$treatment_vals[2])
          ))

    }

    # run match algorithm
    obj_matchit_i <- safely_matchit(
      data = match_candidates_i,
      formula = formula,
      ...
    )

    if(!is.null(obj_matchit_i$error)) {
      # print the errors
      message("Error for trial_start = ", trial_start, ": ", obj_matchit_i$error$message)
      # skip to the next trial
      next
    } else {
      # output the matchit output so user has option to use matchit s3 methods
      obj_matchit[[as.character(trial_start)]] <- obj_matchit_i$result
    }

  }

  class(obj_matchit) <- "seqtrial_matchit"

  return(obj_matchit)

}

#' View a balance summary of a `seqtrial_matchit` object for each trial
#'
#' @param seqtrial_matchit_obj a `seqtrial_matchit` object; the output of a call to [seqtrial_matchit()]
#' @param \dots additional arguments passed to [`summary.matchit`]
#'
#' @return For `seqtrial_matchit` objects, a list with the following components:
#'
#' \item{call}{the original call to [matchit()]}
#' \item{nn}{a tibble of the sample sizes in the original (unmatched) and
#' matched samples per trial}
#' \item{sum.all}{if `un = TRUE`, a tibble of balance statistics for each
#' covariate in the original (unmatched) sample per trial}
#' \item{sum.matched}{a tibble of balance statistics for each covariate in the
#' matched sample per trial}
#' \item{reduction}{if `improvement = TRUE`, a tibble of the percent
#' reduction in imbalance for each covariate in the matched sample}
#'
#' @examples
#' data("dummydata")
#'
#' # for treated vs untreated comparison
#' data_seqtrial_tu <- seqtrial_formatter(
#' id_var ="id",
#' time_var = "time",
#' treated_var = "treated",
#' outcome_vars = "death",
#' keep_vars = c("age_grp", "biomarker_grp", "sex"),
#' data = dummydata
#' )
#' # exact 1:1 matching on age_grp, biomarker_grp and sex
#' data_seqtrial_tu_matched <- seqtrial_matchit(
#' data_seqtrial_tu,
#' formula = treated ~ 1,
#' exact = c("age_grp", "biomarker_grp", "sex"),
#' replace=FALSE,
#' ratio=1,
#' m.order = "random",
#' estimand = "ATT",
#' )
#' summary(data_seqtrial_tu_matched)
#'
#' @exportS3Method summary seqtrial_matchit
summary.seqtrial_matchit <- function(seqtrial_matchit_obj, ...) {

  # TODO maybe extend to handle match.subclass objects?

  res <- lapply(
    seqtrial_matchit_obj,
    function(x) summary(x, ...)
  )
  names(res) <- names(seqtrial_matchit_obj)

  print(res)

  # combine across trials
  res_combine_names <- c("nn", "sum.all", "sum.matched", "reduction")
  combine_summaries <- function(i) {
    lapply(
      names(res),
      function(x) tibble::as_tibble(res[[x]][[i]], rownames = "descr") %>%
        dplyr::mutate(trial_start = as.integer(x), .before=1)
    ) %>%
      dplyr::bind_rows()
  }
  res_combine <- c(
    "call" = res[[1]]$call,
    lapply(res_combine_names, combine_summaries)
  )
  names(res_combine) <- c("call", res_combine_names)

  class(res_combine) <- "summary.seqtrial_matchit"

  return(res_combine)

}

# TODO create functions to summarise/plot matching success
# TODO create function to plot balance across matched groups
# TODO create function to add time-to-event data to a seqtrial_matched_fmt object
