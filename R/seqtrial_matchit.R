#' Formatting a dataset for the sequential trials approach
#'
#' @param id_var The id variable given as a character vector.
#' @param time_var The time variable given as a character vector. Must be a numeric variable.
#' @param treated_var The variable indicating the time of treatment initiation. Must be coercible to 0 / 1. Takes the value 0 until the time of treatment initiation, and 1 on and after the time of treatment initiation.
#' @param treatment_var The treatment variable given as a character vector. Must be NA before treatment initiated. Must be NULL if comparison is treated / untreated.
#' @param treatment_vals Values to define the treatment groups in a "treatment A vs treatment B" comparison. Must match the nonmissing values in treated_var. Must be NULL if treatment_var is NULL.
#' @param matching_vars The matching variables given as a character vector. Must be factor variables (this will be updated when matching strategies extended).
#' @param other_vars Other variables that may be used as covariates in the model or to assess the balance after matching. Can be NULL.
#' @param censor_vars Censoring variables given as a character vector. Variables must be coercible to 0 / 1, taking the value 0 before the censoring event and 1 at the time of the event. Only the first occurrence of the event will be used, so it does not matter whether the value returns to 0 after the event, or remains at 1.
#' @param outcome_vars Outcome variables given as a character vector. Must specify at least one outcome variable. Variables must be coercible to 0 / 1, taking the value 0 before the outcome event and 1 at the time of the event. Only the first occurrence of the event will be used, so it does not matter whether the value returns to 0 after the event, or remains at 1.
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
#' matching_vars = c("age_grp", "biomarker", "sex"),
#' other_vars = NULL,
#' censor_vars = NULL,
#' outcome_vars = death,
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
#' matching_vars = c("age_grp", "biomarker", "sex"),
#' other_vars = NULL,
#' data = dummydata
#' )
seqtrial_formatter <- function(id_var,
                               time_var,
                               treated_var,
                               treatment_var = NULL,
                               treatment_vals = NULL,
                               matching_vars,
                               other_vars = NULL,
                               censor_vars = NULL, # TODO add a censor var to dummydata so this can be tested
                               outcome_vars,
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

  all_vars <- c(id_var, time_var, treated_var, treatment_var, matching_vars, other_vars)

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

  # check matching_vars are factors
  # (for now as only exact matching, but will extend to other matching techniques)
  check_factors <- data %>%
    dplyr::select(dplyr::all_of(matching_vars)) %>%
    purrr::map_lgl(\(x) is.factor(x))
  if (!all(check_factors)) {
    stop(
      paste0(
        paste(names(check_factors)[!check_factors], collapse = " and "),
        " are not factors."
      )
    )
  }


  # check for missing values
  check_missing <- data %>%
    dplyr::select(dplyr::all_of(c(id_var, time_var, treated_var, matching_vars))) %>%
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
  attributes(data)$matching_vars <- matching_vars
  attributes(data)$other_vars <- other_vars

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

  data_reatment_initiation_times <- data %>%
    dplyr::filter(treated == 1) %>%
    dplyr::arrange(id, time) %>%
    dplyr::group_by(id) %>%
    # slice is quicker than taking the min
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(all_of(vars)) %>%
    dplyr::arrange(time)

  return(data_reatment_initiation_times)

}

#' Perform matching according to the sequential trial approach for target trial emulation
#'
#' @param data A data frame with class "seqtrial_fmt".
#'
#' @return A data frame with class "seqtrial_matched_fmt"
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
#' matching_vars = c("age_grp", "biomarker", "sex"),
#' other_vars = NULL,
#' data = dummydata
#' )
#' data_seqtrial_tu_matched <- seqtrial_matchit(data_seqtrial_tu)
#'
#' # for treatment A vs treatment B comparison
#' data_seqtrial_ab <- seqtrial_formatter(
#' id_var ="id",
#' time_var = "time",
#' treated_var = "treated",
#' treatment_var = "treatment",
#' treatment_vals = c("A", "B"),
#' matching_vars = c("age_grp", "biomarker", "sex"),
#' other_vars = NULL,
#' data = dummydata
#' )
#' data_seqtrial_ab_matched <- seqtrial_matchit(data_seqtrial_ab)
seqtrial_matchit <- function(data) {

  # TODO allow grouping dates for trials (e.g. data are daily, but trials are weekly)

  seqtrial_fmt_checks(data)

  # get all the times at which people initiate treatment
  data_treatment_initiation_times <- treatment_initiation_times(data)

  # create empty objects for output
  data_successful_matches <- list()
  previously_matched_ids <- vector(mode = mode(data$id), length=0)

  # create function to catch errors
  safely_matchit <- purrr::safely(MatchIt::matchit)

  message("Sequential trials matching report:")

  for (trial_start in unique(data_treatment_initiation_times$time)) {

    message("---- time = ", trial_start, " ----")

    treated_ids_i <- data_treatment_initiation_times %>%
      dplyr::filter(time == trial_start) %>%
      dplyr::pull(id)

    data_i <- data %>%
      dplyr::filter(time == trial_start)

    if (attributes(data)$comparison == "treated vs untreated") {

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

      group_var <- "treated"
      group_var_labels <- c("untreated" = 0, "treated" = 1)

    }

    if (attributes(data)$comparison == "treatment A vs treatment B") {

      match_candidates_i <- data_i %>%
        dplyr::mutate(across(
          treatment,
          ~as.integer(.x == attributes(data)$treatment_vals[2])
          ))

      group_var <- "treatment"
      group_var_labels <- c(0, 1) # in formatter, define treatment A and treatment B
      names(group_var_labels) <- attributes(data)$treatment_vals

    }

    # TODO in future allow user to specify matching options in seqtrial()
    # TODO sort out message printing - no newline after matchit messages
    # run match algorithm
    obj_matchit_i <-
      safely_matchit(
        formula = stats::as.formula(paste0(group_var, " ~ 1")),
        data = match_candidates_i,
        replace = FALSE,
        estimand = "ATT",
        exact = attributes(data)$matching_vars,
        m.order = "random",
        # verbose = TRUE,
        ratio = 1L # 1:1 match
      )[[1]]

    if(is.null(obj_matchit_i)) {
      message("Skipping trial - no matches found.")
      next
    }

    data_successful_matches[[as.character(trial_start)]] <- dplyr::tibble(
      id = match_candidates_i$id,
      matched = !is.na(obj_matchit_i$subclass),
      match_id = as.integer(as.character(obj_matchit_i$subclass)),
      group = obj_matchit_i$treat,
      weight = obj_matchit_i$weights,
      trial_start = trial_start,
    ) %>%
      dplyr::filter(matched) %>%
      dplyr::mutate(dplyr::across(
        group,
        ~ factor(.x, levels = unname(group_var_labels), labels = names(group_var_labels))
        )) %>%
      dplyr::arrange(group) %>%
      tidyr::pivot_wider(
        names_from = group,
        values_from = id,
      )

    if (attributes(data)$comparison == "treated vs untreated") {

      data_successful_matches[[as.character(trial_start)]] <-
        data_successful_matches[[as.character(trial_start)]] %>%
        dplyr::left_join(
          data_treatment_initiation_times %>%
            dplyr::select(id, untreated_is_treated_time = time),
          by = c("untreated" = "id")
        )

      # ids for individuals matched as a untreated in all trials so far
      previously_matched_ids <- c(
        previously_matched_ids,
        data_successful_matches[[as.character(trial_start)]]$untreated
      )

    }

  }

  # currently a list with a dataset for each trial, bind into one dataset
  data_successful_matches <- dplyr::bind_rows(data_successful_matches)

  # check for duplicates in groups
  for (i in names(group_var_labels)) {
    if (any(duplicated(data_successful_matches[[i]]))) {
      stop("Duplicate values in group: ", i)
    }
  }

  data_matched <- data_successful_matches %>%
    tidyr::pivot_longer(
      cols = names(group_var_labels),
      names_to = "group",
      values_to = "id"
    ) %>%
    # because match_id only unique within trial_start
    dplyr::group_by(trial_start, match_id) %>%
    dplyr::mutate(match_id = dplyr::cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(across(
      group,
      ~factor(.x, levels = names(group_var_labels))
      )) %>%
    dplyr::select(
      id, match_id, group, trial_start,
      any_of("untreated_is_treated_time")
      ) %>%
    tibble::new_tibble(class = "seqtrial_matched_fmt")

  return(data_matched)

}


# TODO create function to summarise matching success
# TODO create function to assess balance across matched groups
# TODO create function to add time-to-event data to a seqtrial_matched_fmt object
