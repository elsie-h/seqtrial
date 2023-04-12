# get all the times at which people initiate treatment
alltreated <- function(id_sym,
                       time_sym,
                       treated_sym,
                       data) {

  data_alltreated <- data %>%
    dplyr::filter(!!treated_sym) %>%
    dplyr::arrange(!!id_sym, !!time_sym) %>%
    dplyr::group_by(!!id_sym) %>%
    # slice is quicker than taking the min
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(!!time_sym)

  return(data_alltreated)

}

#' Formatting a dataset for the sequential trials approach
#'
#' @param id_var The id variable given as a character vector.
#' @param time_var The time variable given as a character vector. Must be a numeric variable.
#' @param treated_var TODO
#' @param treatment_var The treatment variable given as a character vector. Must be NA before treatment initiated.
#' @param treatment_val A value to define the treated group. Must be one of the unique values in treated_var.
#' @param matching_vars The matching variables given as a character vector. Must be a factor variable.
#' @param other_vars TODO
#' @param data A data frame in which these variables exist. All variables must be in this data frame.
#'
#' @return A data frame with class "seqtrial_fmt"
#' @export
#'
#' @examples
#' seqtrial_example <- seqtrial_formatter(
#' id_var ="id",
#' time_var = "time",
#' treated_var = "treated",
#' treatment_var = "treatment",
#' treatment_val = "A",
#' matching_vars = c("age_grp", "biomarker", "sex"),
#' other_vars = NULL,
#' data = dummydata
#' )
seqtrial_formatter <- function(id_var,
                               time_var,
                               treated_var, # 0 until treatment initiation, 1 at all timepoints after
                               treatment_var = NULL, # var specifying type of treatment. If NULL assume treated vs control
                               treatment_val = NULL, # the value of the treatment variable to be coded as 1 (all other values coded as 0)
                               matching_vars,
                               other_vars = NULL, # other variables that might be used when assessing balance across matched groups, or as covariates in models
                               data) {

  # TODO allow grouping dates for trials (e.g. data are daily, but trials are weekly)

  # create symbols for programming
  id_sym <- rlang::sym(id_var)
  time_sym <- rlang::sym(time_var)
  treated_sym <- rlang::sym(treated_var)
  treatment_sym <- rlang::sym(treatment_var)

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
    "If one of `treatment_var` and `treatment_val` specified, both must be specified." =
      is.null(treatment_var) == is.null(treatment_val)
  )

  if (!is.null(treated_var)) {

    treatment_var_levs <- unique((data[[treatment_var]]))
    stopifnot(
      "`treatment_var` must have at least two unique values that are non-missing." =
        length(treatment_var_levs[!is.na(treatment_var_levs)]) >= 2
    )
    stopifnot(
      "`treatment_val` must have length 1." =
        length(treatment_val) == 1
    )
    stopifnot(
      "`treatment_val` must be coerced to a value in `treatment_var`." =
        treatment_val %in% treatment_var_levs
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
      "id" = !!rlang::sym(names_map$id),
      "time" = !!rlang::sym(names_map$time),
      "treated" = !!rlang::sym(names_map$treated),
      "treatment" = !!rlang::sym(names_map$treatment)
    ) %>%
    dplyr::mutate(dplyr::across(treated, as.integer))

  attributes(data)$names_map <- unlist(names_map) # for reverting to original names in seqtrial_unformatter
  attributes(data)$treatment_val <- treatment_val

  class(data) <- append(class(data), "seqtrial_fmt")

  return(data)

}


#' Perform matching according to the sequential trial approach for target trial emulation
#'
#' @param id_var A variable given as a character vector.
#' @param time_var A variable given as a character vector.
#' @param treated_var A variable given as a character vector.
#' @param matching_vars A variable given as a character vector.
#' @param data A data frame in which these variables exist. All variables must be in this data frame.
#'
#' @return A data frame with class
#' @export
#'
#' @examples
#' seqtrial_example <- seqtrial_matchit(
#' id_var ="id",
#' time_var = "time",
#' treated_var = "treated",
#' matching_vars = c("age_grp", "biomarker", "sex"),
#' data = dummydata
#' )
seqtrial_matchit <- function(data) {





  # get all the times at which people initiate treatment
  data_alltreated <- alltreated(id_sym, time_sym, treated_sym, data)

  # create empty objects for output
  data_successfulmatches <- list()
  previouslymatched_ids <- vector(mode = mode(data[[id_var]]), length=0)

  message("Sequential trials matching report:")

  for (trialstart in unique(data_alltreated[[time_var]])) {

    message("---- ", time_var, " = ", trialstart, " ----")

    treated_ids_i <- data_alltreated %>%
      dplyr::filter(!!time_sym == trialstart) %>%
      dplyr::pull(!!id_sym)

    data_i <- data %>%
      dplyr::filter(!!time_sym == trialstart)

    # data for match candidates
    match_candidates_i <- dplyr::bind_rows(
      # treated: people who initiated treatment on trialstart
      # (people can be treated if previously matched as a control)
      data_i %>%
        dplyr::filter(!!id_sym %in% treated_ids_i),
      # control: people who remained untreated on trialstart and have not been
      # matched as a control in a previous trial
      data_i %>%
        dplyr::filter(!(!!treated_sym) & !(!!id_sym %in% previouslymatched_ids))
    )

    # create function to catch errors
    safely_matchit <- purrr::safely(MatchIt::matchit)
    # TODO in future allow user to specify matching options in seqtrial()
    # TODO sort out message printing - no newline after matchit messages
    # run match algorithm
    obj_matchit_i <-
      safely_matchit(
        formula = stats::as.formula(paste0(treated_var, " ~ 1")),
        data = match_candidates_i,
        replace = FALSE,
        estimand = "ATT",
        exact = matching_vars,
        m.order = "random",
        # verbose = TRUE,
        ratio = 1L # 1:1 match
      )[[1]]

    if(is.null(obj_matchit_i)) {
      message("Skipping trial - no matches found.")
      next
    }

    group_labels <- c("control_id" = 0, "treated_id" = 1)

    data_successfulmatches[[as.character(trialstart)]] <- dplyr::tibble(
      !!id_sym := match_candidates_i[[id_var]],
      matched = !is.na(obj_matchit_i$subclass),
      match_id = as.integer(as.character(obj_matchit_i$subclass)),
      !!treated_sym := obj_matchit_i$treat,
      weight = obj_matchit_i$weights,
      trialstart = trialstart,
    ) %>%
      dplyr::filter(matched) %>%
      dplyr::mutate(dplyr::across(
        !!treated_sym,
        \(x) factor(x, levels = unname(group_labels), labels = names(group_labels))
        )) %>%
      tidyr::pivot_wider(
        names_from = !!treated_sym,
        values_from=id,
      ) %>%
      dplyr::left_join(
        data_alltreated %>% dplyr::select(id, controlistreated = !!time_sym),
        by = c("control_id" = "id")
      )

    # ids for individuals matched as a control in all trials so far
    previouslymatched_ids <- c(
      previouslymatched_ids,
      data_successfulmatches[[as.character(trialstart)]]$control_id
      )

  }

  # currently a list with a dataset for each trial, bind into one dataset
  data_successfulmatches <- dplyr::bind_rows(data_successfulmatches)

  # check for duplicates in control_id and treated_id
  for (i in names(group_labels)) {
    if (any(duplicated(data_successfulmatches[[i]]))) {
      stop("Duplicate values in ", i)
    }
  }

  data_matched <- data_successfulmatches %>%
    tidyr::pivot_longer(
      cols = names(group_labels),
      names_to = "treated_tmp",
      values_to = id_var
    ) %>%
    dplyr::transmute(
      !!id_sym,
      match_id, #TODO create unique match_id, currently only unique within a trial
      !!treated_sym := treated_tmp == names(group_labels)[2],
      trialstart,
      controlistreated
    )

  out <- list(
    "data_matched" = data_matched,
    variable_names = list(
      "id_var" = id_var,
      "time_var" = time_var,
      "treated_var" = treated_var,
      "matching_vars" = matching_vars
    ),
    "data" = data
  )

  class(out) <- "seqtrial"

  return(out)

}
