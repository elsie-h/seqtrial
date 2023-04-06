# create function to catch errors
library(rlang)
library(tidyr)
library(purrr)
library(dplyr)
library(MatchIt)

# people initial treatment once then stay on it for the rest of follow-up

# data<-dummydata
# id_var<-"id"
# time_var<-"time"
# treated_var<-"treated"
# treatment_var<-NULL
# matching_vars<-c("age_grp", "biomarker", "sex")

# TODO extend to datasets comparing two treatments, not just treatment vs no treatment

seqtrial <- function(data,
                     id_var,
                     time_var,
                     treated_var,
                     matching_vars) {

  # check data is in correct format
  stopifnot("`data` must have class \"data.frame\"" = "data.frame" %in% class(data))
  # check supplied variables are in `data`
  all_vars <- c(id_var, time_var, treated_var, matching_vars)
  stopifnot("All supplied vars must be in `data`" = all(all_vars %in% names(data)))
  # check time_var is numeric (double or integer)
  # TODO (extend to handle dates)
  stopifnot("`time_var` must be numeric" = is.numeric(data[[time_var]]))
  # check treated_var is logical
  stopifnot("`treated_var` must be logical" = is.logical(data[[treated_var]]))
  # check matching_vars are factors (for now as only exact matching, but will extend to other matching techniques)
  check_factors <- data %>% select(all_of(matching_vars)) %>% map_lgl(\(x) is.factor(x))
  if (!all(check_factors)) {
    stop(
      paste0(
        paste(names(check_factors)[!check_factors], collapse = " and "),
        " are not factors"
      )
    )
  }
  # only keep the specified variables
  data <- data %>% select(all_of(all_vars))
  # check for missing values
  check_missing <- data %>% map_lgl(\(x) any(is.na(x)))
  if (any(check_missing)) {
    stop(
      paste0(
        paste(names(check_missing)[check_missing], collapse = " and "),
        " have missing data"
      )
    )
  }

  # create symbols for programming
  id_sym <- sym(id_var)
  time_sym <- sym(time_var)
  treated_sym <- sym(treated_var)

  # get all the times at which people initiate treatment
  alltreated <- data %>%
    filter(!!treated_sym) %>%
    arrange(!!id_sym, !!time_sym) %>%
    group_by(!!id_sym) %>%
    slice(1) %>%
    arrange(!!time_sym)

  # create empty objects for output
  data_successfulmatches <- list()
  previouslymatched_ids <- vector(mode = mode(data[[id_var]]), length=0)

  message("Sequential trials matching report:")

  for (trialstart in unique(alltreated[[time_var]])) {

    message("---- ", time_var, " = ", trialstart, " ----")

    treated_ids_i <- alltreated %>% filter(!!time_sym == trialstart) %>% pull(!!id_sym)

    data_i <- data %>%
      # only keep rows where time_var = trialstart
      filter(!!time_sym == trialstart)

    # data for match candidates
    match_candidates_i <- bind_rows(
      # treated who initiated treatment on trialstart
      # (people can be treated if previously matched as a control)
      data_i %>%
        filter(!!id_sym %in% treated_ids_i),
      # control
      data_i %>%
        filter(!(!!treated_sym) & !(!!id_sym %in% previouslymatched_ids))
    )

    # TODO in future allow user to specify matching options in seqtrial()
    # run match algorithm
    safely_matchit <- purrr::safely(matchit)
    obj_matchit_i <-
      safely_matchit(
        formula = as.formula(paste0(treated_var, " ~ 1")),
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

    data_successfulmatches[[as.character(trialstart)]] <- tibble(
      !!id_sym := match_candidates_i[[id_var]],
      matched = !is.na(obj_matchit_i$subclass),
      match_id = as.integer(as.character(obj_matchit_i$subclass)),
      !!treated_sym := obj_matchit_i$treat,
      weight = obj_matchit_i$weights,
      trialstart = trialstart,
    ) %>%
      filter(matched) %>%
      mutate(across(
        !!treated_sym,
        \(x) factor(x, levels = unname(group_labels), labels = names(group_labels))
        )) %>%
      pivot_wider(
        names_from = !!treated_sym,
        values_from=id,
      ) %>%
      left_join(
        alltreated %>% select(id, controlistreated = !!time_sym),
        by = c("control_id" = "id")
      )

    # ids for individuals matched as a control in all trials so far
    previouslymatched_ids <- c(previouslymatched_ids, data_successfulmatches[[as.character(trialstart)]]$control_id)

  }

  # currently a list with a dataset for each trial, bind into one dataset
  data_successfulmatches <- bind_rows(data_successfulmatches)

  # check for duplicates in control_id and treated_id
  for (i in names(group_labels)) {
    if (any(duplicated(data_successfulmatches[[i]]))) {
      stop("Duplicate values in ", i)
    }
  }

  # print matching success
  matchingsuccess <- alltreated %>%
    group_by(!!time_sym) %>%
    count(name="n_treated") %>%
    ungroup() %>%
    left_join(
      data_successfulmatches %>%
        group_by(trialstart) %>%
        count(name="n_matched") %>%
        ungroup() %>%
        rename(!!time_sym := trialstart),
      by = time_var
    ) %>%
    mutate(
      n_matched = replace_na(n_matched, 0),
      n_unmatched = n_treated - n_matched
    )

  data_matched <- data_successfulmatches %>%
    pivot_longer(
      cols = names(group_labels),
      names_to = "treated_tmp",
      values_to = id_var
    ) %>%
    transmute(
      !!id_sym,
      match_id, #TODO create unique match_id, currently only unique within a trial
      !!treated_sym := treated_tmp == names(group_labels)[2],
      trialstart,
      controlistreated
    )

  out <- list(
    "matchingsuccess" = matchingsuccess,
    "data_matched" = data_matched
  )

  return(out)

}

# seqtrial(
#   data = dummydata,
#   id_var ="id",
#   time_var = "time",
#   treated_var = "treated",
#   matching_vars = c("age_grp", "biomarker", "sex")
# )
