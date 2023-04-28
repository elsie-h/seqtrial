# TODO needs updating

#' Calculate number of treated people who were matched per trial
#'
#' @param seqtrial_obj  An object with class seqtrial
#'
#' @return A data frame
#' @export
#'
#' @examples
#'
#' seqtrial_example <- seqtrial_matchit(
#' id_var ="id",
#' time_var = "time",
#' treated_var = "treated",
#' matching_vars = c("age_grp", "biomarker", "sex"),
#' data = dummydata
#' )
#' matchsuccess(seqtrial_example)
matchsuccess <- function(seqtrial_obj) {

  stopifnot("seqtrial_obj must have class \"seqtrial\"" = class(seqtrial_obj) == "seqtrial")

  id_sym <- rlang::sym(seqtrial_obj$variable_names$id_var)
  time_sym <- rlang::sym(seqtrial_obj$variable_names$time_var)
  treated_sym <- rlang::sym(seqtrial_obj$variable_names$treated_var)

  # get all the times at which people initiate treatment
  data_alltreated <- alltreated(id_sym, time_sym, treated_sym, seqtrial_obj$data)

  # created named character vector of variables to join by
  join_by_vars <- "trialstart"
  names(join_by_vars) <- seqtrial_obj$variable_names$time_var

  # calculate number of treated people who were matched per trial
  matchsuccess <- data_alltreated %>%
    dplyr::group_by(!!time_sym) %>%
    dplyr::count(name="n_treated") %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
      seqtrial_obj$data_matched %>%
        dplyr::filter(!!treated_sym) %>%
        dplyr::group_by(trialstart) %>%
        dplyr::count(name="n_matched") %>%
        dplyr::ungroup(),
      by = join_by_vars
    ) %>%
    dplyr::mutate(
      n_matched = tidyr::replace_na(n_matched, 0),
      n_unmatched = n_treated - n_matched
    )

  return(matchsuccess)

}
