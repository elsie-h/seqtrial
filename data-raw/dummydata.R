## code to prepare `dummydata` dataset

library(tidyverse)

sex_levs <- c("F", "M")

rbern <- function (n, p = 0.5) {
  stats::runif(n) > (1 - p)
}

set.seed(8765)
dat <- tibble(id = 1:500) %>%
  mutate(
    # number of months follow-up for each individual
    n = 24, #rpois(n=nrow(.), lambda = 10),
    # baseline demographics
    age = floor(rnorm(n = nrow(.), mean = 60, sd = 10)),
    sex = factor(
      sample(sex_levs, size = nrow(.), replace=TRUE),
      levels = sex_levs
    ),
    # baseline biomarker level
    biomarker = rnorm(n = nrow(.), mean = 1000, sd = 200),
    # peoples biomarkers typically increasing slightly over time
    biomarker_mult1 = rnorm(n = nrow(.), mean = 1.02, sd = 0.01)
  ) %>%
  uncount(n, .id = "time") %>%
  mutate(age = age + time%/%12) %>%
  mutate(biomarker_mult1 = rnorm(n = nrow(.), mean = biomarker_mult1, sd = 0.01)) %>%
  group_by(id) %>%
  mutate(across(biomarker, ~.x*cumprod(biomarker_mult1))) %>%
  ungroup() %>%
  # high values more likely to initiate treatment
  mutate(
    treated = if_else(
      biomarker > 1200,
      rbern(n = nrow(.), p=0.75),
      rbern(n = nrow(.), p=0.1)
    )
  ) %>%
  group_by(id) %>%
  mutate(across(treated, ~cumsum(.x) > 0)) %>%
  # create a lagged treatment column
  # FALSE in the last position replaces leading NAs with FALSE
  mutate(lag_treated = lag(treated, n=2, FALSE)) %>%
  ungroup() %>%
  # treatment reduces biomarker with 2 month lag
  mutate(
    biomarker_mult2 = if_else(
      lag_treated,
      rnorm(n = nrow(.), mean = 0.95, sd = 0.01),
      rnorm(n = nrow(.), mean = 1, sd = 0.01)
    )
  ) %>%
  group_by(id) %>%
  mutate(across(biomarker, ~.x*cumprod(biomarker_mult2))) %>%
  ungroup() %>%
  # people accumulate more risk of dying while biomarker high
  mutate(
    risk = if_else(
      biomarker > 1000,
      rbern(n = nrow(.), p=0.7),
      rbern(n = nrow(.), p=0.2)
    )
  ) %>%
  group_by(id) %>%
  # risk accumulates
  mutate(across(risk, ~cumsum(.x))) %>%
  # also people censored
  mutate(censor = rpois(n = 1, lambda = 30)) %>%
  mutate(status = risk > 3) %>%
  mutate(keep = cumsum(status) <= 1) %>%
  ungroup() %>%
  filter(keep & time <= censor)


# dat %>%
#   ggplot(aes(x = time, y = biomarker, colour = status, group = id)) +
#   geom_point(alpha = 0.2) +
#   geom_line(alpha = 0.2) +
#   facet_grid(rows = "treated")


dummydata <- dat %>%
  transmute(
    id, time, sex,
    age_grp = cut(
      age,
      breaks = c(18, seq(40,100,20)),
      include.lowest = TRUE,
      right = FALSE
    ),
    biomarker = cut(
      biomarker,
      breaks = c(200, seq(600, 1200, 200), 2000),
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4
    ),
    status
  )

usethis::use_data(dummydata, overwrite = TRUE)
