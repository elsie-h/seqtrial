## code to prepare `dummydata` dataset

library(tidyverse)

# TODO
# simulate better data
# for now the dummy data is just so that the sequential trial code runs and finds some matches
# survival times need updated to introduce associations with characteristics/treatment

sex_levs <- c("F", "M")

rbern <- function (n, p = 0.5) {
  stats::runif(n) > (1 - p)
}

set.seed(8765)
dummydata0 <- tibble(id = 101:600) %>%
  mutate(
    # max months follow-up for each individual
    n = 24,
    # baseline demographics
    age = floor(rnorm(n = nrow(.), mean = 60, sd = 10)),
    sex = factor(
      sample(sex_levs, size = nrow(.), replace=TRUE),
      levels = sex_levs
    ),
    # baseline biomarker level
    biomarker = rnorm(n = nrow(.), mean = 1000, sd = 200),
    # peoples biomarkers typically increasing slightly over time
    biomarker_mult1 = rnorm(n = nrow(.), mean = 1.02, sd = 0.01),
    # if the individual is treated, which treatment
    # this is for sequential trials that compare two treatments rather
    # than treated vs untreated
    treatment = sample(x = c("A", "B"), replace = TRUE, size = nrow(.))
  ) %>%
  uncount(n, .id = "time") %>%
  # don't bother updating age
  # mutate(age = age + time%/%12) %>%
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
  mutate(across(treatment, ~if_else(treated, .x, NA_character_))) %>%
  # treatment starts reducing biomarker with 2 month lag
  mutate(
    biomarker_mult2 = if_else(
      lag_treated,
      rnorm(n = nrow(.), mean = 0.95, sd = 0.01),
      rnorm(n = nrow(.), mean = 1, sd = 0.01)
    )
  ) %>%
  group_by(id) %>%
  mutate(across(biomarker, ~.x*cumprod(biomarker_mult2))) %>%
  ungroup()

dummydata1 <- dummydata0 %>%
  # TODO
  # for now, survival times unrelated to characteristics
  # update this later to introduce associations
  left_join(
    dummydata0 %>%
      distinct(id) %>%
      mutate(survtime = ceiling(rweibull(n=nrow(.), shape=5, scale=24))),
    by = "id"
  ) %>%
  filter(time <= survtime) %>%
  mutate(death = if_else(time==survtime, TRUE, FALSE))

# dat %>%
#   ggplot(aes(x = time, y = biomarker, colour = status, group = id)) +
#   geom_point(alpha = 0.2) +
#   geom_line(alpha = 0.2) +
#   facet_grid(rows = "treated")

dummydata <- dummydata1 %>%
  transmute(
    id, time, sex, age, biomarker,
    age_grp = cut(
      age,
      breaks = c(18, seq(40,100,20)),
      include.lowest = TRUE,
      right = FALSE
    ),
    biomarker_grp = cut(
      biomarker,
      breaks = c(200, seq(600, 1200, 200), 2000),
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4
    ),
    treated, treatment,
    death
  )

usethis::use_data(dummydata, overwrite = TRUE)
