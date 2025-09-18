suppressPackageStartupMessages({
  librarian::shelf(
    dplyr, tidyr, ggplot2, readr, stringr, tibble, purrr,
    unmarked, vegan, betapart, ggvenn, lme4, emmeans
  )
})

# Petites aides
`%||%` <- function(x, y) if (is.null(x)) y else x
zstd <- function(v){ m <- mean(v,na.rm=TRUE); s <- sd(v,na.rm=TRUE); z <- (v-m)/ifelse(s>0,s,1); z[is.na(z)] <- 0; z }
logit <- function(p) log(p/(1-p))
