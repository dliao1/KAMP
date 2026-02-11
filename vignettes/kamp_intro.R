## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(KAMP)
#library(devtools)
library(tidyverse)
library(spatstat.random)
#devtools::load_all()
set.seed(50)

## -----------------------------------------------------------------------------
data(ovarian_df)
head(ovarian_df)

## -----------------------------------------------------------------------------
ids <- unique(ovarian_df$sample_id)
mark_var <- "immune"

for (id in ids) {
  df_sub <- ovarian_df %>% filter(sample_id == id)
  w <- convexhull.xy(df_sub$x, df_sub$y)
  pp_obj <- ppp(df_sub$x, df_sub$y, window = w, marks = df_sub[[mark_var]])
  
  p <- ggplot(as_tibble(pp_obj), aes(x, y, color = marks)) +
    geom_point(size = 0.6) +
    labs(title = paste("Sample:", id)) +
    theme_minimal()
  
  print(p)
}

## -----------------------------------------------------------------------------
ids <- unique(ovarian_df$sample_id)
univ_data <- ovarian_df %>% filter(sample_id == ids[1])
mark_var <- "immune"
head(univ_data)

## -----------------------------------------------------------------------------
univ_kamp <- kamp(univ_data, 
                  rvals = seq(0, 100, by = 10),
                  univariate = TRUE,
                  mark_var = mark_var,
                  mark1 = "immune")

univ_kamp

## -----------------------------------------------------------------------------
univ_kamp %>%
  ggplot(aes(x = r)) +
  geom_line(aes(y = theo_csr, color = "theo_csr", linetype = "theo_csr"), size = 1) +
  geom_line(aes(y = kamp_csr, color = "kamp_csr", linetype = "kamp_csr"), size = 1) +
  geom_line(aes(y = k, color = "k", linetype = "k"), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(
    values = c(
      "theo_csr" = "black",
      "kamp_csr" = "blue",
      "k" = "red"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "theo_csr" = "solid",
      "kamp_csr" = "solid",
      "k" = "dotted"
    )
  ) +
  labs(
    title = "Univariate KAMP Expectation",
    x = "r",
    y = "Value",
    color = "Series",
    linetype = "Series"
  ) +
  theme_minimal()


## -----------------------------------------------------------------------------
univ_kamp %>%
  ggplot(aes(x = r)) +
  geom_line(aes(y = k - kamp_csr, color = "k - kamp_csr", linetype = "k - kamp_csr"), size = 1) +
  geom_line(aes(y = k - theo_csr, color = "k - theo_csr", linetype = "k - theo_csr"), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(
    values = c(
      "k - theo_csr" = "red",
      "k - kamp_csr" = "blue"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "k - theo_csr" = "solid",
      "k - kamp_csr" = "solid"
    )
  ) +
  labs(
    title = "Univariate KAMP Expectation - differences",
    x = "r",
    y = "Value",
    color = "Series",
    linetype = "Series"
  ) +
  theme_minimal()

## -----------------------------------------------------------------------------
univ_kamp_var <- kamp(univ_data,
                      rvals = seq(0, 100, by = 10),
                      univariate = TRUE,
                      mark_var = mark_var,
                      mark1 = "immune",
                      variance = TRUE)
univ_kamp_var

## -----------------------------------------------------------------------------
#devtools::load_all()

univ_kamp_var %>%
  ggplot(aes(x = r, y = var)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  labs(title = "Univariate KAMP Variance", x = "r", y = "Variance") +
  theme_minimal()

## -----------------------------------------------------------------------------
ids <- unique(ovarian_df$sample_id)
biv_data <- ovarian_df %>% 
  filter(sample_id == ids[1]) #%>%
#filter(phenotype %in% c("helper t cells", "cytotoxic t cells", "other")) %>%
#droplevels()

mark_var <- "phenotype"

head(biv_data)

## -----------------------------------------------------------------------------

biv_kamp <- kamp(df = biv_data,
                 rvals = seq(0, 100, by = 10),
                 univariate = FALSE,
                 mark_var = mark_var,
                 mark1 = "helper t cell",
                 mark2 = "cytotoxic t cell")
head(biv_kamp)

## -----------------------------------------------------------------------------
biv_kamp %>%
  ggplot(aes(x = r)) +
  geom_line(aes(y = theo_csr, color = "theo_csr", linetype = "theo_csr"), size = 1) +
  geom_line(aes(y = kamp, color = "kamp", linetype = "kamp"), size = 1) +
  geom_line(aes(y = k, color = "k", linetype = "k"), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(
    values = c(
      "theo_csr" = "black",
      "kamp" = "blue",
      "k" = "red"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "theo_csr" = "solid",
      "kamp" = "dotted",
      "k" = "dotted"
    )
  ) +
  labs(
    title = "Bivariate KAMP Expectation",
    x = "r",
    y = "Value",
    color = "Series",
    linetype = "Series"
  ) +
  theme_minimal()

## -----------------------------------------------------------------------------
biv_kamp_var <- kamp(df = biv_data,
                     rvals = seq(0, 100, by = 10),
                     univariate = FALSE,
                     mark_var = mark_var,
                     mark1 = "helper t cell",
                     mark2 = "cytotoxic t cell",
                     variance = TRUE)
head(biv_kamp_var)

## -----------------------------------------------------------------------------
biv_kamp_var %>%
  ggplot(aes(x = r, y = var)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  labs(title = "Bivariate KAMP Variance", x = "r", y = "Variance") +
  theme_minimal()

## -----------------------------------------------------------------------------
ids <- unique(ovarian_df$sample_id)
mark_var <- "immune"
univ_data <- ovarian_df %>% filter(sample_id == ids[1])

univ_kamp_lite <- kamp(df = univ_data,
                       rvals = seq(0, 100, by = 10),
                       univariate = TRUE,
                       mark_var = mark_var,
                       mark1 = "immune",
                       thin = TRUE,
                       p_thin = 0.3)
univ_kamp_lite

## -----------------------------------------------------------------------------
univ_kamp_lite %>%
  ggplot(aes(x = r)) +
  geom_line(aes(y = theo_csr, color = "theo_csr", linetype = "theo_csr"), size = 1) +
  geom_line(aes(y = kamp, color = "kamp", linetype = "kamp"), size = 1) +
  geom_line(aes(y = k, color = "k", linetype = "k"), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(
    values = c(
      "theo_csr" = "black",
      "kamp" = "blue",
      "k" = "red"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "theo_csr" = "solid",
      "kamp" = "dotted",
      "k" = "dotted"
    )
  ) +
  labs(
    title = "Univariate KAMP-lite Expectation",
    x = "r",
    y = "Value",
    color = "Series",
    linetype = "Series"
  ) +
  theme_minimal()

## -----------------------------------------------------------------------------
univ_kamp_lite_var <- kamp(df = univ_data,
                           rvals = seq(0, 100, by = 10),
                           univariate = TRUE,
                           mark_var = mark_var,
                           mark1 = "immune",
                           thin = TRUE,
                           p_thin = 0.3,
                           variance = TRUE) # should display a warning message
univ_kamp_lite_var

## -----------------------------------------------------------------------------
ids <- unique(ovarian_df$sample_id)
biv_data <- ovarian_df %>% 
  filter(sample_id == ids[1]) #%>%
#filter(phenotype %in% c("helper t cell", "cytotoxic t cell", "other")) %>%
#droplevels()

mark_var <- "phenotype"

head(biv_data)

## -----------------------------------------------------------------------------
biv_kamp_lite <- kamp(df = biv_data,
                      rvals = seq(0, 100, by = 10),
                      univariate = FALSE,
                      mark_var = mark_var,
                      mark1 = "helper t cell",
                      mark2 = "cytotoxic t cell",
                      thin = TRUE,
                      p_thin = 0.3)
head(biv_kamp_lite)

## -----------------------------------------------------------------------------
biv_kamp_lite %>%
  ggplot(aes(x = r)) +
  geom_line(aes(y = theo_csr, color = "theo_csr", linetype = "theo_csr"), size = 1) +
  geom_line(aes(y = kamp, color = "kamp", linetype = "kamp"), size = 1) +
  geom_line(aes(y = k, color = "k", linetype = "k"), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(
    values = c(
      "theo_csr" = "black",
      "kamp" = "blue",
      "k" = "red"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "theo_csr" = "solid",
      "kamp" = "dotted",
      "k" = "dotted"
    )
  ) +
  labs(
    title = "Bivariate KAMP-lite Expectation",
    x = "r",
    y = "Value",
    color = "Series",
    linetype = "Series"
  ) +
  theme_minimal()

## -----------------------------------------------------------------------------
biv_kamp_lite_var <- kamp(df = biv_data,
                          rvals = seq(0, 100, by = 10),
                          univariate = FALSE,
                          mark_var = mark_var,
                          mark1 = "helper t cell",
                          mark2 = "cytotoxic t cell",
                          variance = TRUE,
                          thin = TRUE,
                          p_thin = 0.3) # should display a warning message
head(biv_kamp_lite_var)

## -----------------------------------------------------------------------------
biv_kamp_lite_var %>%
  ggplot(aes(x = r, y = var)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  labs(title = "Bivariate KAMP lite Variance", x = "r", y = "Variance") +
  theme_minimal()

