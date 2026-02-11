## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(kamp)
#library(devtools)
library(tidyverse)
library(spatstat.random)
#devtools::load_all()
set.seed(50)

## -----------------------------------------------------------------------------
univ_data <- sim_pp_data(lambda_n = 200, 
                    abundance = 0.3,
                    markvar1 = "immune",
                    markvar2 = "background",
                    distribution = "hom",
                    clust = FALSE)


## -----------------------------------------------------------------------------
as_tibble(univ_data) %>%
    ggplot(aes(x,y, color = marks)) +
    geom_point()

## -----------------------------------------------------------------------------
biv_data <- sim_pp_data_biv(lambda_n = 200, 
                    abundance = 0.3,
                    markvar1 = "immune1",
                    markvar2 = "immune2",
                    markvar3 = "background",
                    distribution = "inhom",
                    clust = TRUE)


## -----------------------------------------------------------------------------
as_tibble(biv_data) %>%
    ggplot(aes(x,y, color = marks)) +
    geom_point()

## -----------------------------------------------------------------------------
data(ovarian_df)
head(ovarian_df)

## -----------------------------------------------------------------------------
ids <- unique(ovarian_df$sample_id)
marksvar <- "immune"

for (id in ids) {
  df_sub <- ovarian_df %>% filter(sample_id == id)
  w <- convexhull.xy(df_sub$x, df_sub$y)
  pp_obj <- ppp(df_sub$x, df_sub$y, window = w, marks = df_sub[[marksvar]])
  
  p <- ggplot(as_tibble(pp_obj), aes(x, y, color = marks)) +
    geom_point(size = 0.6) +
    labs(title = paste("Sample:", id)) +
    theme_minimal()
  
  print(p)
}

## -----------------------------------------------------------------------------

# We can use the kamp function to calculate the KAMP CSR for the univariate data
# Lets first subset to one id

ids <- unique(ovarian_df$sample_id)
marksvar <- "immune"
univ_data <- ovarian_df %>% filter(sample_id == ids[1])
w <- convexhull.xy(univ_data$x, univ_data$y)
pp_univ_data <- ppp(univ_data$x, univ_data$y, window = w, marks = univ_data[[marksvar]])



univ_kamp <- kamp(pp_univ_data, 
                  rvals = c(0, 0.1, 0.2),
                  univariate = TRUE,
                  marksvar1 = "immune",
                  marksvar2 = "background")

univ_kamp

