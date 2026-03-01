suppressPackageStartupMessages(library(KAMP))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(spatstat.random))
suppressPackageStartupMessages(library(microbenchmark))
suppressPackageStartupMessages(library(spatstat.geom))
suppressPackageStartupMessages(library(spatstat.explore))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(here))


source(here::here("vignettes", "kamp_sim_helpers.R"))

# objective here is to compare the runtime of Rcpp true vs Rcpp false for variance and also with Kperm :)

n_values <- 500
abundance_values <- 0.3
distribution_values <- "inhom"
clust_values <- c(TRUE)
correction <- c("trans", "iso")
univariate <- FALSE
seed_start = 500
n_rep <- 1 # just in case for testing!
nperm <- 5


param_grid <- expand.grid(n = n_values,
                          abundance = abundance_values,
                          distribution = distribution_values,
                          clust = clust_values,
                          correction = correction,
                          univariate = univariate,
                          rep = seq_len(n_rep),
                          seed_start = seed_start,
                          stringsAsFactors = FALSE)

wd = getwd()

if(substring(wd, 4, 8) == "Users"){
  doLocal = TRUE
}else{
  doLocal = FALSE
}

if (doLocal == TRUE) {
  num_cores <- 8
  cl <- makeCluster(num_cores)

  clusterEvalQ(cl, {
    library(dplyr)
    library(tibble)
    library(purrr)
    library(spatstat)
  })

  clusterExport(
    cl,
    varlist = c(
      "param_grid",
      "sim_pp_data_biv",
      "k_expectation",
      "kinhom_expectation",
      "kamp",
      "perm_expectation"
    ),
    envir = environment()
  )


  results_list <- parLapply(cl, 1:nrow(param_grid), function(i) {
    params <- param_grid[i, ]

    sim_data <- sim_pp_data_biv(
      lambda_n = params$n,
      abundance = params$abundance,
      markvar1 = "immune1",
      markvar2 = "immune2",
      markvar3 = "background",
      distribution = params$distribution,
      clust = params$clust
    )

    # EXPECTATION
    t_k <- system.time({
      k_result <- k_expectation(
        sim_data,
        rvals = seq(0, 0.5, by = 0.05),
        mark1 = "immune1",
        mark2 = "immune2",
        correction = params$correction,
        univariate = params$univariate
      )
    })

    t_kinhom <- system.time({
      kinhom_result <- kinhom_expectation(
        sim_data,
        rvals = seq(0, 0.5, by = 0.05),
        mark1 = "immune1",
        mark2 = "immune2",
        correction = params$correction,
        univariate = params$univariate
      )
    })


    # VARIANCE

    # KAMP
    t_kamp_Rcpp_true_var <- system.time({
      kamp_Rcpp_true_var_result <- kamp(
        sim_data,
        rvals = seq(0, 0.5, by = 0.05),
        mark_var = "marks",
        mark1 = "immune1",
        mark2 = "immune2",
        correction = params$correction,
        univariate = params$univariate,
        thin = FALSE,
        background = "background",
        Rcpp = TRUE,
        variance = TRUE
      )
    })

    t_kamp_Rcpp_false_var <- system.time({
      kamp_Rcpp_false_var_result <- kamp(
        sim_data,
        rvals = seq(0, 0.5, by = 0.05),
        mark_var = "marks",
        mark1 = "immune1",
        mark2 = "immune2",
        correction = params$correction,
        univariate = params$univariate,
        thin = FALSE,
        background = "background",
        Rcpp = FALSE,
        variance = TRUE
      )
    })

    # KAMP LITE
    t_kamp_lite_Rcpp_true_var <- system.time({
      kamp_lite_Rcpp_true_var_result <- kamp(
        sim_data,
        rvals = seq(0, 0.5, by = 0.05),
        mark_var = "marks",
        mark1 = "immune1",
        mark2 = "immune2",
        correction = params$correction,
        univariate = params$univariate,
        thin = TRUE,
        p_thin = 0.5,
        background = "background",
        Rcpp = TRUE,
        variance = TRUE
      )
    })

    t_kamp_lite_Rcpp_false_var <- system.time({
      kamp_lite_Rcpp_false_var_result <- kamp(
        sim_data,
        rvals = seq(0, 0.5, by = 0.05),
        mark_var = "marks",
        mark1 = "immune1",
        mark2 = "immune2",
        correction = params$correction,
        univariate = params$univariate,
        thin = TRUE,
        p_thin = 0.5,
        background = "background",
        Rcpp = FALSE,
        variance = TRUE
      )
    })

    # PERM
    t_kperm_var <- system.time({
      kperm_var_result <- perm_variance(
        sim_data,
        rvals = seq(0, 0.5, by = 0.05),
        mark1 = "immune1",
        mark2 = "immune2",
        univariate = params$univariate,
        nperm = nperm
      )
    })

    list(
      results = list(
        k = k_result,
        kinhom = kinhom_result,
        kamp_Rcpp_true_var = kamp_Rcpp_true_var_result,
        kamp_lite_Rcpp_true_var = kamp_lite_Rcpp_true_var_result,
        kamp_Rcpp_false_var = kamp_Rcpp_false_var_result,
        kamp_lite_Rcpp_false_var = kamp_lite_Rcpp_false_var_result,
        kperm_var = kperm_var_result
      ),
      times = tibble(
        n = params$n,
        abundance = params$abundance,
        clust = params$clust,
        correction = params$correction,
        rep = params$rep,
        method = c("K",
                   "Kinhom",
                   "KAMP_Rcpp_true_var",
                   "KAMP_lite_Rcpp_true_var",
                   "KAMP_Rcpp_false_var",
                   "KAMP_lite_Rcpp_false_var",
                   "Kperm_var"),
        elapsed = c(
          t_k["elapsed"],
          t_kinhom["elapsed"],
          t_kamp_Rcpp_true_var["elapsed"],
          t_kamp_lite_Rcpp_true_var["elapsed"],
          t_kamp_Rcpp_false_var["elapsed"],
          t_kamp_lite_Rcpp_false_var["elapsed"],
          t_kperm_var["elapsed"]
        )
      )
    )
  }
  )

  stopCluster(cl)


} else { # this is on the cluster
  scenario <- as.numeric(commandArgs(trailingOnly=TRUE))

  params <- param_grid[scenario, ]
  print(scenario)

  SEED.START <- params$seed_start

  # for loop looping through all the reps
  results_list = vector("list", length = n_rep)
  for (rep in 1:n_rep) {

    seed.iter = (SEED.START - 1)*n_rep + rep
    print(seed.iter)
    set.seed(seed.iter)

    sim_data <- sim_pp_data_biv(
      lambda_n = params$n,
      abundance = params$abundance,
      markvar1 = "immune1",
      markvar2 = "immune2",
      markvar3 = "background",
      distribution = params$distribution,
      clust = params$clust
    )

    # EXPECTATION
    t_k <- system.time({
      k_result <- k_expectation(
        sim_data,
        rvals = seq(0, 0.5, by = 0.05),
        mark1 = "immune1",
        mark2 = "immune2",
        correction = params$correction,
        univariate = params$univariate
      )
    })

    t_kinhom <- system.time({
      kinhom_result <- kinhom_expectation(
        sim_data,
        rvals = seq(0, 0.5, by = 0.05),
        mark1 = "immune1",
        mark2 = "immune2",
        correction = params$correction,
        univariate = params$univariate
      )
    })

    # VARIANCE

    #KAMP
    t_kamp_Rcpp_true_var <- system.time({
      kamp_Rcpp_true_var_result <- kamp(
        sim_data,
        rvals = seq(0, 0.5, by = 0.05),
        mark_var = "marks",
        mark1 = "immune1",
        mark2 = "immune2",
        correction = params$correction,
        univariate = params$univariate,
        thin = FALSE,
        background = "background",
        Rcpp = TRUE,
        variance = TRUE
      )
    })

    t_kamp_Rcpp_false_var <- system.time({
      kamp_Rcpp_false_var_result <- kamp(
        sim_data,
        rvals = seq(0, 0.5, by = 0.05),
        mark_var = "marks",
        mark1 = "immune1",
        mark2 = "immune2",
        correction = params$correction,
        univariate = params$univariate,
        thin = FALSE,
        background = "background",
        Rcpp = FALSE,
        variance = TRUE
      )
    })

    # KAMP LITE
    t_kamp_lite_Rcpp_true_var <- system.time({
      kamp_lite_Rcpp_true_var_result <- kamp(
        sim_data,
        rvals = seq(0, 0.5, by = 0.05),
        mark_var = "marks",
        mark1 = "immune1",
        mark2 = "immune2",
        correction = params$correction,
        univariate = params$univariate,
        thin = TRUE,
        p_thin = 0.5,
        background = "background",
        Rcpp = TRUE,
        variance = TRUE
      )
    })

    t_kamp_lite_Rcpp_false_var <- system.time({
      kamp_lite_Rcpp_false_var_result <- kamp(
        sim_data,
        rvals = seq(0, 0.5, by = 0.05),
        mark_var = "marks",
        mark1 = "immune1",
        mark2 = "immune2",
        correction = params$correction,
        univariate = params$univariate,
        thin = TRUE,
        p_thin = 0.5,
        background = "background",
        Rcpp = FALSE,
        variance = TRUE
      )
    })



    # PERM
    t_kperm_var <- system.time({
      kperm_var_result <- perm_variance(
        sim_data,
        rvals = seq(0, 0.5, by = 0.05),
        mark1 = "immune1",
        mark2 = "immune2",
        univariate = params$univariate,
        nperm = nperm
      )
    })

    results_list[[rep]] <- list(
      results = list(
        k = k_result,
        kinhom = kinhom_result,
        kamp_Rcpp_true_var = kamp_Rcpp_true_var_result,
        kamp_lite_Rcpp_true_var = kamp_lite_Rcpp_true_var_result,
        kamp_Rcpp_false_var = kamp_Rcpp_false_var_result,
        kamp_lite_Rcpp_false_var = kamp_lite_Rcpp_false_var_result,
        kperm_var = kperm_var_result
      ),
      times = tibble(
        n = params$n,
        abundance = params$abundance,
        clust = params$clust,
        correction = params$correction,
        rep = rep,
        method = c("K",
                   "Kinhom",
                   "KAMP_Rcpp_true_var",
                   "KAMP_lite_Rcpp_true_var",
                   "KAMP_Rcpp_false_var",
                   "KAMP_lite_Rcpp_false_var",
                   "Kperm_var"),
        elapsed = c(
          t_k["elapsed"],
          t_kinhom["elapsed"],
          t_kamp_Rcpp_true_var["elapsed"],
          t_kamp_lite_Rcpp_true_var["elapsed"],
          t_kamp_Rcpp_false_var["elapsed"],
          t_kamp_lite_Rcpp_false_var["elapsed"],
          t_kperm_var["elapsed"]
        )
      )
    )

  }

}

filename = paste0(here::here("vignettes", "output_Rcpp"), "/", scenario, ".RDA")
save(results_list, file = filename)
