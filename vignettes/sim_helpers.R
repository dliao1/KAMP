k_expectation <- function(ppp_obj,
                          rvals,
                          mark1,
                          mark2 = NULL,
                          correction = "trans",
                          univariate = TRUE) {
  k <- NULL
  if (univariate == TRUE) {
    k = Kcross(ppp_obj, i = mark1, j = mark1,
               r = rvals,
               correction = correction) %>%
      as_tibble() %>%
      rename(k = correction) %>%
      mutate(method = "k",
             fundiff = k - theo,
             correction = correction) %>%
      select(r, theo_csr = theo, k, fundiff, method, correction)

  } else {
    k = Kcross(ppp_obj,
               i = mark1,
               j = mark2,
               r = rvals,
               correction = correction) %>%
      as_tibble() %>%
      rename(k = correction) %>%
      mutate(method = "k",
             fundiff = k - theo,
             correction = correction) %>%
      select(r, theo_csr = theo, k, fundiff, method, correction)
  }

  return(k)

}


kinhom_expectation <- function(ppp_obj,
                               rvals,
                               mark1,
                               mark2 = NULL,
                               correction = "trans",
                               univariate = TRUE) {
  kinhom <- NULL
  if (univariate == TRUE) {
    kinhom = Kinhom(subset(ppp_obj, marks == mark1),
                    r = rvals,
                    correction = correction) %>%
      as_tibble() %>%
      rename(kinhom = correction) %>%
      mutate(method = "kinhom",
             fundiff = kinhom - theo,
             correction = correction) %>%
      select(r, theo_csr = theo, kinhom, fundiff, method, correction)

  } else {
    kinhom = Kcross.inhom(ppp_obj,
                          i = mark1,
                          j = mark2,
                          r = rvals,
                          correction = correction) %>%
      as_tibble() %>%
      rename(kinhom = correction) %>%
      mutate(method = "kinhom",
             fundiff = kinhom - theo,
             correction = correction) %>%
      select(r, theo_csr = theo, kinhom, fundiff, method, correction)
  }

  return(kinhom)
}

perm_expectation <- function(ppp_obj,
                             rvals,
                             mark1,
                             mark2 = NULL,
                             correction = "trans",
                             univariate = TRUE,
                             nperm = 10) {
  kperm <- NULL
  if (univariate == TRUE) {

    k = Kcross(ppp_obj, i = mark1, j = mark1,
               r = rvals,
               correction = correction)

    k = k %>%
      as_tibble() %>%
      rename(correction_k = correction) %>%
      mutate(method = "k",
             fundiff = correction_k - theo,
             initial_correction = correction) %>%
      select(r, theo_csr = theo, correction_k, fundiff, method, initial_correction)

    kf = function(obj){
      kdf = Kcross(i = mark1,
                   j = mark1,
                   r = rvals,
                   correction = correction)

      as_tibble(kdf) %>%
        filter(r %in% rvals) %>%
        rename(correction_k = correction) %>%
        select(r, correction_k)
    }

    perms = rlabel(ppp_obj, nsim = nperm)
    kperm = map_dfr(perms, kf)
    kperm = kperm %>%
      group_by(r) %>%
      summarise(csr = mean(correction_k)) %>% # empirical CSR
      ungroup() %>%
      mutate(correction_k = k$correction_k,
             fundiff = correction_k - csr,
             method = "perm",
             initial_correction = k$initial_correction,
             theo_csr = k$theo_csr) %>%
      select(r, theo_csr, kperm_csr = csr, correction_k, fundiff, method, initial_correction)

  } else {

    k = Kcross(ppp_obj,
               i = mark1,
               j = mark2,
               r = rvals,
               correction = correction)

    k = k %>%
      as_tibble() %>%
      rename(correction_k = correction) %>%
      mutate(method = "k",
             fundiff = correction_k - theo,
             initial_correction = correction) %>%
      select(r, theo_csr = theo, correction_k, fundiff, method, initial_correction)

    kf = function(obj){
      kdf = Kcross(obj,
                   i = mark1,
                   j = mark2,
                   r = rvals,
                   correction = correction)

      as_tibble(kdf) %>%
        filter(r %in% rvals) %>%
        rename(correction_k = correction) %>%
        select(r, correction_k)
    }

    perms = rlabel(ppp_obj, nsim = nperm)
    kperm = map_dfr(perms, kf)
    kperm = kperm %>%
      group_by(r) %>%
      summarise(csr = mean(correction_k)) %>% # empirical CSR
      ungroup() %>%
      mutate(correction_k = k$correction_k,
             fundiff = correction_k - csr,
             method = "perm",
             initial_correction = k$initial_correction,
             theo_csr = k$theo_csr) %>%
      select(r, theo_csr, kperm_csr = csr, k = correction_k, fundiff, method, initial_correction)

  }

  return(kperm)
}

perm_variance <- function(ppp_obj,
                                  rvals,
                                  mark1,
                                  mark2 = NULL,
                                  correction = "trans",
                                  univariate = TRUE,
                                  nperm = 10) {
  kperm <- NULL
  if (univariate == TRUE) {

    k = Kcross(ppp_obj, i = mark1, j = mark1,
               r = rvals,
               correction = correction)

    k = k %>%
      as_tibble() %>%
      rename(correction_k = correction) %>%
      mutate(method = "k",
             fundiff = correction_k - theo,
             initial_correction = correction) %>%
      select(r, theo_csr = theo, correction_k, fundiff, method, initial_correction)

    kf = function(obj){
      kdf = Kcross(i = mark1,
                   j = mark1,
                   r = rvals,
                   correction = correction)

      as_tibble(kdf) %>%
        filter(r %in% rvals) %>%
        rename(correction_k = correction) %>%
        select(r, correction_k) %>%
        mutate(khat = k$correction_k)
    }

    perms = rlabel(ppp_obj, nsim = nperm)
    kperm = map_dfr(perms, kf)
    kperm = kperm %>%
      group_by(r) %>%
      summarise(csr = mean(correction_k),
                var = var(correction_k),
                pvalue = sum(correction_k >= khat)/nperm
                ) %>% # empirical CSR
      ungroup() %>%
      mutate(correction_k = k$correction_k,
             fundiff = correction_k - csr,
             method = "perm",
             initial_correction = k$initial_correction,
             theo_csr = k$theo_csr) %>%
      select(r, theo_csr, kperm_csr = csr, correction_k, fundiff, method, initial_correction)

    res = kperm %>%
      mutate(khat = rep(k$correction_k, times = 1))


    res2 = res %>%
      mutate(Z = (khat - kperm_csr) / sqrt(var), # kperm_csr = kperm expectation/new csr
             pvalue = pnorm(-Z),
             method = "kperm approx")

    kperm = bind_rows(res, res2) %>%
      mutate(fundiff = khat - kperm_csr) %>%
      select(r, theo_csr, kperm_csr, khat, fundiff, var, Z, pvalue, method, initial_correction)

  } else {

    k = Kcross(ppp_obj,
               i = mark1,
               j = mark2,
               r = rvals,
               correction = correction)

    k = k %>%
      as_tibble() %>%
      rename(correction_k = correction) %>%
      mutate(method = "k",
             fundiff = correction_k - theo,
             initial_correction = correction) %>%
      select(r, theo_csr = theo, correction_k, fundiff, method, initial_correction)

    kf = function(obj){
      kdf = Kcross(obj,
                   i = mark1,
                   j = mark2,
                   r = rvals,
                   correction = correction)

      as_tibble(kdf) %>%
        filter(r %in% rvals) %>%
        rename(correction_k = correction) %>%
        select(r, correction_k) %>%
        mutate(khat = k$correction_k)
    }

    perms = rlabel(ppp_obj, nsim = nperm)
    kperm = map_dfr(perms, kf)
    kperm = kperm %>%
      group_by(r) %>%
      summarise(csr = mean(correction_k),
                var = var(correction_k),
                pvalue = sum(correction_k >= khat)/nperm
      ) %>% # empirical CSR
      ungroup() %>%
      mutate(correction_k = k$correction_k,
             fundiff = correction_k - csr,
             method = "perm",
             initial_correction = k$initial_correction,
             theo_csr = k$theo_csr) %>%
      select(r, theo_csr, kperm_csr = csr, correction_k, fundiff, method, initial_correction)

    res = kperm %>%
      mutate(khat = rep(k$correction_k, times = 1))


    res2 = res %>%
      mutate(Z = (khat - kperm_csr) / sqrt(var), # kperm_csr = kperm expectation/new csr
             pvalue = pnorm(-Z),
             method = "kperm approx")

    kperm = bind_rows(res, res2) %>%
      mutate(fundiff = khat - kperm_csr) %>%
      select(r, theo_csr, kperm_csr, khat, fundiff, var, Z, pvalue, method, initial_correction)
  }

  return(kperm)
}

