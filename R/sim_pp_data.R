#' @title Simulates spatial data
#' @description Simulates a spatial point pattern with two cell types: background and immune.
#'
#' @param lambda_n Number of total cells in image
#' @param abundance Percentage intensity for marker positive cells
#' @param clust Determines whether an image is simulated with or without clustering (TRUE/FALSE)
#' @param cell_type1 Marker positive cell type (default is "immune")
#' @param cell_type2 Marker negative cell type (default is "background")
#' @param distribution Determines whether the image is homogeneous ("hom") or inhomogeneous ("inhom")
#'
#' @returns A point pattern object of class "ppp" from the spatstat package where the two cell types are background and immune.
#'
#' @importFrom spatstat.explore Kcross Kest
#' @importFrom spatstat.geom area.owin ppp as.owin
#' @importFrom dplyr mutate select rename filter
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_point
#' @importFrom spatstat.random rmpoispp
#' @importFrom scSpatialSIM CreateSimulationObject GenerateSpatialPattern GenerateCellPositivity GenerateHoles CreateSpatialList
#'
#' @export
#'
#' @examples
#' if (requireNamespace("dplyr", quietly = TRUE) &&
#'    requireNamespace("ggplot2", quietly = TRUE) &&
#'    requireNamespace("tibble", quietly = TRUE) &&
#'    requireNamespace("magrittr", quietly = TRUE)) {
#'   pp_obj <- sim_pp_data(lambda_n = 200, abundance = 0.3)
#'   magrittr::`%>%`(
#'   tibble::as_tibble(pp_obj),
#'   ggplot2::ggplot(ggplot2::aes(x, y, color = marks))
#'   ) +
#'   ggplot2::geom_point()
#' }
#'
sim_pp_data <- function(lambda_n,
                        abundance, # needs to be divisible by 5
                        cell_type1 = "immune",
                        cell_type2 = "background",
                        distribution = "hom",
                        clust = FALSE){

  wm <- spatstat.geom::owin(xrange = c(0, 1), yrange = c(0, 1))
  pp_obj = NULL

  if (clust == FALSE) {
    #lambda_immune = round((lambda_n * abundance)/(1-abundance))
    #lambda_background = lambda_n - lambda_immune

    lambda_immune <- round(lambda_n * abundance)
    lambda_background <- lambda_n * (1 - abundance)


    if(distribution == "hom"){
      # homogeneous background and immune
      pp_obj = rmpoispp(c(lambda_immune, lambda_background), types = c(cell_type1, cell_type2),
                        win = wm)
    }else if(distribution == "inhom"){
      # inhomogeneous background and inhomogeneous immune
      lams <- list(function(x,y){lambda_immune*5*x^2},
                   function(x,y){lambda_background*5*x^2}
                   )

      pp_obj = rmpoispp(lams, types = c(cell_type1, cell_type2),
                        win = wm)
    }

  } else {
    # larger window for clustering
    wm <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
    sim_object = CreateSimulationObject(sims = 1, cell_types = 1, window = wm)
    sim_object = GenerateSpatialPattern(sim_object,
                                        lambda = lambda_n/100)

    sim_object = GenerateCellPositivity(sim_object,
                                        k = 25,
                                        sdmin = .7, sdmax = .71,
                                        step_size = 0.1, cores = 1,
                                        probs = c(0.0, 1))

    # Generate holes for inhomogeneity
    if(distribution == "inhom"){
      sim_object = GenerateHoles(sim_object, step_size = 0.1, cores = 1)
    }

    pp = CreateSpatialList(sim_object, single_df = TRUE) %>%
      rename(immune = `Cell 1 Assignment`)

    phat = sum(pp$immune)/nrow(pp)

    if(phat > abundance){
      nhat = nrow(pp)
      nthin = round((phat - abundance) * nhat)

      indices = which(pp$immune == 1)
      pp$immune[sample(indices, nthin)] <- 0
    }

    if(distribution == "inhom"){
      pp = pp %>%
        rename(hole = `Hole Assignment`) %>%
        filter(hole == "Keep") %>%
        select(-hole)

    }

    pp = pp %>%
      mutate(immune = ifelse(immune == 0, cell_type2, cell_type1))

    pp_obj = spatstat.geom::ppp(pp$x, pp$y, window = wm,  marks = factor(pp$immune))

  }

  return(pp_obj)
}

#' Simulates bivariate point process data
#'
#' @title Simulates bivariate point process data
#' @description Simulates a spatial point pattern with three cell types: background, immune1, and immune2.
#'
#' @param lambda_n Number of total cells in image
#' @param abundance Percentage intensity for marker positive cells
#' @param cell_type1 Marker positive cell type (default is "immune1")
#' @param cell_type2 Marker positive cell type (default is "immune2")
#' @param cell_type3 Marker negative cell type (default is "background")
#' @param distribution Determines whether the image is homogeneous ("hom") or inhomogeneous ("inhom")
#' @param clust Determines whether an image is simulated with or without clustering (TRUE/FALSE)
#'
#' @returns
#'
#' #' @importFrom spatstat.explore Kcross Kest
#' @importFrom spatstat.geom area.owin ppp as.owin
#' @importFrom dplyr mutate select rename filter
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_point
#' @importFrom spatstat.random rmpoispp
#' @importFrom scSpatialSIM CreateSimulationObject GenerateSpatialPattern GenerateCellPositivity GenerateHoles CreateSpatialList
#'
#' @export
#'
#' @examples
#' if (requireNamespace("dplyr", quietly = TRUE) &&
#'    requireNamespace("ggplot2", quietly = TRUE) &&
#'    requireNamespace("tibble", quietly = TRUE) &&
#'    requireNamespace("magrittr", quietly = TRUE)) {
#'   pp_obj <- sim_pp_data_bivariate(lambda_n = 200, abundance = 0.3)
#'   magrittr::`%>%`(
#'   tibble::as_tibble(pp_obj),
#'   ggplot2::ggplot(ggplot2::aes(x, y, color = marks))
#'   ) +
#'   ggplot2::geom_point()
#' }
#'
sim_pp_data_bivariate <- function(lambda_n,
                                  abundance,
                                  cell_type1 = "immune1",
                                  cell_type2 = "immune2",
                                  cell_type3 = "background",
                                  distribution = "hom",
                                  clust = FALSE){

  wm <- spatstat.geom::owin(xrange = c(0, 1), yrange = c(0, 1))
  pp_obj = NULL

  lambda_immune <- round(lambda_n * abundance)
  lambda_background <- lambda_n * (1 - abundance)

  if (clust == FALSE) {
    if(distribution == c("inhom")){
      lams <- list(function(x,y){lambda_immune*5*x^2},
                   function(x,y){lambda_immune*5*x^2},
                   function(x,y){lambda_background*5*x^2}
      )

      pp_obj = rmpoispp(lams,
                        types = c("immune1", "immune2", "background"),
                        win = wm)
    } else {
      pp_obj = rmpoispp(c(lambda_immune, lambda_immune, lambda_background),
                        types = c("immune1", "immune2", "background"),
                        win = wm)
    }


  } else {
    wm <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
    sim_object = CreateSimulationObject(sims = 1, cell_types = 2, window = wm)
    sim_object = GenerateSpatialPattern(sim_object,
                                        lambda = lambda_n/100)

    sim_object = GenerateCellPositivity(sim_object,
                                        k = 25,
                                        sdmin = .7, sdmax = .71,
                                        step_size = 0.1, cores = 1,
                                        probs = c(0.0, 1, 1))

    if (distribution == "inhom") {
      sim_object = GenerateHoles(sim_object, step_size = 0.1, cores = 1)
    }

    pp = CreateSpatialList(sim_object, single_df = TRUE) %>%
      rename(immune1 = `Cell 1 Assignment`,
             immune2 = `Cell 2 Assignment`)

    phat1 = sum(pp$immune1)/nrow(pp)
    phat2 = sum(pp$immune2)/nrow(pp)

    if (phat1 > abundance){
      nhat = nrow(pp)
      nthin = round((phat1 - abundance) * nhat)
      indices = which(pp$immune1 == 1)
      pp$immune1[sample(indices, nthin)] <- 0
    }

    if (phat2 > abundance){
      nhat = nrow(pp)
      nthin = round((phat2 - abundance) * nhat)
      indices = which(pp$immune2 == 1)
      pp$immune2[sample(indices, nthin)] <- 0
    }

    if (distribution == "inhom"){
      pp = pp %>%
        rename(hole = `Hole Assignment`) %>%
        filter(hole == "Keep") %>%
        select(-hole)
    }

    pp = pp %>%
      mutate(immune = case_when(immune1 == 1 ~ "immune1",
                                immune2 == 1 ~ "immune2",
                                TRUE ~ "background"))
    pp_obj = spatstat.geom::ppp(pp$x, pp$y, window = wm,  marks = factor(pp$immune))
  }

  return(pp_obj)
}
