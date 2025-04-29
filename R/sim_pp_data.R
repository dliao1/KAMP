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
#' @importFrom scSpatialSIM CreateSimulationObject
#' @importFrom scSpatialSIM GenerateSpatialPattern
#' @importFrom scSpatialSIM GenerateCellPositivity
#' @importFrom scSpatialSIM GenerateHoles
#' @importFrom scSpatialSIM CreateSpatialList
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


    if(distribution == c("inhom")){
      lams <- list(function(x,y){lambda_immune*5*x^2},
                   function(x,y){lambda_background*5*x^2}
      )
    }

    if(distribution == "hom"){
      # homogeneous background and immune
      pp_obj = rmpoispp(c(lambda_immune, lambda_background), types = c(cell_type1, cell_type2),
                        win = wm)
    }else if(distribution == "inhom"){
      # inhomogeneous background and inhomogeneous immune
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
