#' Title
#'
#' @param lambda_n intensity for background cells
#' @param abundance intensity for marker positive cells
#' @param type defines the distribution of the point process- homogeneous, inhomogeneous, or clustered
#' @param clust should an image be simulated with or without holes
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
get_sim_data <- function(lambda_n,
                         abundance, # needs to be divisible by 5
                         type = c("hom", "inhom"),
                         clust = FALSE){

  wm <- spatstat.geom::owin(xrange = c(0, 1), yrange = c(0, 1))
  pp_obj = NULL

  if (clust == FALSE) {
    #lambda_immune = round((lambda_n * abundance)/(1-abundance))
    #lambda_background = lambda_n - lambda_immune

    lambda_immune <- round(lambda_n * abundance)
    lambda_background <- lambda_n * (1 - abundance)


    if(type %in% c("inhom")){
      lams <- list(function(x,y){lambda_immune*5*x^2},
                   function(x,y){lambda_background*5*x^2}
      )
    }

    if(type == "hom"){
      # homogeneous background and immune
      pp_obj = rmpoispp(c(lambda_immune, lambda_background), types = c("immune", "background"),
                        win = wm)
    }else if(type == "inhom"){
      # inhomogeneous background and inhomogeneous immune
      pp_obj = rmpoispp(lams, types = c("immune", "background"),
                        win = wm)
    }

    as_tibble(pp_obj) %>%
      ggplot(aes(x,y, color = marks)) +
      geom_point()

  } else {
    # larger window
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
    if(type == "inhom"){
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

    if(type == "inhom"){
      pp = pp %>%
        rename(hole = `Hole Assignment`) %>%
        filter(hole == "Keep") %>%
        select(-hole)

    }

    pp = pp %>%
      mutate(immune = ifelse(immune == 0, "background", "immune"))

    pp_obj = spatstat.geom::ppp(pp$x, pp$y, window = wm,  marks = factor(pp$immune))

  }

  return(pp_obj)
}
