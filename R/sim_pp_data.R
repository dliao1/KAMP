#' @title Simulates spatial data
#' @description Simulates a spatial point pattern with two cell types: background and immune.
#'
#' @param lambda_n Number of total cells in image
#' @param abundance Percentage intensity for marker positive cells
#' @param clust Determines whether an image is simulated with or without clustering (TRUE/FALSE)
#' @param markvar1 Marker positive cell type (default is "immune")
#' @param markvar2 Marker negative cell type (default is "background")
#' @param distribution Determines whether the image is homogeneous ("hom") or inhomogeneous ("inhom")
#'
#' @returns A point pattern object of class "ppp" from the spatstat package where the two cell types are background and immune.
#'
#' @importFrom spatstat.explore Kcross Kest
#' @importFrom spatstat.geom area.owin ppp as.owin
#' @importFrom dplyr mutate select rename filter case_when
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
                        markvar1 = "immune",
                        markvar2 = "background",
                        distribution = "hom",
                        clust = FALSE){

  if (distribution %in% c("inhom", "hom") == FALSE) {
    stop("distribution must be either 'inhom' or 'hom'")
  }

  if (abundance < 0 || abundance > 1) {
    stop("abundance must be between 0 and 1")
  }

  if (lambda_n <= 0) {
    stop("lambda_n must be greater than 0")
  }

  if (markvar1 == markvar2) {
    stop("markvar1 and markvar2 must be different")
  }



  wm <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  A  <- spatstat.geom::area.owin(wm)

  pp_obj = NULL

  if (clust == FALSE) {
    #lambda_immune = round((lambda_n * abundance)/(1-abundance))
    #lambda_background = lambda_n - lambda_immune

    lambda_immune <- (lambda_n * abundance) / A
    lambda_background <- (lambda_n * (1 - abundance)) / A


    if(distribution == "hom"){
      # homogeneous background and immune
      pp_obj = rmpoispp(c(lambda_immune, lambda_background), types = c(markvar1, markvar2),
                        win = wm)
    }else if(distribution == "inhom"){
      # inhomogeneous background and inhomogeneous immune
      lams <- list(
        function(x, y) { lambda_immune     * 3 * (x/10)^2 },
        function(x, y) { lambda_background * 3 * (x/10)^2 }
      )

      pp_obj = rmpoispp(lams, types = c(markvar1, markvar2),
                        win = wm)
    }

  } else {
    N_target <- lambda_n

    keep_target <- 0.4

    sim_object = CreateSimulationObject(sims = 1, cell_types = 1, window = wm)
    sim_object <- GenerateSpatialPattern(
      sim_object,
      lambda = (N_target / keep_target) / A
    )

    sim_object <- GenerateCellPositivity(
      sim_object,
      k = 25,
      sdmin = 0.7, sdmax = 0.71,
      step_size = 0.1, cores = 1,
      probs = c(0, 1)
    )

    # Generate holes for inhomogeneity
    if(distribution == "inhom"){
      sim_object <- GenerateHoles(
        sim_object,
        use_window = TRUE,
        density_heatmap = TRUE,
        hole_prob = c(0.9, 0.995),
        sdmin = 1.25,
        sdmax = 2.25,
        step_size = 0.05,
        overwrite = TRUE
      )
    }

    pp = CreateSpatialList(sim_object, single_df = TRUE) %>%
      rename(immune = `Cell 1 Assignment`)

    # apply holes first for inhom only
    if (distribution == "inhom") {
      pp <- pp %>%
        dplyr::rename(hole = `Hole Assignment`) %>%
        dplyr::filter(hole == "Keep") %>%
        dplyr::select(-hole)
    }

    # enforce exactly N_target points (downsample if needed)
    if (nrow(pp) > N_target) {
      pp <- pp[sample.int(nrow(pp), N_target), ]
    }

    phat = sum(pp$immune)/nrow(pp)

    if(phat > abundance){
      nhat = nrow(pp)
      nthin = round((phat - abundance) * nhat)

      indices = which(pp$immune == 1)
      pp$immune[sample(indices, nthin)] <- 0
    }


    pp = pp %>%
      mutate(immune = ifelse(immune == 0, markvar2, markvar1))

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
#' @param markvar1 Marker positive cell type (default is "immune1")
#' @param markvar2 Marker positive cell type (default is "immune2")
#' @param markvar3 Marker negative cell type (default is "background")
#' @param distribution Determines whether the image is homogeneous ("hom") or inhomogeneous ("inhom")
#' @param clust Determines whether an image is simulated with or without clustering (TRUE/FALSE)
#'
#' @returns A point pattern object of class "ppp" from the spatstat package where the three cell types are background, immune1, and immune2.
#'
#' @importFrom spatstat.explore Kcross Kest
#' @importFrom spatstat.geom area.owin ppp as.owin
#' @importFrom dplyr mutate select rename filter case_when
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
#'   pp_obj <- sim_pp_data_biv(lambda_n = 200, abundance = 0.3)
#'   magrittr::`%>%`(
#'   tibble::as_tibble(pp_obj),
#'   ggplot2::ggplot(ggplot2::aes(x, y, color = marks))
#'   ) +
#'   ggplot2::geom_point()
#' }
#'
sim_pp_data_biv <- function(lambda_n,
                            abundance,
                            markvar1 = "immune1",
                            markvar2 = "immune2",
                            markvar3 = "background",
                            distribution = "hom",
                            clust = FALSE){

  if (distribution %in% c("inhom", "hom") == FALSE) {
    stop("distribution must be either 'inhom' or 'hom'")
  }

  if (markvar1 == markvar2) {
    stop("markvar1 and markvar2 must be different")
  }

  if (markvar1 == markvar3 || markvar2 == markvar3) {
    stop("markvar1 and markvar3, or markvar2 and markvar3 must be different")
  }

  if (abundance < 0 || abundance > 1) {
    stop("abundance must be between 0 and 1")
  }

  if (lambda_n <= 0) {
    stop("lambda_n must be greater than 0")
  }


  wm <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  A  <- spatstat.geom::area.owin(wm)

  pp_obj = NULL

  lambda_immune <- (lambda_n * abundance / 2) / A
  lambda_background  <- (lambda_n * (1 - abundance)) / A

  if (clust == FALSE) {
    if(distribution == "inhom"){
      lams <- list(
        function(x,y){ lambda_immune     * 3 * (x/10)^2 },
        function(x,y){ lambda_immune     * 3 * (x/10)^2 },
        function(x,y){ lambda_background * 3 * (x/10)^2 }
      )

      pp_obj = rmpoispp(lams,
                        types = c(markvar1, markvar2, markvar3),
                        win = wm)
    } else {
      pp_obj = rmpoispp(c(lambda_immune, lambda_immune, lambda_background),
                        types = c(markvar1, markvar2, markvar3),
                        win = wm)
    }


  } else {
    N_target <- lambda_n
    keep_target <- 0.6

    sim_object = CreateSimulationObject(sims = 1, cell_types = 2, window = wm)
    sim_object = GenerateSpatialPattern(sim_object,
                                        lambda = (N_target / keep_target) / A)

    sim_object = GenerateCellPositivity(sim_object,
                                        k = 25,
                                        sdmin = .7, sdmax = .71,
                                        step_size = 0.1, cores = 1,
                                        probs = c(0.0, 1, 1))

    if (distribution == "inhom") {
      sim_object <- GenerateHoles(
        sim_object,
        use_window = TRUE,
        density_heatmap = TRUE,
        hole_prob = c(0.9, 0.995),
        sdmin = 1.25, # bigger holes
        sdmax = 2.25,
        step_size = 0.05,
        overwrite = TRUE
      )
    }

    pp = CreateSpatialList(sim_object, single_df = TRUE) %>%
      rename(immune1 = `Cell 1 Assignment`,
             immune2 = `Cell 2 Assignment`)



    if (distribution == "inhom"){
      pp = pp %>%
        rename(hole = `Hole Assignment`) %>%
        filter(hole == "Keep") %>%
        select(-hole)
    }

    if (nrow(pp) > N_target) {
      pp <- pp[sample.int(nrow(pp), N_target), ]
    }


    phat1 = sum(pp$immune1)/nrow(pp)
    phat2 = sum(pp$immune2)/nrow(pp)

    target <- abundance/2
    if (phat1 > target){
      nhat <- nrow(pp)
      nthin <- round((phat1 - target) * nhat)
      idx <- which(pp$immune1 == 1)
      pp$immune1[sample(idx, nthin)] <- 0
    }

    if (phat2 > target){
      nhat <- nrow(pp)
      nthin <- round((phat2 - target) * nhat)
      idx <- which(pp$immune2 == 1)
      pp$immune2[sample(idx, nthin)] <- 0
    }


    pp = pp %>%
      mutate(immune = case_when(immune1 == 1 ~ markvar1,
                                immune2 == 1 ~ markvar2,
                                TRUE ~ markvar3))
    pp_obj = spatstat.geom::ppp(pp$x, pp$y, window = wm,  marks = factor(pp$immune))
  }

  return(pp_obj)
}
