#' Adaptive spatial blocking
#' 
#' @title Adaptive spatial blocking for B-KAMP
#' 
#' @description Partition a marked pattern into rectangular blocks using grid
#'   representation, two block identification phases, and block complement.
#'   
#' @note The observation window must be rectangular. Returned blocks retain their
#'   original marks. Unassigned points are excluded with a warning.
#'   
#' @param ppp_obj A marked `ppp` with a rectangular observation window.
#' @param univariate Logical; TRUE for one target type, FALSE for two.
#' @param mark1 Variable used to mark the points in the point pattern object for the first type.
#' @param mark2 Variable used to mark the points in the point pattern object for the second type (optional, only used if `univariate` is FALSE).
#' @param split_max Maximum point count per initial grid cell. NULL uses `sqrt(n)`, where n denotes the total number of points in `ppp_obj`.
#' @param n_min Minimum target/background counts, or type1/type2/background counts.
#'   NULL uses the paper's square-root thresholds based on full-pattern counts.
#' @param rho1 Phase I bound on longer/shorter side ratio. NULL uses the full-window ratio.
#' @param rho2 Phase II aspect ratio bound. Defaults to Inf; must be >= rho1.
#' @param max_grid_dim Maximum number of grid cells per axis to search. Defaults to 100.
#' 
#' @returns A list whose elements are marked `ppp` blocks.
#' 
#' @importFrom spatstat.geom Window npoints
#' @importFrom tibble tibble
#' @importFrom stats pnorm
#' 
#' @export
#' @examples
#' set.seed(100)
#' win <- spatstat.geom::owin(c(0, 1), c(0, 1))
#' pp <- spatstat.geom::ppp(runif(500), runif(500), window = win,
#'   marks = factor(sample(c("immune", "background"), 500, replace = TRUE)))
#' blocks <- bkamp_block(pp, mark1 = "immune")
#' length(blocks)
bkamp_block <- function(ppp_obj, univariate = TRUE,
                         mark1,
                         mark2 = NULL,
                         split_max = NULL, n_min = NULL, rho1 = NULL, rho2 = Inf,
                         max_grid_dim = 100L) {
  check_bkamp_block_options(univariate, split_max, n_min, rho1, rho2,
                            max_grid_dim)
  check_bkamp_input(ppp_obj, univariate, mark1, mark2)
  setup <- bkamp_setup(ppp_obj, univariate, mark1, mark2,
                       split_max, n_min, rho1, rho2)

  # Grid representation: at most split_max points per cell (sqrt(n) by default).
  grid <- bkamp_grid(ppp_obj, setup$split_max, max_grid_dim)
  grid <- lapply(grid, bkamp_marks, univariate = univariate,
                 target_mark = mark1, type1_mark = mark1, type2_mark = mark2)
  ref_window <- spatstat.geom::Window(grid[[1]])
  ref_wd <- diff(ref_window$xrange)
  ref_hi <- diff(ref_window$yrange)
  base_tol <- sqrt(.Machine$double.eps) * 10

  # Identify initially valid cells
  assigned <- vapply(grid, bkamp_counts, logical(1),
                     n_min = setup$n_min, univariate = univariate)
  
  if (max(ref_wd, ref_hi) / min(ref_wd, ref_hi) > setup$rho1) assigned[] <- FALSE
  
  success_blocks <- lapply(which(assigned), function(i) list(ppp = grid[[i]], ids = i))

  # Phase I: repeatedly identify blocks under the initial aspect ratio bound.
  phase1 <- bkamp_phase(setup$rho1, grid, setup$n_min, success_blocks, assigned,
                        ref_wd, ref_hi, base_tol, univariate = univariate)

  # Phase II: identify further blocks under the relaxed aspect ratio bound.
  phase2 <- bkamp_phase(setup$rho2, grid, setup$n_min,
                        phase1$success_blocks, phase1$assigned,
                        ref_wd, ref_hi, base_tol, univariate = univariate)

  # Phase III: expand retained blocks to incorporate adjacent residual cells.
  partition <- bkamp_complement(grid, phase2$success_blocks, phase2$assigned, base_tol)
  if (length(partition$success_blocks) == 0L) {
    stop("No blocks satisfy the count and shape constraints.", call. = FALSE)
  }

  # Return ppp blocks with their original marks and unique point ownership.
  blocks <- lapply(partition$success_blocks, function(x) {
    bkamp_subset(ppp_obj, spatstat.geom::Window(x$ppp),
                  ppp_obj$window$xrange[2], ppp_obj$window$yrange[2])
  })

  for (i in seq_along(blocks)) {
    check_bkamp_pattern(blocks[[i]], univariate, mark1, mark2,
                        label = paste0("Block ", i))
  }
  
  n_residual <- ppp_obj$n - sum(vapply(blocks, spatstat.geom::npoints, numeric(1)))
  if (n_residual > 0L) {
    warning(sprintf("%d points remain in residual cells and are excluded from B-KAMP.",
                    n_residual), call. = FALSE)
  }
  blocks
}
