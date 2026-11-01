#' B-KAMP for scalable clustering inference
#' 
#' @title B-KAMP
#' 
#' @description Identify adaptive blocks and combine their standardized KAMP
#'   statistics using a weighting scheme. Larger statistics indicate
#'   clustering or colocalization relative to random labeling.
#'   
#' @note Blocks are identified once and reused at all radii. Undefined block
#'   statistics give combined NA results; radius 0 is typically undefined.
#'   
#' @param ppp_obj A marked `ppp` with a rectangular observation window.
#' @param univariate A logical value indicating whether to compute univariate KAMP (default is TRUE).
#' @param mark1 Variable used to mark the points in the point pattern object for the first type.
#' @param mark2 Variable used to mark the points in the point pattern object for the second type (optional, only used if `univariate` is FALSE).
#' @param split_max Maximum point count per initial grid cell. NULL uses `sqrt(n)`, where n denotes the total number of points in `ppp_obj`.
#' @param n_min Minimum target/background counts, or type1/type2/background counts.
#'   NULL uses the paper's square-root thresholds based on full-pattern counts.
#' @param rho1 Phase I bound on longer/shorter side ratio. NULL uses the full-window ratio.
#' @param rho2 Phase II aspect ratio bound. Defaults to Inf; must be >= rho1.
#' @param max_grid_dim Maximum number of grid cells per axis to search. Defaults to 100.
#' @param rvals Relative radii. Each value is multiplied by each block's shorter
#'   side to obtain the absolute radius used for within-block K calculations.
#' @param correction Type of edge correction. Defaults to translational.
#'   
#' @returns A tibble with the following columns:
#' \describe{
#'   \item{r}{The relative radius at which K was calculated within each block.}
#'   \item{z}{Weighted statistic combining block-level KAMP z-statistics.}
#'   \item{pvalue}{right-sided normal p-value, `pnorm(-z)`.}
#' }
#' 
#' @seealso [bkamp_block()], [bkamp_Rcpp()]
#' 
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @importFrom stats pnorm
#' 
#' @export
#' @examples
#' # one image from the ovarian cancer data in the package
#' ids <- unique(ovarian_df$sample_id)
#' img <- ovarian_df[ovarian_df$sample_id == ids[4], ]
#' pp <- spatstat.geom::ppp(img$x, img$y,
#'                          window = spatstat.geom::owin(range(img$x), range(img$y)),
#'                          marks = img$immune)
#'
#' bkamp(pp, rvals = c(0.05, 0.1, 0.15), mark1 = "immune")
bkamp <- function(ppp_obj, rvals = c(0, .05, .075, .1, .15, .2),
                    univariate = TRUE, correction = "trans",
                    mark1,
                    mark2 = NULL,
                    split_max = NULL, n_min = NULL, rho1 = NULL, rho2 = Inf,
                    max_grid_dim = 100L) {
  correction <- check_bkamp_inference(rvals, univariate, correction)
  
  blocks <- bkamp_block(ppp_obj = ppp_obj, univariate = univariate,
                         mark1 = mark1, mark2 = mark2, split_max = split_max,
                         n_min = n_min, rho1 = rho1, rho2 = rho2,
                         max_grid_dim = max_grid_dim)
  
  purrr::map_dfr(rvals,
    ~bkamp_helper(blocks = blocks, rvalue = .x, univariate = univariate,
                   correction = correction, mark1 = mark1, mark2 = mark2),
    .progress = TRUE)
}
