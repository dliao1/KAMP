#' Construct the adaptive regular grid satisfying the maximum cell count
#' 
#' @description Refine a regular grid over the rectangular observation window
#'   until every cell satisfies `split_max`.
#'   
#' @param ppp_obj A marked `ppp` with a rectangular observation window.
#' @param split_max Maximum point count per cell, controlling grid refinement.
#' @param max_grid_dim Maximum number of cells per axis to search. Defaults to 100.
#' 
#' @returns A list of marked `ppp` grid cells with distinct point ownership.
#' @keywords internal
bkamp_grid <- function(ppp_obj, split_max, max_grid_dim = 100L) {
  win <- spatstat.geom::Window(ppp_obj)
  
  for (N in seq.int(2L, max_grid_dim)) {
    xb <- seq(win$xrange[1], win$xrange[2], length.out = N + 1L)
    yb <- seq(win$yrange[1], win$yrange[2], length.out = N + 1L)
    
    # Each point belongs to one grid cell, including points on internal grid lines.
    col <- findInterval(ppp_obj$x, xb, rightmost.closed = TRUE)
    row <- findInterval(ppp_obj$y, yb, rightmost.closed = TRUE)
    ids <- (row - 1L) * N + col
    
    if (max(tabulate(ids, nbins = N * N)) > split_max) next
    
    grid <- vector("list", N * N)
    for (i in seq_len(N)) {
      for (j in seq_len(N)) {
        
        id <- (i - 1L) * N + j
        cell_win <- spatstat.geom::owin(xb[c(j, j + 1L)], yb[c(i, i + 1L)])
        keep <- ids == id
        grid[[id]] <- spatstat.geom::ppp(ppp_obj$x[keep], ppp_obj$y[keep],
          window = cell_win, marks = spatstat.geom::marks(ppp_obj)[keep], checkdup = FALSE)
      }
    }
    return(grid)
  }
  stop(sprintf("No grid up to %dx%d satisfies split_max; increase max_grid_dim.",
               max_grid_dim, max_grid_dim), call. = FALSE)
}

#' Logical index of points in `pp` inside `win`, assigning shared edges to the
#' rightmost and topmost cell only.
#' @noRd
bkamp_inside <- function(pp, win, xmax, ymax) {
  pp$x >= win$xrange[1] & pp$y >= win$yrange[1] &
    (pp$x < win$xrange[2] | (win$xrange[2] == xmax & pp$x <= xmax)) &
    (pp$y < win$yrange[2] | (win$yrange[2] == ymax & pp$y <= ymax))
}

#' Subset `pp` to points inside `win` using `bkamp_inside` ownership rules,
#' returning a `ppp` with the original marks.
#' @noRd
bkamp_subset <- function(pp, win, xmax, ymax) {
  keep <- bkamp_inside(pp, win, xmax, ymax)
  spatstat.geom::ppp(pp$x[keep], pp$y[keep], window = win,
                   marks = spatstat.geom::marks(pp)[keep], checkdup = FALSE)
}

#' Remap user-supplied marks to internal labels (`"target"`, `"background"`,
#' or `"target1"`, `"target2"`, `"background"`) for consistent blocking logic.
#' @noRd
bkamp_marks <- function(pp, univariate = TRUE,
                         target_mark = "target", type1_mark = "target1",
                         type2_mark = "target2") {
  mk <- as.character(spatstat.geom::marks(pp))
  
  # Compare against the original labels simultaneously to avoid name collisions.
  mapped <- if (!univariate) {
    ifelse(mk == type1_mark, "target1",
           ifelse(mk == type2_mark, "target2", "background"))
  } else {
    ifelse(mk == target_mark, "target", "background")
  }
  spatstat.geom::marks(pp) <- factor(mapped,
    levels = if (!univariate) c("target1", "target2", "background") else c("target", "background"))
  pp
}

#' Return TRUE if `pp` satisfies the minimum target and background count
#' constraints in `n_min` (univariate: target, background; bivariate: type1, type2, background).
#' @noRd
bkamp_counts <- function(pp, n_min, univariate = TRUE) {
  if (!univariate) {
    target1 <- sum(pp$marks == "target1")
    target2 <- sum(pp$marks == "target2")
    bg <- sum(pp$marks == "background")
    
    return((target1 >= n_min[1]) &
             (target2 >= n_min[2]) &
             (bg >= n_min[3]))
  }
  
  target <- sum(pp$marks == "target")
  bg <- sum(pp$marks == "background")
  return((target >= n_min[1]) & (bg >= n_min[2]))
  
}

#' Return TRUE if the area of `win_target` not covered by `win_actual` is
#' within `tol`, i.e. the two windows are the same rectangle.
#' @noRd
bkamp_window <- function(win_target, win_actual, tol) {
  unfilled_area <- spatstat.geom::area.owin(
    spatstat.geom::setminus.owin(win_target, win_actual)
  )
  unfilled_area <= tol
}

#' Return the directions ("left", "right", "up", "down") in which `win_b`
#' touches `win_a`, requiring shared-edge overlap greater than `tol`.
#' @noRd
bkamp_adjacent <- function(win_a, win_b, tol) {
  dirs <- c()
  
  overlap_y <- (min(win_a$yrange[2], win_b$yrange[2]) -
                  max(win_a$yrange[1], win_b$yrange[1])) > tol
  overlap_x <- (min(win_a$xrange[2], win_b$xrange[2]) -
                  max(win_a$xrange[1], win_b$xrange[1])) > tol
  
  if (abs(win_a$xrange[1] - win_b$xrange[2]) < tol && overlap_y) {
    dirs <- c(dirs, "left")
  }
  if (abs(win_a$xrange[2] - win_b$xrange[1]) < tol && overlap_y) {
    dirs <- c(dirs, "right")
  }
  if (abs(win_a$yrange[1] - win_b$yrange[2]) < tol && overlap_x) {
    dirs <- c(dirs, "down")
  }
  if (abs(win_a$yrange[2] - win_b$yrange[1]) < tol && overlap_x) {
    dirs <- c(dirs, "up")
  }
  
  unique(dirs)
}

#' Compute default `n_min` (sqrt-thresholds), `split_max`, `rho1`, and `rho2`
#' from full-pattern counts when the user supplies NULL for any of them.
#' @noRd
bkamp_setup <- function(pp, univariate, mark1, mark2, split_max, n_min, rho1, rho2) {
  mk <- spatstat.geom::marks(pp)
  m1 <- sum(mk == mark1)
  m2 <- if (univariate) 0L else sum(mk == mark2)
  counts <- if (univariate) c(m1, pp$n - m1) else c(m1, m2, pp$n - m1 - m2)
  if (any(counts < 1L)) {
    stop("B-KAMP requires target and background points in the full pattern.", call. = FALSE)
  }
  if (is.null(n_min)) {
    roots <- sqrt(counts)
    n_min <- c(roots[-length(roots)], min(roots))
  }
  n_min <- ceiling(n_min)
  if (any(n_min > counts)) {
    stop("n_min exceeds the corresponding full-pattern counts.", call. = FALSE)
  }
  if (is.null(split_max)) split_max <- sqrt(pp$n)
  side <- c(diff(pp$window$xrange), diff(pp$window$yrange))
  if (is.null(rho1)) rho1 <- max(side) / min(side)
  if (rho2 < rho1) stop("rho2 must be at least rho1.", call. = FALSE)
  list(split_max = split_max, n_min = n_min, rho1 = rho1, rho2 = rho2)
}

#' Validate that `mark1` (and `mark2` for bivariate) are present in `ppp_obj`,
#' stopping with an informative message on failure.
#' @noRd
check_bkamp_marks <- function(ppp_obj, univariate, mark1, mark2 = NULL,
                              label = NULL) {
  prefix <- if (is.null(label)) "" else paste0(label, ": ")
  fail <- function(msg) stop(paste0(prefix, msg), call. = FALSE)
  if (!inherits(ppp_obj, "ppp")) fail("Input must be a ppp object.")
  if (!is.logical(univariate) || length(univariate) != 1L || is.na(univariate)) {
    fail("univariate must be TRUE or FALSE.")
  }
  mk <- spatstat.geom::marks(ppp_obj)
  if (is.null(mk) || length(mk) == 0L) {
    fail("The point pattern object does not have marks.")
  }
  if (!is.atomic(mk) || !is.null(dim(mk)) || anyNA(mk)) {
    fail("Marks must be a vector without missing values.")
  }
  valid_mark <- function(x) !is.null(x) && is.atomic(x) &&
    length(x) == 1L && !is.na(x)
  if (!valid_mark(mark1) || !mark1 %in% mk) {
    fail("mark1 is not a mark in the point pattern object.")
  }
  if (!is.null(mark2) && (!valid_mark(mark2) || !mark2 %in% mk)) {
    fail("mark2 is not a mark in the point pattern object.")
  }
  if (univariate) {
    if (!is.null(mark2)) message(prefix, "mark2 is not used in univariate KAMP. It will be ignored.")
  } else {
    if (is.null(mark2)) fail("Both mark1 and mark2 must be specified for bivariate KAMP.")
    if (mark1 == mark2) fail("mark1 and mark2 cannot be the same for bivariate KAMP.")
  }
  invisible(mk)
}

#' Call `check_bkamp_marks` and optionally warn when target count < 5 or
#' total points exceed 100,000.
#' @noRd
check_bkamp_pattern <- function(ppp_obj, univariate, mark1, mark2 = NULL,
                                report_size = TRUE, label = NULL) {
  mk <- check_bkamp_marks(ppp_obj, univariate, mark1, mark2, label)
  prefix <- if (is.null(label)) "" else paste0(label, ": ")
  n <- spatstat.geom::npoints(ppp_obj)
  m1 <- sum(mk == mark1)
  m2 <- if (univariate) 0L else sum(mk == mark2)
  if (report_size) {
    if (m1 < 5L) {
      message(prefix, "Less than 5 target cells marked as '", mark1,
              "'. This may lead to unreliable results.")
    }
    if (n > 100000L) {
      message(prefix, "Point pattern has more than 100,000 points. At this time, it is not recommended to use edge correction with KAMP on datasets of this size. Results may be unreliable. Consider using KAMP-lite (thin = TRUE) or subsetting your data.")
    }
  }
  invisible(c(n = n, m1 = m1, m2 = m2))
}

#' Identify a valid rectangle using row-wise histograms
#' @description Find a count-valid rectangle with the smallest point count
#'   using the points in the available cells.
#' @param grid_mat Binary availability matrix; rows run from top to bottom.
#' @param max_aspect_ratio Maximum longer/shorter side ratio of a candidate.
#' @param ppp_obj A point pattern object of class "ppp" from the spatstat package.
#' @param n_min Target/background minimum counts, or type1/type2/background counts.
#' @param x_coords,y_coords Increasing grid boundary coordinates.
#' @param ref_wd,ref_hi Width and height of a grid cell.
#' @param univariate A logical value indicating whether to compute univariate KAMP (default is TRUE).
#' @returns The selected rectangle bounds, or NULL if no candidate is valid.
#' @keywords internal
bkamp_rectangle <- function(grid_mat, max_aspect_ratio, ppp_obj, n_min,
                                x_coords, y_coords, ref_wd, ref_hi,
                                univariate = TRUE) {
  if (sum(grid_mat) == 0) return(NULL)
  
  nr <- nrow(grid_mat)
  nc <- ncol(grid_mat)
  heights <- rep(0, nc)
  
  min_total_points <- Inf
  best_rect <- NULL
  
  for (i in seq_len(nr)) {
    heights <- ifelse(grid_mat[i, ] == 1, heights + 1, 0)
    stack <- list()
    
    for (j in seq_len(nc + 1)) {
      h <- if (j <= nc) heights[j] else 0
      
      while (length(stack) > 0 && heights[stack[[length(stack)]]$idx] >= h) {
        popped <- stack[[length(stack)]]
        stack[[length(stack)]] <- NULL
        
        height <- popped$h
        width <- if (length(stack) == 0) j - 1 else j - stack[[length(stack)]]$idx - 1
        
        if (height == 0 || width == 0) next
        
        aspect_ratio <- max(width * ref_wd, height * ref_hi) / min(width * ref_wd, height * ref_hi)
        if (aspect_ratio > max_aspect_ratio) next
        
        start_col <- if (length(stack) == 0) 1 else stack[[length(stack)]]$idx + 1
        x_start <- x_coords[start_col]
        x_end <- x_coords[start_col + width]
        r_start_in_grid <- i - height + 1
        y_start <- y_coords[nr - (r_start_in_grid + height) + 2]
        y_end <- y_coords[nr - r_start_in_grid + 2]
        
        cand_win <- spatstat.geom::owin(xrange = c(x_start, x_end), yrange = c(y_start, y_end))
        points_in_cand <- bkamp_subset(ppp_obj, cand_win,
          xmax = max(x_coords), ymax = max(y_coords))
        
        if (!bkamp_counts(points_in_cand, n_min = n_min, univariate = univariate)) next
        
        if (points_in_cand$n < min_total_points) {
          min_total_points <- points_in_cand$n
          best_rect <- c(row = r_start_in_grid, col = start_col, h = height, w = width)
        }
      }
      
      stack[[length(stack) + 1]] <- list(idx = j, h = h)
    }
  }
  
  return(best_rect)
}

#' Execute an adaptive block identification phase
#' @description Repeatedly extract valid rectangles under one aspect ratio bound.
#' @param max_aspect_ratio Aspect ratio bound for the current phase.
#' @param grid List of grid cells carrying temporary blocking marks.
#' @param n_min Minimum counts required in each candidate block.
#' @param success_blocks Existing block records, each containing `ppp` and cell `ids`.
#' @param assigned Logical vector identifying cells already assigned to blocks.
#' @param ref_wd,ref_hi Width and height of a grid cell.
#' @param base_tol Numerical tolerance for matching grid boundaries.
#' @param univariate A logical value indicating whether to compute univariate KAMP (default is TRUE).
#' @returns A list containing updated `success_blocks` and `assigned`.
#' @keywords internal
bkamp_phase <- function(max_aspect_ratio, grid, n_min,
                                    success_blocks, assigned,
                                    ref_wd, ref_hi, base_tol,
                                    univariate = TRUE) {
  global_win <- do.call(spatstat.geom::boundingbox, lapply(grid, spatstat.geom::Window))
  resolution_col <- min(sapply(grid, function(b) diff(spatstat.geom::Window(b)$xrange)))
  resolution_row <- min(sapply(grid, function(b) diff(spatstat.geom::Window(b)$yrange)))
  
  ncols <- round(diff(global_win$xrange) / resolution_col)
  nrows <- round(diff(global_win$yrange) / resolution_row)
  x_coords <- seq(from = global_win$xrange[1], to = global_win$xrange[2], length.out = ncols + 1)
  y_coords <- seq(from = global_win$yrange[1], to = global_win$yrange[2], length.out = nrows + 1)
  
  grid_mat <- matrix(1, nrow = nrows, ncol = ncols)
  s_wins <- lapply(success_blocks, function(list) spatstat.geom::Window(list$ppp))
  
  for (win in s_wins) {
    c_start <- floor((win$xrange[1] - global_win$xrange[1]) / resolution_col) + 1
    c_end <- floor((win$xrange[2] - global_win$xrange[1] - base_tol) / resolution_col) + 1
    r_start <- floor((win$yrange[1] - global_win$yrange[1]) / resolution_row) + 1
    r_end <- floor((win$yrange[2] - global_win$yrange[1] - base_tol) / resolution_row) + 1
    
    c_start <- max(1, c_start)
    c_end <- min(ncols, c_end)
    r_start <- max(1, r_start)
    r_end <- min(nrows, r_end)
    
    if (c_start <= c_end && r_start <= r_end) {
      grid_mat[(nrows - r_end + 1):(nrows - r_start + 1), c_start:c_end] <- 0
    }
  }
  
  unassigned_ids_current_phase <- which(!assigned)
  if (length(unassigned_ids_current_phase) == 0) {
    return(list(success_blocks = success_blocks, assigned = assigned))
  }
  
  sub_blocks <- grid[unassigned_ids_current_phase]
  x_range <- range(unlist(lapply(sub_blocks, function(x) x$window$xrange)))
  y_range <- range(unlist(lapply(sub_blocks, function(x) x$window$yrange)))
  combined_window <- spatstat.geom::owin(xrange = x_range, yrange = y_range)
  
  unassigned_ppp <- do.call(
    spatstat.geom::superimpose,
    c(sub_blocks, list(W = combined_window))
  )
  
  found_blocks <- list()
  
  repeat {
    candidate_info <- bkamp_rectangle(
      grid_mat = grid_mat,
      max_aspect_ratio = max_aspect_ratio,
      ppp_obj = unassigned_ppp,
      n_min = n_min,
      x_coords = x_coords,
      y_coords = y_coords,
      ref_wd = ref_wd,
      ref_hi = ref_hi,
      univariate = univariate
    )
    
    if (is.null(candidate_info)) break
    
    r_final <- candidate_info[1]
    c_final <- candidate_info[2]
    h_final <- candidate_info[3]
    w_final <- candidate_info[4]
    
    x_start <- x_coords[c_final]
    x_end <- x_coords[c_final + w_final]
    y_start <- y_coords[nrows - (r_final + h_final - 1) + 1]
    y_end <- y_coords[nrows - r_final + 2]
    
    new_found_win <- spatstat.geom::owin(xrange = c(x_start, x_end), yrange = c(y_start, y_end))
    found_blocks[[length(found_blocks) + 1]] <- new_found_win
    
    grid_mat[r_final:(r_final + h_final - 1), c_final:(c_final + w_final - 1)] <- 0
    keep <- !bkamp_inside(unassigned_ppp, new_found_win,
                          xmax = max(x_coords), ymax = max(y_coords))
    unassigned_ppp <- unassigned_ppp[keep]
  }
  
  if (length(found_blocks) > 0) {
    
    for (found_win in found_blocks) {
      contained_ids <- unassigned_ids_current_phase[sapply(grid[unassigned_ids_current_phase], function(b) {
        small_win <- spatstat.geom::Window(b)
        all(small_win$xrange >= found_win$xrange[1] - base_tol &
              small_win$xrange <= found_win$xrange[2] + base_tol) &&
          all(small_win$yrange >= found_win$yrange[1] - base_tol &
                small_win$yrange <= found_win$yrange[2] + base_tol)
      })]
      
      if (length(contained_ids) > 0) {
        sub_blocks <- grid[contained_ids]
        x_range <- range(unlist(lapply(sub_blocks, function(x) x$window$xrange)))
        y_range <- range(unlist(lapply(sub_blocks, function(x) x$window$yrange)))
        combined_window <- spatstat.geom::owin(xrange = x_range, yrange = y_range)
        
        new_ppp <- do.call(
          spatstat.geom::superimpose,
          c(sub_blocks, list(W = combined_window))
        )
        
        success_blocks[[length(success_blocks) + 1]] <- list(ppp = new_ppp, ids = contained_ids)
        assigned[contained_ids] <- TRUE
      }
    }
  }
  
  return(list(success_blocks = success_blocks, assigned = assigned))
}

#' Complement retained blocks with adjacent residual cells
#' @description Greedily expand blocks while keeping them rectangular and disjoint.
#' @param grid List of grid cells carrying temporary blocking marks.
#' @param success_blocks Retained block records, each containing `ppp` and cell `ids`.
#' @param assigned Logical vector identifying cells already assigned to blocks.
#' @param base_tol Numerical tolerance for adjacency and rectangle checks.
#' @returns A list containing updated `success_blocks` and `assigned`.
#' @keywords internal
bkamp_complement <- function(grid, success_blocks, assigned, base_tol) {
  repeat {
    merge_made_in_this_pass <- FALSE
    leftover_ids <- which(!assigned)
    
    if (length(leftover_ids) == 0) break
    
    for (orphan_idx in leftover_ids) {
      local_candidate_merges <- list()
      
      for (success_idx in seq_along(success_blocks)) {
        orphan_win <- spatstat.geom::Window(grid[[orphan_idx]])
        success_win <- spatstat.geom::Window(success_blocks[[success_idx]]$ppp)
        
        if (length(bkamp_adjacent(orphan_win, success_win, tol = base_tol)) == 0) {
          next
        }
        
        cand_bbox <- spatstat.geom::boundingbox(orphan_win, success_win)
        
        ids_needed_to_fill <- which(sapply(grid, function(b) {
          w <- spatstat.geom::Window(b)
          (w$xrange[1] >= cand_bbox$xrange[1] - base_tol) &&
            (w$xrange[2] <= cand_bbox$xrange[2] + base_tol) &&
            (w$yrange[1] >= cand_bbox$yrange[1] - base_tol) &&
            (w$yrange[2] <= cand_bbox$yrange[2] + base_tol)
        }))
        
        current_success_ids <- success_blocks[[success_idx]]$ids
        valid_ids_for_merge <- c(current_success_ids, leftover_ids)
        
        if (!all(ids_needed_to_fill %in% valid_ids_for_merge)) next
        
        sub_blocks <- grid[ids_needed_to_fill]
        x_range <- range(unlist(lapply(sub_blocks, function(x) x$window$xrange)))
        y_range <- range(unlist(lapply(sub_blocks, function(x) x$window$yrange)))
        combined_window <- spatstat.geom::owin(xrange = x_range, yrange = y_range)
        
        merged_for_check <- do.call(
          spatstat.geom::superimpose,
          c(sub_blocks, list(W = combined_window))
        )
        
        if (!bkamp_window(cand_bbox, spatstat.geom::Window(merged_for_check), tol = base_tol)) {
          next
        }
        
        new_width <- diff(cand_bbox$xrange)
        new_height <- diff(cand_bbox$yrange)
        
        score_ratio <- min(new_width, new_height) / max(new_width, new_height)
        score_points <- merged_for_check$n
        
        local_candidate_merges[[length(local_candidate_merges) + 1]] <- list(
          score_ratio = score_ratio,
          score_points = score_points,
          success_idx_to_update = success_idx,
          final_ids = ids_needed_to_fill,
          final_ppp = merged_for_check
        )
      }
      
      if (length(local_candidate_merges) > 0) {
        best_score_ratio <- -1
        best_score_points <- Inf
        best_candidate_idx <- -1
        
        for (i in seq_along(local_candidate_merges)) {
          cand <- local_candidate_merges[[i]]
          
          if (cand$score_ratio > best_score_ratio + base_tol) {
            best_score_ratio <- cand$score_ratio
            best_score_points <- cand$score_points
            best_candidate_idx <- i
          } else if (abs(cand$score_ratio - best_score_ratio) < base_tol) {
            if (cand$score_points < best_score_points) {
              best_score_points <- cand$score_points
              best_candidate_idx <- i
            }
          }
        }
        
        best_merge <- local_candidate_merges[[best_candidate_idx]]
        s_idx <- best_merge$success_idx_to_update
        success_blocks[[s_idx]]$ppp <- best_merge$final_ppp
        success_blocks[[s_idx]]$ids <- best_merge$final_ids
        assigned[best_merge$final_ids] <- TRUE
        
        merge_made_in_this_pass <- TRUE
        break
      }
    }
    
    if (!merge_made_in_this_pass) break
  }
  
  list(success_blocks = success_blocks, assigned = assigned)
}

#' Stop with a descriptive error if any of `split_max`, `n_min`, `rho1`,
#' `rho2`, or `max_grid_dim` are out of range or of the wrong type.
#' @noRd
check_bkamp_block_options <- function(univariate, split_max, n_min, rho1,
                                      rho2, max_grid_dim) {
  for (nm in c("univariate")) {
    value <- get(nm)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(nm, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  scalar <- function(x, lower, finite = TRUE) is.numeric(x) && length(x) == 1L &&
    !is.na(x) && x >= lower && (!finite || is.finite(x))
  if (!is.null(split_max) && !scalar(split_max, 1)) {
    stop("split_max must be finite and at least 1.", call. = FALSE)
  }
  if (!is.null(rho1) && !scalar(rho1, 1, FALSE)) {
    stop("rho1 must be at least 1.", call. = FALSE)
  }
  if (!scalar(rho2, 1, FALSE) || (!is.null(rho1) && rho2 < rho1)) {
    stop("rho2 must be at least rho1 and at least 1.", call. = FALSE)
  }
  if (!scalar(max_grid_dim, 2) || max_grid_dim != floor(max_grid_dim)) {
    stop("max_grid_dim must be an integer of at least 2.", call. = FALSE)
  }
  if (!is.null(n_min) && (!is.numeric(n_min) ||
      length(n_min) != (if (univariate) 2L else 3L) ||
      any(!is.finite(n_min)) || any(n_min < 1))) {
    stop("n_min must contain 2 or 3 positive finite minimum counts.", call. = FALSE)
  }
  invisible(NULL)
}

#' Call `check_bkamp_marks` and additionally require a positive-area
#' rectangular observation window.
#' @noRd
check_bkamp_input <- function(ppp_obj, univariate, mark1, mark2) {
  check_bkamp_marks(ppp_obj, univariate, mark1, mark2)
  win <- spatstat.geom::Window(ppp_obj)
  if (win$type != "rectangle" || spatstat.geom::area.owin(win) <= 0) {
    stop("B-KAMP requires a positive-area rectangular observation window.", call. = FALSE)
  }
  invisible(NULL)
}

#' Validate `rvals` and `correction`, and normalise full-name aliases
#' (e.g. `"translational"` -> `"trans"`) before passing to KAMP functions.
#' @noRd
check_bkamp_inference <- function(rvals, univariate, correction) {
  if (!is.logical(univariate) || length(univariate) != 1L || is.na(univariate)) {
    stop("univariate must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(rvals) || length(rvals) == 0L || any(!is.finite(rvals)) ||
      any(rvals < 0) || anyDuplicated(rvals)) {
    stop("rvals must contain distinct nonnegative finite relative radii.", call. = FALSE)
  }
  allowed <- c("trans", "translational", "iso", "isotropic", "none")
  if (univariate) allowed <- c(allowed, "border")
  if (!is.character(correction) || length(correction) != 1L ||
      is.na(correction) || !correction %in% allowed) {
    stop("Unsupported edge correction for this analysis type.", call. = FALSE)
  }
  switch(correction, translational = "trans", isotropic = "iso", correction)
}

#' Compute unnormalized weights `n_i / p_i` for each block,
#' where `p_i` is the combined target abundance within the block.
#' @noRd
bkamp_weights <- function(blocks, univariate, mark1, mark2) {
  if (!is.list(blocks) || length(blocks) == 0L ||
      !all(vapply(blocks, inherits, logical(1), "ppp"))) {
    stop("blocks must be a nonempty list of ppp objects.", call. = FALSE)
  }
  
  vapply(seq_along(blocks), function(i) {
    pp <- blocks[[i]]
    
    check_bkamp_input(pp, univariate, mark1, mark2)
    
    mk <- spatstat.geom::marks(pp)
    m1 <- sum(mk == mark1)
    m2 <- if (univariate) 0L else sum(mk == mark2)
    
    if (pp$n < 4L || (univariate && m1 < 2L)) {
      stop("Block ", i, ": KAMP variance requires at least 4 total points and, ",
           "for univariate analysis, at least 2 target points.", call. = FALSE)
    }
    
    abundance <- (m1 + m2) / pp$n
    pp$n / abundance
  }, numeric(1))
}


#' Combine block-level z-statistics
#' (`sum(w_i * z_i) / sqrt(sum(w_i^2))`) and return a tibble with `r`, `z`, `pvalue`.
#' @noRd
bkamp_combine <- function(z, weights, rvals) {
  z <- matrix(z, nrow = length(rvals), ncol = length(weights))
  combined <- as.vector(z %*% weights) / sqrt(sum(weights^2))
  tibble::tibble(r = rvals, z = combined, pvalue = stats::pnorm(-combined))
}

#' Compute B-KAMP at one relative radius
#' 
#' @description Standardize within-block KAMP statistics and combine them.
#' 
#' @param blocks Nonempty list of marked rectangular `ppp` blocks.
#' @param rvalue A single radius
#'   In B-KAMP, this is relative to each block's shorter side and is multiplied
#'   by that length to obtain the absolute radius used for KAMP calculations.
#' @param univariate A logical value indicating whether to compute univariate KAMP (default is TRUE).
#' @param correction Type of edge correction. Defaults to translational.
#' @param mark1 Variable used to mark the points in the point pattern object for the first type.
#' @param mark2 Variable used to mark the points in the point pattern object for the second type (optional, only used if `univariate` is FALSE).
#' 
#' @returns A one-row tibble with relative radius `r`, combined `z`, and `pvalue`.
#' @keywords internal
bkamp_helper <- function(blocks, rvalue, univariate = TRUE, correction = "trans",
                           mark1,
                           mark2 = NULL) {
  correction <- check_bkamp_inference(rvalue, univariate, correction)
  
  if (length(rvalue) != 1L) stop("rvalue must be a single relative radius.", call. = FALSE)
  
  weights <- bkamp_weights(blocks, univariate, mark1, mark2)
  
  z <- vapply(seq_along(blocks), function(i) {
    
    pp <- blocks[[i]]
    radius <- rvalue * min(diff(pp$window$xrange), diff(pp$window$yrange))
    
    result <- if (univariate) {
      kamp_variance_helper(ppp_obj = pp, rvalue = radius,
                           correction = correction, mark1 = mark1)
    } else {
      kamp_variance_biv_helper(ppp_obj = pp, rval = radius,
                               correction = correction, mark1 = mark1, mark2 = mark2)
    }
    
    result$kamp / sqrt(result$var)
  }, numeric(1))
  bkamp_combine(z, weights, rvalue)
}

#' Compute B-KAMP across relative radii using KAMP's Rcpp functions
#' 
#' @description Calculate all requested radii once per block and combine block statistics.
#' 
#' @param blocks Nonempty list of marked rectangular `ppp` blocks.
#' @param rvals A vector of radii at which to calculate the KAMP variance.
#'   In B-KAMP, these are relative to each block's shorter side and are multiplied
#'   by that length to obtain the absolute radii used for KAMP calculations.
#' @param univariate A logical value indicating whether to compute univariate KAMP (default is TRUE).
#' @param correction Type of edge correction. Defaults to translational.
#' @param mark1 Variable used to mark the points in the point pattern object for the first type.
#' @param mark2 Variable used to mark the points in the point pattern object for the second type (optional, only used if `univariate` is FALSE).
#' 
#' @returns A tibble with relative radius `r`, combined `z`, and `pvalue`.
#' @keywords internal
bkamp_Rcpp_helper <- function(blocks, rvals = c(0, .05, .075, .1, .15, .2),
                                univariate = TRUE, correction = "trans",
                                mark1,
                                mark2 = NULL) {
  correction <- check_bkamp_inference(rvals, univariate, correction)
  
  weights <- bkamp_weights(blocks, univariate, mark1, mark2)
  

  order_r <- order(rvals)
  relative <- rvals[order_r]
  
  z <- vapply(seq_along(blocks), function(i) {
    pp <- blocks[[i]]
    radii <- relative * min(diff(pp$window$xrange), diff(pp$window$yrange))
    
    # The existing bivariate Rcpp wrapper passes NULL weights to C++ when
    # no pair is present. Handle the zero-variance case before that call.
    result <- if (min(spatstat.geom::nndist(pp)) > max(radii)) {
      tibble::tibble(
        kamp = rep(0, length(radii)),
        var = rep(0, length(radii))
      )
    } else if (univariate) {
      kamp_variance_Rcpp(ppp_obj = pp, rvals = radii,
                         correction = correction, mark1 = mark1)
    } else {
      kamp_variance_biv_Rcpp(ppp_obj = pp, rvals = radii,
                             correction = correction, mark1 = mark1, mark2 = mark2)
    }
    
    variance <- result$var
    variance[!is.finite(variance) | variance <= 0] <- NA_real_
    
    result$kamp / sqrt(variance)
  }, numeric(length(relative)))
  
  result <- bkamp_combine(z, weights, relative)
  
  result[match(rvals, relative), , drop = FALSE]
}
