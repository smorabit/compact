# =============================================================================
# HeatmapDistanceDiff() — diverging heatmap for a post-minus-pre distance matrix
#
# Companion to compact::HeatmapDistance(). That function puts two magnitude
# matrices on a shared SEQUENTIAL ramp, which is right for what it does. A
# difference matrix is signed and centred on zero, so it needs a DIVERGING scale
# with the neutral midpoint pinned at 0.
#
# Invariants held in every scale mode and at any limits:
#   * zero is always the neutral grey colour
#   * sign always maps to hue (blue = negative, red = positive)
#   -> an all-negative matrix stays entirely on the blue arm; the least-negative
#      cell is never rescaled into grey.
#
# Usage:  source("HeatmapDistanceDiff.R")   # after library(compact)
# =============================================================================


#' Shared colour limits across several difference matrices
#'
#' Compute one `limits` pair covering a set of difference matrices, so the same
#' colour means the same delta in every figure. Pass the result to each
#' `HeatmapDistanceDiff(limits = ...)` call.
#'
#' @param mats A list of difference matrices / data.frames.
#' @param quantile_clip Quantile of the pooled values used for the extremes.
#'   Default 1 (true min/max). Use e.g. 0.98 so one outlier module does not
#'   flatten every other figure.
#' @param symmetric If `TRUE` (default) returns `c(-M, M)` with
#'   `M = max(|lo|, hi)`, so equal colour intensity means equal magnitude on
#'   both arms. `FALSE` returns the raw `c(lo, hi)`, giving a lopsided bar that
#'   uses the full ramp.
#' @return Numeric length-2 vector, `c(lower, upper)`.
DiffLimits <- function(mats, quantile_clip = 1, symmetric = TRUE) {
  if (is.matrix(mats) || is.data.frame(mats)) mats <- list(mats)
  vals <- unlist(lapply(mats, function(m) {
    m <- as.matrix(m)
    diag(m) <- NA_real_
    as.numeric(m[upper.tri(m)])
  }), use.names = FALSE)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) stop("No finite values across the supplied matrices.")

  if (quantile_clip >= 1) {
    lo <- min(vals); hi <- max(vals)
  } else {
    p  <- (1 - quantile_clip) / 2
    lo <- as.numeric(stats::quantile(vals, p))
    hi <- as.numeric(stats::quantile(vals, 1 - p))
  }
  lo <- min(lo, 0); hi <- max(hi, 0)
  if (symmetric) {
    M <- max(abs(lo), hi)
    if (M <= 0) M <- 1e-8
    return(c(-M, M))
  }
  c(lo, hi)
}


# -----------------------------------------------------------------------------
# Internal: build the value -> colour-position mapping.
#
# Returns a transform into colour-position space, the position-space limits, the
# palette interpolation points that pin neutral at zero, and colourbar
# breaks/labels in the ORIGINAL units so a nonlinear scale self-documents.
# -----------------------------------------------------------------------------
.diff_colour_scale <- function(vals, mode, quantile_clip = 1,
                               linthresh = NULL, min_limit = 0,
                               limits = NULL, n_colours = 11L,
                               fmt = function(x) formatC(x, format = "g", digits = 3)) {

  vals <- vals[is.finite(vals)]
  if (!length(vals)) stop("No finite values to scale.")
  q <- function(x, p) as.numeric(stats::quantile(x, probs = p, na.rm = TRUE))

  fixed <- !is.null(limits)
  if (fixed) {
    if (!is.numeric(limits) || length(limits) != 2L || anyNA(limits)) {
      stop("`limits` must be a numeric vector of length 2, c(lower, upper).")
    }
    limits <- sort(limits)
    lo <- limits[1]; hi <- limits[2]
    if (lo > 0 || hi < 0) {
      stop("`limits` must bracket zero (lower <= 0 <= upper); a difference scale ",
           "has to keep the neutral colour on zero. Got c(", lo, ", ", hi, ").")
    }
    # keep a sliver of the unused arm so the palette stays strictly monotone
    if (hi <= 0) hi <-  abs(lo) / 1000
    if (lo >= 0) lo <- -abs(hi) / 1000
  }

  L_all <- if (fixed) max(abs(lo), hi) else max(q(abs(vals), quantile_clip), min_limit)
  if (!is.finite(L_all) || L_all <= 0) L_all <- 1e-8

  neg <- vals[vals < 0]
  pos <- vals[vals > 0]

  # ---- auto ------------------------------------------------------------------
  # With fixed limits the point is cross-figure comparability, which also
  # requires the same TRANSFORM in every figure. auto would be free to pick
  # differently per module, so it resolves to symmetric. Pass scale_mode
  # explicitly to override.
  if (mode == "auto") {
    if (fixed) {
      mode <- "symmetric"
      auto_note <- "auto->symmetric (fixed limits)"
    } else {
      nz <- abs(vals)[abs(vals) > 0]
      if (!length(nz)) {
        mode <- "symmetric"
        auto_note <- "auto->symmetric (all differences are zero)"
      } else {
        small <- max(q(nz, 0.10), L_all / 1e4)
        ratio <- L_all / small
        mode  <- if (ratio > 50) "symlog" else "symmetric"
        auto_note <- sprintf("auto->%s (range ratio %.0f×)", mode, ratio)
      }
    }
  } else {
    auto_note <- NULL
  }

  # ---- transform + position-space limits + breaks ----------------------------
  if (mode == "symmetric") {
    M     <- L_all
    trans <- function(x) x / M
    p_lo  <- if (fixed) lo / M else -1
    p_hi  <- if (fixed) hi / M else  1
    bvals <- unique(c(p_lo * M, p_lo * M / 2, 0, p_hi * M / 2, p_hi * M))
    note  <- "linear, symmetric"

  } else if (mode == "asymmetric") {
    Ln <- if (fixed) abs(lo) else if (length(neg)) max(q(abs(neg), quantile_clip), min_limit) else 0
    Lp <- if (fixed)     hi  else if (length(pos)) max(q(pos,      quantile_clip), min_limit) else 0
    if (Ln <= 0) Ln <- if (Lp > 0) Lp else 1e-8
    if (Lp <= 0) Lp <- Ln
    trans <- function(x) ifelse(x < 0, x / Ln, x / Lp)
    p_lo  <- -1; p_hi <- 1
    bvals <- c(-Ln, -Ln / 2, 0, Lp / 2, Lp)
    note  <- sprintf("linear, arms scaled independently (−%s / +%s) — magnitude NOT comparable across zero",
                     fmt(Ln), fmt(Lp))

  } else if (mode == "symlog") {
    if (is.null(linthresh)) linthresh <- max(q(abs(vals), 0.10), L_all / 1000)
    if (!is.finite(linthresh) || linthresh <= 0) linthresh <- L_all / 1000
    linthresh <- min(linthresh, L_all / 2)

    den   <- log1p(L_all / linthresh)
    trans <- function(x) ifelse(x == 0, 0, sign(x) * log1p(abs(x) / linthresh) / den)
    p_lo  <- if (fixed) trans(lo) else -1
    p_hi  <- if (fixed) trans(hi) else  1

    n_dec <- max(0, ceiling(log10(L_all / linthresh)))
    mags  <- unique(c(linthresh * 10^(0:n_dec), L_all))
    mags  <- sort(mags[mags <= L_all * 1.000001])
    if (length(mags) > 3) mags <- mags[round(seq(1, length(mags), length.out = 3))]
    bvals <- c(-rev(mags), 0, mags)
    note  <- sprintf("signed log, linthresh %s", fmt(linthresh))

  } else {
    stop("Unknown scale_mode: ", mode)
  }

  if (!is.null(auto_note)) note <- paste0(note, " | ", auto_note)
  if (fixed) note <- paste0(note, " | fixed limits [", fmt(lo), ", ", fmt(hi), "]")

  # ---- pin the neutral colour on zero ----------------------------------------
  # `values` is in rescaled [0,1] space over the fill limits. Placing the middle
  # palette colour at the position corresponding to 0 keeps grey on zero even
  # when the bar is lopsided.
  if (n_colours %% 2 != 1) stop("Palette must have an odd number of colours.")
  neutral01 <- (0 - p_lo) / (p_hi - p_lo)
  neutral01 <- min(max(neutral01, 1e-4), 1 - 1e-4)
  half      <- (n_colours - 1L) %/% 2L
  values01  <- c(seq(0, neutral01, length.out = half + 1L),
                 seq(neutral01, 1, length.out = half + 1L)[-1])

  # ---- colourbar breaks, dropped if outside the visible range ----------------
  bpos <- trans(bvals)
  keep <- bpos >= p_lo - 1e-9 & bpos <= p_hi + 1e-9
  # drop ticks from a nudged, effectively-unused arm (e.g. limits = c(-0.5, 0))
  keep <- keep & (bvals == 0 | abs(bvals) >= L_all / 500)
  bvals <- bvals[keep]; bpos <- bpos[keep]
  ord <- order(bpos); bvals <- bvals[ord]; bpos <- bpos[ord]
  dup <- duplicated(round(bpos, 6)); bvals <- bvals[!dup]; bpos <- bpos[!dup]

  list(
    trans    = trans,
    p_lo     = p_lo,
    p_hi     = p_hi,
    values01 = values01,
    breaks   = bpos,
    labels   = vapply(bvals, fmt, character(1)),
    mode     = mode,
    note     = note,
    L_all    = L_all
  )
}


#' Diverging heatmap of a difference matrix (post - pre)
#'
#' @param df_diff Signed difference matrix or data.frame (post - pre).
#' @param df_original Optional pre-perturbation matrix. Required when
#'   `relative = TRUE`; also used to report effect size in the subtitle.
#' @param title,subtitle Plot title / subtitle. `subtitle = NULL` auto-generates
#'   one reporting magnitude and the scale in use.
#' @param custom_order Character vector giving row/column order.
#' @param relative If `TRUE`, plot 100 * (post - pre) / pre. Needs `df_original`.
#' @param limits Optional `c(lower, upper)` fixing the colour range in plotted
#'   units. Must bracket zero. Supply the SAME value to several calls to make
#'   figures directly comparable side by side; see [DiffLimits()] to derive one
#'   pair covering a set of matrices. Values outside the range are squished to
#'   the end colour rather than dropped. When set, `quantile_clip` and
#'   `min_limit` are ignored, and `scale_mode = "auto"` resolves to
#'   `"symmetric"` so every figure shares one transform as well as one range.
#' @param scale_mode One of `"auto"`, `"symmetric"`, `"asymmetric"`, `"symlog"`.
#'   \describe{
#'     \item{symmetric}{One linear scale. Magnitude comparable everywhere.}
#'     \item{asymmetric}{Each arm linear to its own extreme, so a small positive
#'       arm stays visible beside a large negative one. Magnitude NOT comparable
#'       across zero — the bar labels both endpoints and the subtitle says so.}
#'     \item{symlog}{One signed-log scale. Comparable across sign AND across
#'       orders of magnitude — the right choice when values span e.g. -40 to
#'       +0.05, where per-arm linear scaling still crushes everything between.}
#'     \item{auto}{symlog when max|v| exceeds ~50x the 10th percentile of |v|,
#'       symmetric otherwise. Reports its choice in the subtitle.}
#'   }
#' @param linthresh symlog only: below this magnitude the scale is linear.
#'   `NULL` = 10th percentile of |v|.
#' @param min_limit Floor on the auto-fitted arm limits, in plotted units. Stops
#'   tiny differences being magnified into a dramatic-looking result. Ignored
#'   when `limits` is supplied.
#' @param quantile_clip Quantile of |value| used to set auto-fitted extremes.
#'   Default 1 (true max); e.g. 0.95 stops one outlier pair flattening the map.
#'   Ignored when `limits` is supplied.
#' @param show_values Draw the value in each cell. `NULL` = show when <= 12 groups.
#' @param digits Rounding for cell labels.
#' @param diverging_palette Colours low -> mid -> high. Must be odd length.
#'
#' @return A ggplot object.
HeatmapDistanceDiff <- function(
    df_diff,
    df_original       = NULL,
    title             = "Change in energy distance (post − pre)",
    subtitle          = NULL,
    custom_order      = NULL,
    relative          = FALSE,
    limits            = NULL,
    scale_mode        = c("auto", "symmetric", "asymmetric", "symlog"),
    linthresh         = NULL,
    min_limit         = 0,
    quantile_clip     = 1,
    show_values       = NULL,
    digits            = 3,
    diverging_palette = NULL
) {

  scale_mode <- match.arg(scale_mode)
  for (pkg in c("ggplot2", "reshape2", "scales")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.")
  }

  # ---- diverging ramp: blue <-> neutral grey <-> red -------------------------
  # Equal step count per arm, lightness monotone outward from the midpoint.
  # Poles separated by dE 18.5 under protanopia, 31.1 normal vision.
  if (is.null(diverging_palette)) {
    diverging_palette <- c(
      "#0d366b", "#1c5cab", "#2a78d6", "#86b6ef", "#cde2fb",  # negative arm
      "#f0efec",                                               # neutral zero
      "#fbd9d6", "#f3a9a4", "#e8706c", "#d03b3b", "#8f2020"    # positive arm
    )
  }
  if (length(diverging_palette) %% 2 != 1) {
    stop("`diverging_palette` must be odd length so a neutral midpoint sits at zero.")
  }

  if (is.data.frame(df_diff)) df_diff <- as.matrix(df_diff)
  if (!is.matrix(df_diff)) stop("`df_diff` must be a matrix or data.frame.")
  if (nrow(df_diff) != ncol(df_diff)) stop("`df_diff` must be square.")
  if (is.null(rownames(df_diff)) || is.null(colnames(df_diff))) {
    stop("`df_diff` must have rownames and colnames.")
  }
  if (!is.null(df_original)) {
    if (is.data.frame(df_original)) df_original <- as.matrix(df_original)
    if (!identical(dim(df_original), dim(df_diff))) {
      stop("`df_original` and `df_diff` must have the same dimensions.")
    }
  }

  # ---- ordering --------------------------------------------------------------
  if (!is.null(custom_order)) {
    miss_r <- setdiff(custom_order, rownames(df_diff))
    miss_c <- setdiff(custom_order, colnames(df_diff))
    if (length(miss_r)) stop("custom_order labels not in rows: ", paste(miss_r, collapse = ", "))
    if (length(miss_c)) stop("custom_order labels not in columns: ", paste(miss_c, collapse = ", "))
    df_diff <- df_diff[custom_order, custom_order]
    if (!is.null(df_original)) df_original <- df_original[custom_order, custom_order]
  }
  groups <- rownames(df_diff)

  # diagonal is structurally zero — drop it so it cannot anchor the eye
  diag(df_diff) <- NA_real_
  if (!is.null(df_original)) diag(df_original) <- NA_real_

  # ---- absolute or relative --------------------------------------------------
  if (relative) {
    if (is.null(df_original)) stop("`relative = TRUE` requires `df_original`.")
    denom <- df_original
    denom[!is.finite(denom) | denom == 0] <- NA_real_
    plot_mat  <- 100 * df_diff / denom
    value_lab <- "% change"
    cell_fmt  <- function(x) paste0(formatC(x, format = "f", digits = 1), "%")
    leg_fmt   <- function(x) paste0(formatC(x, format = "g", digits = 3), "%")
  } else {
    plot_mat  <- df_diff
    value_lab <- "Delta distance"
    cell_fmt  <- function(x) formatC(x, format = "f", digits = digits)
    leg_fmt   <- function(x) formatC(x, format = "g", digits = 3)
  }

  vals <- plot_mat[is.finite(plot_mat)]
  if (!length(vals)) stop("No finite values to plot.")

  sc <- .diff_colour_scale(
    vals, mode = scale_mode, quantile_clip = quantile_clip,
    linthresh = linthresh, min_limit = min_limit, limits = limits,
    n_colours = length(diverging_palette), fmt = leg_fmt
  )

  # ---- subtitle that keeps the reader honest about magnitude -----------------
  if (is.null(subtitle)) {
    bits <- paste0("max |Delta| = ", cell_fmt(max(abs(vals))))
    if (!is.null(df_original)) {
      base_med <- stats::median(df_original[is.finite(df_original)])
      max_raw  <- max(abs(df_diff[is.finite(df_diff)]))
      if (is.finite(base_med) && base_med > 0) {
        bits <- paste0(bits, sprintf("  -  %.2f%% of median baseline distance (%s)",
                                     100 * max_raw / base_med,
                                     formatC(base_med, format = "f", digits = 2)))
      }
    }
    if (!is.null(limits) && max(abs(vals)) > max(abs(limits))) {
      bits <- paste0(bits, "  |  values beyond the fixed range are clipped")
    }
    subtitle <- paste0(bits, "\nscale: ", sc$note)
  }

  # ---- upper triangle, long format -------------------------------------------
  upper <- plot_mat
  upper[lower.tri(upper, diag = TRUE)] <- NA_real_
  long <- reshape2::melt(upper, varnames = c("Var1", "Var2"), na.rm = TRUE)
  long$Var1 <- factor(as.character(long$Var1), levels = groups)
  long$Var2 <- factor(as.character(long$Var2), levels = groups)

  long$pos <- sc$trans(long$value)   # colour position
  long$lab <- cell_fmt(long$value)   # cell text is ALWAYS the true value
  # label ink stays a text token; flips to white only on the darkest cells
  arm_frac <- ifelse(long$pos < 0,
                     long$pos / ifelse(sc$p_lo == 0, -1e-9, sc$p_lo),
                     long$pos / ifelse(sc$p_hi == 0,  1e-9, sc$p_hi))
  long$ink <- ifelse(arm_frac > 0.72, "#ffffff", "#0b0b0b")

  if (is.null(show_values)) show_values <- length(groups) <= 12

  p <- ggplot2::ggplot(long, ggplot2::aes(x = Var2, y = Var1)) +
    ggplot2::geom_tile(ggplot2::aes(fill = pos), colour = "white", linewidth = 0.6) +
    ggplot2::scale_fill_gradientn(
      colours  = diverging_palette,
      values   = sc$values01,
      limits   = c(sc$p_lo, sc$p_hi),
      breaks   = sc$breaks,
      labels   = sc$labels,
      oob      = scales::squish,
      na.value = "transparent",
      name     = value_lab
    ) +
    ggplot2::scale_x_discrete(limits = groups, drop = FALSE) +
    ggplot2::scale_y_discrete(limits = groups, drop = FALSE) +
    ggplot2::coord_fixed() + # NOTE let ggplot stretch the heatmap panel to better use the available width/height, so it should reduce a lot of the empty space. The tradeoff is that the tiles may no longer be perfectly square—they can become rectangular depending on your PDF dimensions.
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x         = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1,
                                                  size = 9, colour = "#52514e"),
      axis.text.y         = ggplot2::element_text(size = 9, colour = "#52514e"),
      plot.title          = ggplot2::element_text(size = 12, colour = "#0b0b0b", hjust = 0),
      plot.subtitle       = ggplot2::element_text(size = 8, colour = "#898781", hjust = 0),
      plot.title.position = "plot",
      legend.title        = ggplot2::element_text(size = 8.5, colour = "#52514e"),
      legend.text         = ggplot2::element_text(size = 8,   colour = "#898781"),
      panel.grid          = ggplot2::element_blank(),
      panel.border        = ggplot2::element_blank(),
      panel.background    = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background     = ggplot2::element_rect(fill = "white", colour = NA)
    )

  if (show_values) {
    p <- p +
      ggplot2::geom_text(ggplot2::aes(label = lab, colour = ink),
                         size = 2.5, show.legend = FALSE) +
      ggplot2::scale_colour_identity()
  }

  p
}


#' Pre / post / difference triptych
#'
#' Sequential shared scale on pre and post (via compact::HeatmapDistance),
#' diverging zero-centred scale on the difference.
#'
#' @param df_pre,df_post Distance matrices.
#' @param label Perturbation name, used in the post panel title.
#' @param custom_order Row/column order.
#' @param custom_palette Sequential palette passed to HeatmapDistance().
#' @param min_val,max_val Optional fixed limits for the pre/post panels, so
#'   those are also comparable across figures.
#' @param ... Passed to HeatmapDistanceDiff() (limits, scale_mode, ...).
#' @return A patchwork object.
HeatmapDistanceTriptych <- function(df_pre, df_post, label = "",
                                    custom_order = NULL,
                                    custom_palette = NULL,
                                    min_val = NULL, max_val = NULL, ...) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }
  if (is.data.frame(df_pre))  df_pre  <- as.matrix(df_pre)
  if (is.data.frame(df_post)) df_post <- as.matrix(df_post)

  pair <- HeatmapDistance(
    df_original     = df_pre,
    df_perturbed    = df_post,
    title_original  = "Pre-perturbation",
    title_perturbed = paste0("Post: ", label),
    custom_palette  = custom_palette,
    custom_order    = custom_order,
    min_val         = min_val,
    max_val         = max_val
  )

  diff_panel <- HeatmapDistanceDiff(
    df_post - df_pre,
    df_original  = df_pre,
    title        = "Difference (post − pre)",
    custom_order = custom_order,
    ...
  )

  patchwork::wrap_plots(pair, diff_panel, nrow = 1, widths = c(2, 1))
}
