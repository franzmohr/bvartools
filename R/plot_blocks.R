
# One figure per block of coefficients.
#
# The panels of a model form a grid with one row per equation and one column per
# regressor. Drawn as a single figure, a model of any size produces panels too
# small to read -- a GVAR sub-model with two dozen regressors leaves each panel
# narrower than its own axis labels, and 'layout' then fails outright with
# "figure margins too large". Each block is therefore a figure of its own, and a
# block with more columns than 'max_cols' is split over several figures.
#
# 'blocks' is a list whose elements carry a 'title', the column 'labels' of the
# block, and a 'panel' function. The latter is called once per panel with an
# index that runs over the block column by column, which is the order the panels
# of a grid are drawn in.

.plot_blocks <- function(blocks, row_names, title, max_cols = 6, lab_size = .05) {

  k <- length(row_names)

  for (block in blocks) {

    n_col <- length(block[["labels"]])
    if (n_col == 0) {
      next
    }

    # Pages of nearly equal width, so that a block of seven columns becomes two
    # figures of four and three rather than one of six and one of a single
    # column.
    n_pages <- ceiling(n_col / max_cols)
    pages <- split(seq_len(n_col), sort(rep_len(seq_len(n_pages), n_col)))

    for (page in pages) {

      n_page <- length(page)

      mat <- matrix(NA_integer_, k + 2, n_page + 1)
      mat[1, ] <- 1
      mat[-1, 1] <- c(0, 2:(k + 1))
      mat[2, -1] <- (k + 1) + 1:n_page
      mat[-(1:2), -1] <- matrix(1:(k * n_page) + k + n_page + 1, k, n_page)
      graphics::layout(mat,
                       widths = c(lab_size, rep((1 - lab_size) / n_page, n_page)),
                       heights = c(.07, lab_size, rep((1 - lab_size) / k, k)))

      # Title
      header <- block[["title"]]
      if (n_pages > 1) {
        header <- paste0(header, " (", which(vapply(pages, identical, logical(1), page)),
                         " of ", n_pages, ")")
      }
      graphics::par(mar = c(0, 0, 0, 0))
      graphics::plot.new()
      graphics::text(0.5, 0.65, labels = title, cex = 1.2)
      graphics::text(0.5, 0.25, labels = header, cex = 1)

      # Fill rows
      graphics::par(mar = c(3, 0, 0, 0))
      for (i in row_names) {
        graphics::plot.new(); graphics::text(0.5, 0.5, labels = i, adj = 0.5)
      }

      # Fill columns
      graphics::par(mar = c(0, 0, 0, 0))
      for (i in block[["labels"]][page]) {
        graphics::plot.new(); graphics::text(0.5, 0.5, labels = i, adj = 0.5)
      }

      graphics::par(mar = c(3, 2.1, .5, 1))
      for (i in page) {
        for (j in 1:k) {
          block[["panel"]]((i - 1) * k + j)
        }
      }
    }
  }

  return(invisible(NULL))
}


# The panel of one coefficient, as a path of credible bands for a time varying
# model and as a summary of the draws for a constant one.

.plot_coefficient_panel <- function(draws, tvp, type, ci_low, ci_high, show_zero_y) {

  if (tvp) {
    stats::ts.plot(t(apply(draws, 2, stats::quantile, probs = c(ci_low, .5, ci_high))),
                   xlab = "")
    if (show_zero_y) {
      graphics::abline(h = 0)
    }
  } else {
    if (type == "hist") {
      graphics::hist(draws, plot = TRUE, main = NA)
    }
    if (type == "trace") {
      stats::ts.plot(draws, xlab = "")
    }
    if (type == "boxplot") {
      graphics::boxplot(draws)
    }
  }

  return(invisible(NULL))
}
