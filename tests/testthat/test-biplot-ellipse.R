# A small, wide group has an ellipse that extends past the axis limits.
# The polygon must stay whole instead of being split by NA vertices.

test_that("biplot ellipse is not split by axis limits", {
    set.seed(1)
    m <- cbind(
        matrix(rnorm(2 * 300, sd = c(6, 5)), nrow = 2),
        matrix(rnorm(2 * 6, mean = c(-8, 3), sd = c(8, 1.5)), nrow = 2))
    colnames(m) <- paste0("s", seq_len(ncol(m)))
    rownames(m) <- c("g1", "g2")
    meta <- data.frame(grp = rep(c("A", "B"), c(300, 6)), row.names = colnames(m))

    pc <- pca(m, metadata = meta, center = FALSE)
    p <- biplot(pc, colby = "grp", ellipse = TRUE, ellipseType = "norm",
        xlim = c(-20, 20), ylim = c(-12, 12))

    built <- suppressWarnings(ggplot2::ggplot_build(p))
    is_poly <- vapply(built$plot$layers,
        function(l) inherits(l$geom, "GeomPolygon"), logical(1))
    polys <- built$data[is_poly]
    expect_length(polys, 1)
    expect_equal(length(unique(polys[[1]]$group)), 2)
    expect_false(anyNA(polys[[1]]$x))
    expect_false(anyNA(polys[[1]]$y))
})

test_that("biplot ellipse is not split by y-axis limits", {
    set.seed(3)
    m <- cbind(
        matrix(rnorm(2 * 300, sd = c(6, 5)), nrow = 2),
        matrix(rnorm(2 * 6, mean = c(3, -8), sd = c(1.5, 8)), nrow = 2))
    colnames(m) <- paste0("s", seq_len(ncol(m)))
    rownames(m) <- c("g1", "g2")
    meta <- data.frame(grp = rep(c("A", "B"), c(300, 6)), row.names = colnames(m))

    pc <- pca(m, metadata = meta, center = FALSE)
    p <- biplot(pc, colby = "grp", ellipse = TRUE, ellipseType = "norm",
        xlim = c(-30, 30), ylim = c(-12, 12))

    built <- suppressWarnings(ggplot2::ggplot_build(p))
    is_poly <- vapply(built$plot$layers,
        function(l) inherits(l$geom, "GeomPolygon"), logical(1))
    expect_false(anyNA(built$data[is_poly][[1]]$y))
})

test_that("biplot keeps reversed limits reversed", {
    set.seed(2)
    m <- matrix(rnorm(2 * 50, sd = 4), nrow = 2,
        dimnames = list(c("g1", "g2"), paste0("s", 1:50)))
    pc <- pca(m, center = FALSE)
    p <- biplot(pc, lab = NULL, xlim = c(5, -5), ylim = c(-5, 5))
    sc <- suppressWarnings(ggplot2::ggplot_build(p))$layout$panel_scales_x[[1]]
    tr <- if (is.function(sc$get_transformation)) sc$get_transformation() else sc$trans
    expect_equal(tr$name, "reverse")
})

test_that("biplot limits survive a later coord and clip at the panel", {
    set.seed(2)
    m <- matrix(rnorm(2 * 50, sd = 4), nrow = 2,
        dimnames = list(c("g1", "g2"), paste0("s", 1:50)))
    pc <- pca(m, center = FALSE)
    p <- biplot(pc, lab = NULL, xlim = c(-5, 5), ylim = c(-5, 5))
    expect_equal(p$coordinates$clip, "on")

    rng <- function(p) {
        pp <- suppressWarnings(ggplot2::ggplot_build(p))$layout$panel_params[[1]]
        list(x = pp$x.range, y = pp$y.range)
    }
    expect_equal(rng(p)$x, c(-5.5, 5.5))
    expect_equal(rng(p + ggplot2::coord_fixed())$x, c(-5.5, 5.5))
    expect_equal(rng(p + ggplot2::coord_fixed())$y, c(-5.5, 5.5))
    expect_equal(rng(p + ggplot2::coord_flip())$y, c(-5.5, 5.5))
})
