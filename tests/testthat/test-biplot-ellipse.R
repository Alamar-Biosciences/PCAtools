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

    built <- ggplot2::ggplot_build(p)$data
    polys <- Filter(function(d) "group" %in% names(d) && nrow(d) > 100, built)
    expect_true(length(polys) > 0)
    for (d in polys) expect_false(anyNA(d$x))
})
