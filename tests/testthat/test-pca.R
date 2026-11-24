# Tests the pca() function.
# library(testthat); library(PCAtools); source("setup.R"); source("test-pca.R")

test_that("pca settings work as expected", {
    basic <- pca(lcounts)
    expect_equal(sum(basic$variance), 100)
    expect_equal(as.matrix(basic$rotated), prcomp(t(lcounts))$x)

    scaled <- pca(lcounts, scale=TRUE)
    expect_equal(as.matrix(scaled$rotated), prcomp(t(lcounts), scale=TRUE)$x)

    not_centered <- pca(lcounts, center=FALSE, rank=10)
    expect_equal(
      as.matrix(not_centered$rotated),
      prcomp(t(lcounts), center=FALSE, rank=10)$x
    )
    
    trunc <- pca(lcounts, rank=10)
    expect_equal(as.matrix(trunc$rotated), prcomp(t(lcounts), rank=10)$x)

    expect_identical(pca(t(lcounts), transposed=TRUE), basic)
})

test_that("removal of low-variance genes works", {
    v <- apply(lcounts, 1, var)
    keep <- order(v, decreasing=TRUE)[seq_len(nrow(lcounts)/2)]
    ref <- pca(lcounts[keep,])
    alt <- pca(lcounts, removeVar=0.5)
    expect_equal(ref$variance, alt$variance)

    # Works properly at extremes.
    ref <- pca(lcounts)
    alt <- pca(lcounts, removeVar=0)
    expect_equal(ref$variance, alt$variance)

    ref <- pca(lcounts[which.max(v),,drop=FALSE])
    alt <- pca(lcounts, removeVar=1)
    expect_equal(ref$variance, alt$variance)
})

test_that("percentage of variance calculations are correct", {
    basic <- pca(lcounts, rank=min(dim(lcounts)))
    expect_equal(sum(basic$variance), 100)

    scaled <- pca(lcounts, rank=min(dim(lcounts)), scale=TRUE)
    expect_equal(sum(scaled$variance), 100)
    expect_false(isTRUE(all.equal(scaled, basic)))

    not_centered <- pca(lcounts, center=FALSE, rank=10)
    expect_equal(
      not_centered$variance,
      head(summary(prcomp(t(lcounts), center=FALSE, rank=10))$importance[2,], 10) * 100,
      tolerance=1e-4, # summary.prcomp rounds to 5th decimal place.
      check.attributes=TRUE
    )

    removed <- pca(lcounts, rank=min(dim(lcounts)), removeVar=0.5)
    expect_equal(sum(removed$variance), 100)
    expect_equal(nrow(removed$loadings), nrow(lcounts)/2)

    both <- pca(lcounts, rank=min(dim(lcounts)), removeVar=0.5, scale=TRUE)
    expect_equal(sum(both$variance), 100)
    expect_equal(nrow(both$loadings), nrow(lcounts)/2)
    expect_false(isTRUE(all.equal(both, removed)))
})

test_that("loading PCAtools does not load SparseArray", {
    # SparseArray 1.10.x leaves entries on R's protect stack when it loads.
    # The child process can only load an installed build, so check that build.
    skip_if(requireNamespace("pkgload", quietly=TRUE) &&
        pkgload::is_dev_package("PCAtools"))
    lib <- dirname(getNamespaceInfo("PCAtools", "path"))
    code <- sprintf('loadNamespace("PCAtools", lib.loc="%s"); cat("SparseArray" %%in%% loadedNamespaces())', lib)
    out <- system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(code)),
        stdout=TRUE)
    expect_identical(tail(out, 1), "FALSE")
})

test_that("a constant variable adds no PC when scaled", {
    x <- lcounts
    x[1, ] <- 5.123
    scaled <- pca(x, scale=TRUE)
    nonnull <- scaled$sdev > 1e-8
    expect_equal(unname(unlist(scaled$loadings[1, nonnull])), rep(0, sum(nonnull)))
    expect_equal(sum(scaled$variance), 100 * (nrow(x) - 1) / nrow(x))
})

test_that("rank above the number of components warns and is capped", {
    expect_warning(capped <- pca(lcounts, rank=1e6), "more singular values")
    expect_equal(capped, pca(lcounts))
})

test_that("pca works with alternative matrix representations", {
    library(Matrix)
    X <- as(lcounts, "dgeMatrix")
    ref <- pca(lcounts, rank=5)
    basic <- pca(X, rank=5)
    expect_equal(ref, basic)
})
