context("frame utils tests")

library(checkmate)

test_that("empty contsant_matrix() and frame conversion are correct", {
    emptymtx <- constant_matrix(0, list(orange = character(0),
                                        apple = character(0)))
    expect_equal(names(dimnames(emptymtx)), c("orange", "apple"))
    expect_equal(dim(emptymtx), c(0L, 0L))
    emptydf <- matrix2frame(emptymtx)
    expect_data_frame(emptydf, nrows = 0)
    expect_names(colnames(emptydf), must.include=c("orange", "apple", "w"))
})
