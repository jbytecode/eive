library(testthat)
library(eive)

test_that("eive.cga classical example", {
    
    require("eive")

    testdata <- generate.eive.data(
        n = 100,
        e.sd = 1,
        delta.sd = 1,
        seed = 12345,
        useotherx = TRUE
    )
    
    expect_equal(dim(testdata), c(100, 3))
    expect_true(
        all.equal(
            names(testdata), 
            c("xdelta", "xother", "y")
        )
    )
})
