
testthat::test_that("Data generated by makeExampleData has matching dimensions", {
    dat <- makeExampleData()
    testthat::expect_equal(nrow(dat$X), dat$n)
    testthat::expect_equal(length(dat$y), dat$n)
    testthat::expect_equal(ncol(dat$X), dat$p)
    testthat::expect_equal(length(dat$annot), dat$p)
    testthat::expect_equal(length(dat$beta), dat$p)
    testthat::expect_equal(length(dat$s), dat$p)
    testthat::expect_equal(length(dat$beta_tilde), dat$p)
    testthat::expect_equal(dat$g, length(dat$gammas))
    testthat::expect_equal(dat$g, length(dat$pis))
})


testthat::test_that("Data generated by makeExampleDataWithUnequalGroups has matching dimensions", {
    dat <- makeExampleDataWithUnequalGroups()
    testthat::expect_equal(nrow(dat$X), dat$n)
    testthat::expect_equal(length(dat$y), dat$n)
    testthat::expect_equal(ncol(dat$X), dat$p)
    testthat::expect_equal(length(dat$annot), dat$p)
    testthat::expect_equal(length(dat$beta), dat$p)
    testthat::expect_equal(length(dat$s), dat$p)
    testthat::expect_equal(length(dat$beta_tilde), dat$p)
    testthat::expect_equal(dat$g, length(dat$gammas))
    testthat::expect_equal(dat$g, length(dat$pis))
})

testthat::test_that("coef returns the coefficients with/without intercept as specified",{
    dat <- makeExampleData()
    fit <- graper(dat$X, dat$y, dat$annot)
    testthat::expect_equal(as.numeric(coef(fit, include_intercept = FALSE)), as.numeric(fit$EW_beta))
    testthat::expect_equal(as.numeric(coef(fit, include_intercept = TRUE)), c(fit$intercept, as.numeric(fit$EW_beta)))
})

testthat::test_that("getPIP returns PIPs of correct length",{
    dat <- makeExampleData()
    fit <- graper(dat$X, dat$y, dat$annot)
    testthat::expect_equal(length(as.numeric(getPIPs(fit))), dat$p)
    testthat::expect_equal(as.numeric(getPIPs(fit)), as.numeric(fit$EW_s))
})

testthat::test_that("predict returns correct number of samples",{
    dat <- makeExampleData()
    fit <- graper(dat$X, dat$y, dat$annot)
    preds <- predict(fit, dat$X)
    testthat::expect_equal(length(as.numeric(preds)), dat$n)
})

testthat::test_that("graper returns a graper object",{
    dat <- makeExampleData()
    fit <- graper(dat$X, dat$y, dat$annot)
    testthat::expect_is(fit, "graper")
})
