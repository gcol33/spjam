# Tests for formula parsing

test_that("spjam_formula creates valid object", {
  f <- spjam_formula(y ~ x1 + x2)
  expect_s3_class(f, "spjam_formula")
  expect_equal(f$ecological, y ~ x1 + x2, ignore_attr = TRUE)
  expect_equal(deparse(f$sampling), "~1")
})

test_that("spjam_formula accepts sampling formula", {
  f <- spjam_formula(y ~ x1, sampling = ~ z1 + z2)
  expect_s3_class(f, "spjam_formula")
  expect_equal(f$sampling, ~ z1 + z2, ignore_attr = TRUE)
})

test_that("spjam_formula rejects non-formula input", {
  expect_error(spjam_formula("y ~ x"), "must be a formula")
  expect_error(spjam_formula(y ~ x, sampling = "z ~ 1"), "must be a formula")
})
