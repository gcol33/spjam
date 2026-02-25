# Tests for family specifications

test_that("spjam_family creates valid object", {
  fam <- spjam_family("poisson")
  expect_s3_class(fam, "spjam_family")
  expect_equal(fam$family, "poisson")
  expect_equal(fam$link, "log")
})

test_that("spjam_family uses correct default links", {
  expect_equal(spjam_family("poisson")$link, "log")
  expect_equal(spjam_family("binomial")$link, "logit")
  expect_equal(spjam_family("negbin")$link, "log")
})

test_that("spjam_family allows custom link", {
  fam <- spjam_family("binomial", link = "probit")
  expect_equal(fam$link, "probit")
})

test_that("spjam_family detects zero-inflation", {
  expect_false(spjam_family("poisson")$zero_inflated)
  expect_true(spjam_family("zi_poisson")$zero_inflated)
  expect_true(spjam_family("zi_negbin")$zero_inflated)
})
