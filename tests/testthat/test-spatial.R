# Tests for spatial specifications

test_that("spatial_gp creates valid object", {
  sp <- spatial_gp()
  expect_s3_class(sp, "spjam_spatial")
  expect_s3_class(sp, "spjam_spatial_gp")
  expect_equal(sp$type, "gp")
  expect_equal(sp$nu, 1.5)
})

test_that("spatial_spde creates valid object", {
  sp <- spatial_spde()
  expect_s3_class(sp, "spjam_spatial")
  expect_s3_class(sp, "spjam_spatial_spde")
  expect_equal(sp$type, "spde")
})

test_that("spatial_car creates valid object", {
  sp <- spatial_car(type = "icar")
  expect_s3_class(sp, "spjam_spatial")
  expect_s3_class(sp, "spjam_spatial_car")
  expect_equal(sp$type, "car")
  expect_equal(sp$car_type, "icar")
})

test_that("spatial_bym2 creates valid object", {
  sp <- spatial_bym2()
  expect_s3_class(sp, "spjam_spatial")
  expect_s3_class(sp, "spjam_spatial_bym2")
  expect_equal(sp$type, "bym2")
})
