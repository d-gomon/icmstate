test_that("bbase_singletime similar to JOPS::bbase", {
  ic_bbase <- bbase_singletime(0.09, 0, 1, nseg = 20, bdeg = 3)
  ic_bbase2 <- bbase_singletime(0.1, 0, 1, nseg = 20, bdeg = 3)
  ic_bbase3 <- bbase_singletime(0.09, 0, 1, nseg = 100, bdeg = 5)
  
  jops_bbase <- JOPS::bbase(0.09, 0, 1, nseg = 20, bdeg = 3)
  jops_bbase2 <- JOPS::bbase(0.1, 0, 1, nseg = 20, bdeg = 3)
  jops_bbase3 <- JOPS::bbase(0.09, 0, 1, nseg = 100, bdeg = 5)
  expect_equal(ic_bbase, jops_bbase, ignore_attr = TRUE)
  expect_equal(ic_bbase2, jops_bbase2, ignore_attr = TRUE)
  expect_equal(ic_bbase3, jops_bbase3, ignore_attr = TRUE)
})
