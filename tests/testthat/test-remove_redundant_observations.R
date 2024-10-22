test_that("Remove duplicate entries without specifying tmat", {
  gd <- data.frame(id = c(1, 1, 1, 1, 1, 1, 1),
                   state = c(1, 1, 1, 2, 2, 2, 3),
                   time = c(1, 2, 3, 4, 5, 6, 7))
  gd_rm <- remove_redundant_observations(gd)
  expect_equal(gd_rm, data.frame(id = c(1, 1, 1, 1, 1),
                                 state = c(1, 1, 2, 2, 3),
                                 time = c(1, 3, 4, 6, 7)),
               ignore_attr = "row.names")
})




test_that("Remove absorbing states with tmat, and keep without tmat", {
  gd <- data.frame(id = c(1, 1, 1, 1, 1, 1),
                   state = c(1, 3, 3, 3, 3, 3),
                   time = c(1, 2, 3, 4, 5, 6))
  tmat <- mstate::transMat(list(c(2, 3), c(3), c()))
  gd_rm1 <- remove_redundant_observations(gd)
  gd_rm2 <- remove_redundant_observations(gd, tmat)
  
  expect_equal(gd_rm1, data.frame(id = c(1, 1, 1),
                                  state = c(1, 3, 3),
                                  time = c(1, 2, 6)),
               ignore_attr = "row.names")
  expect_equal(gd_rm2, data.frame(id = c(1, 1),
                                  state = c(1, 3),
                                  time = c(1, 2)),
               ignore_attr = "row.names")
})

