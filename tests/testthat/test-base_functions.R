
#Testing subseteq functions + benchmarking

A <- matrix(c(1,3,4,7), ncol = 2, byrow = TRUE)

B <- matrix(c(0,4,1,3,1,4,0,3,5,8), ncol = 2, byrow = TRUE)

test_that("subset closed/open intervals", {
  expect_equal(AsubsetB(A, B), matrix(c(rep(1, 4), rep(0, 6)), nrow = 2, byrow = TRUE))
  expect_equal(AsubsetB(A, B, B.left.open = TRUE), matrix(c(c(1, 0, 0, 1), rep(0, 6)), nrow = 2, byrow = TRUE))
  expect_equal(AsubsetB(A, B, B.right.open = TRUE), matrix(c(c(1, 0, 1, 0), rep(0, 6)), nrow = 2, byrow = TRUE))
  expect_equal(AsubsetB(A, B, B.left.open = TRUE, B.right.open = TRUE), matrix(c(c(1, 0, 0, 0), rep(0, 6)), nrow = 2, byrow = TRUE))
})

AsubsetB(A, B)

# library(rbenchmark)
# benchmark("simple" = {
#   AsubsetB(A, B)
# }, replications = 10000)



#Testing larger eq functions + benchmarking
A <- matrix(c(1, 3, 2, 4, 5, 6), ncol = 2, byrow = TRUE)
b <- c(0, 3, 2, 1)

test_that("Larger equal subset infinite interval", {
  expect_equal(Alargerb(A, b), matrix(c(1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1), ncol = 4, byrow = TRUE))
})

Alargerb(A, b)


#Testing event contained in half-open interval + benchmarking

a <- c(1,4,2,5)

B <- matrix(c(1, 2, 2, 5, 4, 7), ncol = 2, byrow = TRUE)

test_that("Event time in half open interval", {
  expect_equal(ainB(a, B), matrix(c(rep(0, 4), 1, 0, 1, 0, 0, 0, 1, 1), ncol = 3, byrow = TRUE))
})

ainB(a, B)



#Testing event times larger equal other event times


a <- c(1,5, 4, 2)
b <- c(1, 3, 7, 6)

test_that("Event time greater equal other event time", {
  expect_equal(ageqb(a, b), matrix(c(1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0), byrow = TRUE, ncol = 4))
})

ageqb(a, b)



#Test intersection half-open and event times

A <- matrix(c(1, 3, 4, 7, 2, 5), ncol = 2, byrow = TRUE)
b <- c(0, 1, 2, 5, 6, 9)


test_that("Intersection", {
  expect_equal(Aintersectb(A, b), c(1,2, 5, 6))
  expect_equal(Aintersectb(A, b, A.left.open = TRUE), c(2, 5, 6))
})

Aintersectb(A, b)


# library(rbenchmark)
# benchmark("simple" = {
#   Aintersectb(A, b)
# }, 
# replications = 1000)






