test_that("Non-truncated interval censored data, simple example", {
  int_test <- data.frame(L = c(1,3,2,1,4), R = c(3,6,3,3,5), cens = c(1,1,1,1,1))
  q <- supportHudgens(int_test)
  g <- q$graph
  supp <- q$support
  
  #Expected support is [2,3] & [4,5]
  expect_equal(supp, matrix(c(2,4,3,5), nrow = 2), ignore_attr = TRUE)
  expect_equal(V(g)$L, c(1,3,2,1,4))
  expect_equal(V(g)$R, c(3,6,3,3,5))
  #V(g)$name
  #V(g)$name
  #V(g)$trunc
}
)

test_that("Truncated interval censored data, simple example (Turnbull)", {
  int_turn <- data.frame(L = c(2,3,0.5,1,0.75,6), R = c(Inf,5,Inf,4,Inf,7), cens = c(0,1,0,1,0,1))
  s <- supportHudgens(int_turn)
  
  #Expected support (from Turnbull) is [1,2] and [3,4] and [6,7]
  expect_equal(s$support, matrix(c(1,3,6,2,4,7), nrow = 3), ignore_attr = TRUE)  
  
  #plot(s$graph)
}
)

test_that("Truncated interval censored data, reduction possible through Hudgens (2005) Lemma 3", {
  int_reduc <- data.frame(L = c(3,4, 1, 2), R = c(Inf, 6, Inf, 5), cens = c(0, 1, 0, 1))
  t <- supportHudgens(int_reduc)
  #We expect support (From Hudgens) [4,5] instead of [2,3] and [4,5]
  expect_equal(t$support, matrix(c(4,5), nrow = 1), ignore_attr = TRUE) 
  
  int_reduc2 <- data.frame(L = c(3,4,5,6,1,2), R = c(Inf, 7, Inf, 9, Inf, 8), cens = c(0, 1, 0, 1, 0, 1))
  t2 <- supportHudgens(int_reduc2)
  #We expect support [6,7] instead of [2,3] and [4,5] and [6,7]
  expect_equal(t2$support, matrix(c(6,7), nrow = 1), ignore_attr = TRUE) 
}
)

test_that("Existence of NPMLE (Hudgens examples)", {
  int_exis <- data.frame(L = c(3, 4, 1, 2), R = c(Inf, 6, Inf, 5), cens = c(0, 1, 0, 1), id = c(1,1,2,2))
  l <- supportHudgens(int_exis, existence = TRUE)
  #In this simple example, we expect the NPMLE to exist
  expect_true(l$exist_mle)
  #NPMLE exists if directed graph is strongly connected
  expect_equal(l$exist_mle, igraph::is_connected(l$dir_graph, mode = "strong"))
  
  #Hudgens (2005) example 1
  int_hudg1 <- data.frame(L = c(1,2,3,5,3,7), R = c(Inf, 4, Inf, 6, Inf, 8), 
                          cens = c(0, 1, 0, 1, 0, 1), id = c(1, 1, 2, 2, 3, 3))
  HE1 <- supportHudgens(int_hudg1, existence = TRUE)
  #We expect example 1 not to exist
  expect_false(HE1$exist_mle)
  
  #Hudgens (2005) example 2
  int_hudg2 <- data.frame(L = c(1, 2, 1, 5, 3, 7), R = c(Inf, 4, Inf, 6, Inf, 8), 
                          cens = c(0, 1, 0, 1, 0, 1), id = c(1, 1, 2, 2, 3, 3))
  HE2 <- supportHudgens(int_hudg2, existence = TRUE)
  #We expect example 2 to exist
  expect_true(HE2$exist_mle)
  
  #Hudgens (2005) example 3
  int_hudg3 <- data.frame(L = c(1, 2, 1, 4, 3, 6), R = c(Inf, 5, Inf, 7, Inf, 7), 
                          cens = c(0, 1, 0, 1, 0, 1), id = c(1, 1, 2, 2, 3, 3))
  HE3 <- supportHudgens(int_hudg3, existence = TRUE)
  #We expect example 2 to exist
  expect_true(HE3$exist_mle)
  
  #Hudgens (2005) example 4
  int_hudg4 <- data.frame(L = c(1, 2, 1, 4, 5, 6, 5, 8), R = c(Inf, 3, Inf, 9, Inf, 7, Inf, 9), 
                          cens = c(0, 1, 0, 1, 0, 1, 0, 1), id = c(1, 1, 2, 2, 3, 3, 4, 4))
  HE4_noreduc <- supportHudgens(int_hudg4, reduction = FALSE, existence = TRUE)
  HE4 <- supportHudgens(int_hudg4, reduction = TRUE, existence = TRUE)
  #Without reduction we don't expect the MLE to exist (directed graph not strongly connected)
  #expect_false(HE4_noreduc$exist_mle)
  #With reduction we expect the MLE to exist (directed graph strongly connected)
  #expect_true(HE4$exist_mle)
})







