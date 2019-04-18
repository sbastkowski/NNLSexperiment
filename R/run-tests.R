library(testthat)
source("Tests/tests.R")

test_results <- test_dir("Tests", reporter="summary")
