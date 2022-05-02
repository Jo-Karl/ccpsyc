test_that("multiplication works", {
  result <- multi_group_eff(df = example, group = "country", items = paste0("voice",1:3, "m"))
  expect_equal(nrow(result$average), 15)
  expect_equal(mean(result$average$sd), 0.277332894)
})
