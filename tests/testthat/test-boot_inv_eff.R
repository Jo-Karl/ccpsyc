test_that("multiplication works", {
  two_country <- dplyr::filter(example, country %in% c("NZ" , "BRA"))
  suppressWarnings(boot_eff <- boot_inv_eff(n = 10, n_sample = 200, df = two_country, group = "country", items = paste0("voice",1:3, "m")))
  expect_equal(mean(boot_eff$av), 0.210652667)
})
