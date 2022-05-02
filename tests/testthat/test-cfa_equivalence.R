test_that("Extracted CFI works", {
  model <- "voice =~ voice1m + voice2m + voice3m
            help =~ help1m + help2m + help3m"
  expect_equal(round(equival(x = model, dat = example, group = "country")$DELTACFI,3) , c(NA, 0.010, 0.029))
})
