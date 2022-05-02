test_that("Extracted CFI works", {
  model <- "voice =~ voice1m + voice2m + voice3m
            help =~ help1m + help2m + help3m"
  expect_equal(equival(x = model, dat = example, group = "country")$DELTACFI , c(NA, 0.009853305, 0.028655287))
})
