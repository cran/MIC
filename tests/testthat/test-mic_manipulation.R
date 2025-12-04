test_that("test get_mic", {
  meta_dt <- data.frame(genome_id = c("a1",
                                      "b2",
                                      "f6",
                                      "c1"),
                        gent_mic = suppressWarnings(
                          AMR::as.mic(c("4", "<0.5", "1.5", "NA"))
                          )
                        )
  ids <- c("c1", "b2", "f6", "n4")
  expected_mics <- suppressWarnings(
    AMR::as.mic(c("NA", "<0.5", "1.5", "NA"))
    )
  expect_equal(get_mic(meta_dt,
                       ids,
                       ab_col = "gent_mic",
                       id_col = "genome_id",
                       simplify = TRUE),
               expected_mics)

})

test_that("test mic_censor", {
  censor_rules <- list("B_ESCHR_COLI" = list(
    "AMK" = list(min = 2, max = 32),
    "CHL" = list(min = 4, max = 64),
    "GEN" = list(min = 1, max = 16),
    "CIP" = list(min = 0.015, max = 4),
    "MEM" = list(min = 0.016, max = 16),
    "AMX" = list(min = 2, max = 64),
    "AMC" = list(min = 2, max = 64),
    "FEP" = list(min = 0.5, max = 64),
    "CAZ" = list(min = 1, max = 128),
    "TGC" = list(min = 0.25, max = 1)
  ))

  suppressMessages(
    {
      expect_s3_class(mic_censor(AMR::as.mic(0.5), "TGC", "B_ESCHR_COLI", censor_rules), "mic")

      expect_equal(mic_censor(AMR::as.mic(128), "AMK", "B_ESCHR_COLI", censor_rules),
                   AMR::as.mic(">32"))
      expect_equal(mic_censor(AMR::as.mic(">128"), "AMK", "B_ESCHR_COLI", censor_rules),
                   AMR::as.mic(">32"))
      expect_equal(mic_censor(">128", "AMK", "B_ESCHR_COLI", censor_rules),
                   AMR::as.mic(">32"))

      # test vectorized
      expect_equal(mic_censor(c(AMR::as.mic(128), AMR::as.mic(0.5)), "AMK", "B_ESCHR_COLI", censor_rules),
                   c(AMR::as.mic(">32"), AMR::as.mic("<=2")))
    }
  )

  # test without rules
  expect_equal(mic_censor(AMR::as.mic(128), max = "32"), AMR::as.mic(">32"))
  expect_equal(mic_censor(AMR::as.mic(0.5), min = "2"), AMR::as.mic("<=2"))
})

test_that("test mic_uncensor", {
  expect_equal(mic_uncensor(">16", method = "scale"), AMR::as.mic(32))
  expect_equal(mic_uncensor(">16", method = "simple"), AMR::as.mic(16))

  expect_equal(mic_uncensor(c("0.5", "<1"), method = "scale"), AMR::as.mic(c("0.5", "0.5")))
  expect_equal(mic_uncensor(c("0.5", "<1"), method = "simple"), AMR::as.mic(c("0.5", "1")))

  ab <- "GEN"
  mo <- "Escherichia coli"

  expect_true(suppressMessages(as.numeric(
    mic_uncensor(">16", method = "bootstrap", ab = ab, mo = mo)) > 16))

  mic1 <- c("0.5", "4", ">8", "2", "<=1", ">16")
  uncensored_mic1 <- suppressMessages(
    mic_uncensor(mic1, method = "bootstrap", ab = ab, mo = mo))

  expect_true(as.numeric(uncensored_mic1[3]) > 8)

  expect_equal(mic_uncensor(NA, method = "bootstrap", ab = ab, mo = mo), AMR::NA_mic_)

  expect_warning(suppressMessages(
    mic_uncensor("<0.004", method = "bootstrap", ab = ab, mo = mo)))
})

test_that("test force_mic", {
  mode <- "closest"
  expect_equal(force_mic("2", method = mode),
               "2")
  expect_equal(force_mic("2.1", method = mode),
               "2")
  expect_equal(force_mic("2.5", method = mode),
               "2")
  expect_equal(force_mic("2.9", method = mode),
               "2")
  expect_equal(force_mic("3", method = mode, prefer = "max"),
               "4")

  mode <- "round up"
  expect_equal(force_mic("2", method = mode),
               "2")
  expect_equal(force_mic("2.1", method = mode),
               "4")
  expect_equal(force_mic("2.5", method = mode),
               "4")
  expect_equal(force_mic("2.9", method = mode),
               "4")
  expect_equal(force_mic("3", method = mode, prefer = "max"),
               "4")

  # keep censoring
  mode <- "closest"
  expect_equal(force_mic("<0.5", method = mode, leq = NULL, geq = NULL),
               "<0.5")
  expect_equal(force_mic("<===0.5", method = mode, leq = NULL, geq = NULL),
               "<=0.5")
  expect_equal(force_mic(">0.5", method = mode, leq = NULL, geq = NULL),
               ">0.5")
  expect_equal(force_mic(">===0.5", method = mode, leq = NULL, geq = NULL),
               ">=0.5")
  expect_equal(force_mic("<0.53", method = mode, leq = NULL, geq = NULL),
               "<0.53")
  expect_equal(force_mic("<0.53", method = "round up", leq = NULL, geq = NULL),
               "<0.53")

  # leq and geq
  mode <- "closest"
  expect_equal(force_mic("<=0.5", method = mode, leq = TRUE),
               "<=0.5")
  expect_equal(force_mic(">0.5", method = mode, geq = TRUE),
               ">=0.5")
  expect_equal(force_mic("<0.5", method = mode, leq = TRUE),
               "<=0.5")
  expect_equal(force_mic(">=0.5", method = mode, geq = FALSE),
               ">0.5")

  mode <- "round up"
  expect_equal(force_mic("<=0.5", method = mode, leq = TRUE),
               "<=0.5")
  expect_equal(force_mic(">0.5", method = mode, geq = TRUE),
               ">=0.5")
  expect_equal(force_mic("<0.5", method = mode, leq = TRUE),
               "<=0.5")
  expect_equal(force_mic(">=0.5", method = mode, geq = FALSE),
               ">0.5")
})

test_that("test mic_s_breakpoint", {
  suppressMessages(
    expect_s3_class(mic_s_breakpoint(mo = "E. coli", ab = "TZP", accept_ecoff = FALSE), "mic")
  )
  suppressMessages(
    expect_s3_class(mic_s_breakpoint(mo = "E. coli", ab = "CHL", accept_ecoff = TRUE), "mic")
  )

  expect_error(
    suppressMessages(
      mic_s_breakpoint(mo = "E. coli", ab = "CHL", accept_ecoff = FALSE)
      )
    )
})

test_that("test mic_r_breakpoint", {
  suppressMessages(
    expect_s3_class(mic_r_breakpoint(mo = "E. coli", ab = "TZP", accept_ecoff = FALSE), "mic")
  )
  suppressMessages(
    expect_s3_class(mic_r_breakpoint(mo = "E. coli", ab = "CHL", accept_ecoff = TRUE), "mic")
  )

  expect_error(
    suppressMessages(
      mic_r_breakpoint(mo = "E. coli", ab = "CHL", accept_ecoff = FALSE)
      )
    )
})

test_that("test fill dilution levels", {
  test_range <- AMR::as.mic(c("0.5", "4", "32"))
  test_range <- droplevels(test_range, as.mic = TRUE)
  filled_range <- fill_dilution_levels(test_range)
  expect_true(all.equal(as.character(test_range),
                        as.character(filled_range)))
  expect_false(all(levels(filled_range) %in% levels(test_range)))
  expect_true(length(levels(filled_range)) >= length(levels(test_range)))

  start <- AMR::as.mic(c("0.002","0.008","0.016","0.03","0.06","0.125","0.25","0.5","1","2","4","8","16","32"))
  target <- AMR::as.mic(c("0.002", "0.004", "0.008","0.016","0.03","0.06","0.125","0.25","0.5","1","2","4","8","16","32"))
  filled_range <- fill_dilution_levels(start)
  expect_true(all.equal(levels(as.character(target)),
                        levels(as.character(filled_range))))
})
