test_that("test compare mic", {
  mic1 <- AMR::as.mic(c("4","8", ">32"))
  mic2 <- AMR::as.mic(c("4", "4", "2"))
  validation <- compare_mic(mic1, mic2)
  expect_s3_class(validation, "mic_validation")
})

test_that("test essential agreement", {
  expect_true(essential_agreement(2, 2))
  expect_true(essential_agreement(2, 4))
  expect_true(essential_agreement(2, ">2"))
  expect_true(essential_agreement("2", "<2"))
  expect_true(essential_agreement(0.002, "<0.002"))
  expect_true(essential_agreement("0.001", "<0.002"))

  expect_false(essential_agreement(2, ">4"))
  expect_false(essential_agreement(256, ">512"))
  expect_false(essential_agreement(128, ">256"))
  expect_false(essential_agreement(512, "128"))
  expect_false(essential_agreement(0.002, "0.006"))
  expect_error(essential_agreement(2, ">1024"))
  expect_false(essential_agreement(2, "<1024"))


  left_check <- c(4, 2, 7, ">32")
  right_check <- c(4, 1, 1, "32")
  expect_equal(essential_agreement(left_check, right_check),
               c(TRUE, TRUE, FALSE, TRUE))

  # this section tests tolerance
  expect_true(essential_agreement("<=2", "0.002", tolerate_censoring = "x"))
  expect_true(essential_agreement("0.002", "<=2", tolerate_censoring = "y"))
  expect_true(essential_agreement("<=2", "0.002", tolerate_censoring = "both"))
  expect_true(essential_agreement("0.002", "<=2", tolerate_censoring = "both"))
  expect_false(essential_agreement("<=2", "0.002", tolerate_censoring = "y"))
  expect_false(essential_agreement("<=2", "0.002", tolerate_censoring = "strict"))

  expect_true(essential_agreement(">32", "256", tolerate_censoring = "x"))
  expect_false(essential_agreement(">32", "256", tolerate_censoring = "y"))
  expect_false(essential_agreement(">32", "256", tolerate_censoring = "strict"))

  expect_true(essential_agreement(">32", ">128", tolerate_censoring = "x"))
  expect_false(essential_agreement(">32", ">16", tolerate_censoring = "x"))
  expect_true(suppressWarnings(
    essential_agreement(">32", ">16", tolerate_censoring = "both")))
  expect_warning(essential_agreement(">32", ">16", tolerate_censoring = "both"))
  expect_true(suppressMessages(
    is.na(essential_agreement(">32", ">16", tolerate_censoring = "strict"))))

  expect_false(essential_agreement(">32", ">128", tolerate_censoring = "y"))
  expect_true(is.na(suppressMessages(
    essential_agreement(">32", ">128", tolerate_censoring = "strict"))))

  expect_true(essential_agreement("<=2", "<=0.125", tolerate_censoring = "x"))
  expect_false(essential_agreement("<=2", "<=4", tolerate_censoring = "x"))
  expect_false(essential_agreement("<=2", "<=0.125", tolerate_censoring = "y"))
  expect_true(is.na(suppressMessages(
    essential_agreement("<=2", "<=0.125", tolerate_censoring = "strict"))))

  expect_false(essential_agreement(">4", "<=0.5", tolerate_censoring = "both"))
  expect_false(essential_agreement(">4", "<=0.5", tolerate_censoring = "x"))
  expect_false(essential_agreement(">4", "<=0.5", tolerate_censoring = "y"))
  expect_false(essential_agreement(">4", "<=0.5", tolerate_censoring = "strict"))
  expect_false(essential_agreement(">4", "<=4", tolerate_censoring = "both"))

  expect_false(essential_agreement("2", "<=0.001", tolerate_censoring = "both"))
  expect_false(essential_agreement("2", "<=16", tolerate_censoring = "x"))
  expect_true(essential_agreement("2", "<=16", tolerate_censoring = "y"))

  expect_false(essential_agreement("2", ">32", tolerate_censoring = "both"))
  expect_false(essential_agreement("2", ">32", tolerate_censoring = "x"))
  expect_false(essential_agreement("2", ">32", tolerate_censoring = "y"))

  expect_false(essential_agreement(">2", "2", tolerate_matched_censoring = "strict"))
  expect_true(essential_agreement(">2", "2", tolerate_matched_censoring = "both"))
  expect_true(essential_agreement(">2", "2", tolerate_matched_censoring = "x"))
  expect_false(essential_agreement(">2", "2", tolerate_matched_censoring = "y"))

  expect_false(essential_agreement("<=0.5", "0.5", tolerate_matched_censoring = "strict"))
  expect_true(essential_agreement("<=0.5", "0.5", tolerate_matched_censoring = "both"))
  expect_true(essential_agreement("<=0.5", "0.5", tolerate_matched_censoring = "x"))
  expect_false(essential_agreement("<=0.5", "0.5", tolerate_matched_censoring = "y"))

  # test equality censors (usually on lower end)
  expect_true(essential_agreement("<=0.5", "1", tolerate_censoring = "x"))
  expect_false(essential_agreement("<=0.5", "1", tolerate_censoring = "y"))
  expect_false(essential_agreement("<=0.5", "1", tolerate_censoring = "strict"))
  expect_false(essential_agreement("<0.5", "1", tolerate_censoring = "both"))
  expect_true(essential_agreement("<=0.5", "0.75", tolerate_censoring = "x"))
  expect_false(essential_agreement("<=0.5", "2", tolerate_censoring = "x"))

  expect_false(essential_agreement(">4", "2", tolerate_censoring = "both"))
  expect_true(essential_agreement(">=4", "2", tolerate_censoring = "both"))

  # test with y censoring
  expect_true(essential_agreement("1", "<=0.5", tolerate_censoring = "both", tolerate_leq = TRUE, tolerate_geq = TRUE))
  expect_false(essential_agreement("1", "<0.5", tolerate_censoring = "both", tolerate_leq = TRUE, tolerate_geq = TRUE))
  expect_false(essential_agreement("1", "<=0.5", tolerate_censoring = "strict", tolerate_leq = TRUE, tolerate_geq = TRUE))
  expect_true(essential_agreement("0.75", "<=0.5", tolerate_censoring = "both", tolerate_leq = TRUE, tolerate_geq = TRUE))
  expect_false(essential_agreement("2", "<=0.5", tolerate_censoring = "both", tolerate_leq = TRUE, tolerate_geq = TRUE))

})

test_that("test compare_mic categorical", {
  gs <- c("0.5", "4", ">8", "2")
  test <- c("0.5", "8", "2", "0.5")
  expect_s3_class(compare_mic(gs, test), "mic_validation")
  expect_equal(summary(compare_mic(gs, test))[["EA_pcent"]], 0.5)

  ab <- c("amoxicillin", "amoxicillin", "gentamicin", "gentamicin")
  mo <- c("Escherichia coli", "Escherichia coli", "Proteus mirabilis", "Proteus mirabilis")
  suppressMessages(
    {
      expect_s3_class(compare_mic(gs, test, ab, mo), "mic_validation")
      expect_equal(summary(compare_mic(gs, test, ab, mo))[[1, "EA_pcent"]], 1)
      expect_equal(summary(compare_mic(gs, test, ab, mo))[[2, "EA_pcent"]], 0)
    }
  )
})

test_that("test compare_mic speed", {
  # this is currently slow
  n <- 10
  many_mics <- sample(mic_range(), n, replace = TRUE)
  ab <- sample(c("AMX", "CIP"), n, replace = TRUE)
  mo <- sample(c("Escherichia coli", "Proteus mirabilis"), n, replace = TRUE)
  suppressMessages(
    {
      expect_s3_class(compare_mic(many_mics, many_mics, ab, mo), "mic_validation")
    }
  )

  # this is fast
  n <- 1000
  many_mics <- sample(mic_range(), n, replace = TRUE)
  ab <- rep("AMX", n)
  mo <- sample(c("Escherichia coli", "Proteus mirabilis"), n, replace = TRUE)
  suppressMessages(
    {
      expect_s3_class(compare_mic(many_mics, many_mics, ab, mo), "mic_validation")
    }
  )
})

test_that("test additional arguments passed from compare_mic to AMR::as.sir", {
  gs <- c("0.5", "4", ">8", "2")
  test <- c("0.5", "8", "2", "0.5")
  ab <- "amoxicillin"
  mo <- "Escherichia coli"
  val <- suppressMessages(compare_mic(gs, test, ab, mo, uti = TRUE))
  expect_s3_class(val, "mic_validation")
})

test_that("test subset mic_validation", {
  gs <- c("0.5", "4", ">8", "2")
  test <- c("0.5", "8", "2", "0.5")
  ab <- c("amoxicillin", "amoxicillin", "gentamicin", "gentamicin")
  mo <- c("Escherichia coli", "Proteus mirabilis", "Proteus mirabilis", "Proteus mirabilis")
  val <- suppressMessages(compare_mic(gs, test, ab, mo))
  expect_s3_class(val, "mic_validation")

  # subset amox
  val_sub <- subset(val, ab == "amoxicillin")
  expect_s3_class(val_sub, "mic_validation")
  expect_equal(nrow(val_sub), 2)

  # subset amox and E. coli
  val_sub <- subset(val, ab == "amoxicillin" & mo == "Escherichia coli")
  expect_s3_class(val_sub, "mic_validation")
  expect_equal(nrow(val_sub), 1)
})

test_that("test compare_sir", {
  gs <- c("4", "16", ">8", "0.5")
  test <- c("0.5", "8", "0.25", ">64")
  ab <- "GEN"
  mo <- "Escherichia coli"
  val <- suppressMessages(compare_mic(gs, test, ab, mo))
  expect_s3_class(val, "mic_validation")
  sum_val <- summary(val)
  expect_equal(sum_val$EA_pcent[sum_val$ab == "GEN"], 0.25)
  expect_equal(sum_val$very_major_error_pcent[sum_val$ab == "GEN"], 50)
})

test_that("test fall back to ECOFFS in compare_mic", {
  # single
  val <- suppressMessages(
    compare_mic(AMR::as.mic(0.5),
                     AMR::as.mic(0.5),
                     "CHL",
                     "Escherichia coli",
                     accept_ecoff = TRUE)
  )
  expect_equal(val$error[val$ab == "CHL"], factor(NA,
                                                  levels = c("M", "vM", "m")))

  val <- suppressMessages(
    compare_mic(AMR::as.mic(64),
                     AMR::as.mic(0.5),
                     "CHL",
                     "Escherichia coli",
                     accept_ecoff = TRUE)
  )
  expect_equal(val$error[val$ab == "CHL"], factor("vM",
                                                  levels = c("M", "vM", "m")))

  # multiple
  gs <- c("0.5", "4", ">64", "2")
  test <- c("0.5", "8", "0.25", ">64")
  ab <- "CHL"
  mo <- "Escherichia coli"
  val <- suppressMessages(
    compare_mic(gs, test, ab, mo, accept_ecoff = TRUE))
  expect_s3_class(val, "mic_validation")
  sum_val <- summary(val)
  expect_equal(sum_val$EA_pcent[sum_val$ab == "CHL"], 0.5)
  expect_equal(sum_val$very_major_error_pcent[sum_val$ab == "CHL"], 25)
  expect_equal(sum_val$minor_error_pcent[sum_val$ab == "CHL"], 0)
  expect_equal(sum_val$major_error_pcent[sum_val$ab == "CHL"], 25)

  # here we test some results that have breakpoints and some that do not
  gs <- c("0.5", "4", ">64", "2")
  test <- c("0.5", "8", "0.25", ">64")
  ab <- c("CIP", "CIP", "CHL", "CHL")
  mo <- "Escherichia coli"
  val <- suppressMessages(
    compare_mic(gs, test, ab, mo, accept_ecoff = TRUE))
  expect_s3_class(val, "mic_validation")

  # all mics have breakpoints but accept_ecoff is TRUE
  ab <- "CIP"
  val <- suppressMessages(
    compare_mic(gs, test, ab, mo, accept_ecoff = TRUE))
  expect_s3_class(val, "mic_validation")

  # all mics have breakpoints but accept_ecoff is FALSE (different ab)
  ab <- c("CIP", "CIP", "AMX", "AMX")
  val <- suppressMessages(
    compare_mic(gs, test, ab, mo, accept_ecoff = TRUE))
  expect_s3_class(val, "mic_validation")
})

test_that("test droplevels mic_validation", {
  t <- AMR::as.mic(c("<=0.25", "0.25", "0.5", "1", "2", "1", "0.5"))
  g <- AMR::as.mic(c("0.004", "0.08", "<=0.25", "0.5", "1", "0.5", "0.5"))

  expect_g <- AMR::as.mic(
    c("0.004", "0.06", "<=0.25", "0.5", "1", "0.5", "0.5"))

  v <- compare_mic(g, t)
  expect_equal(droplevels.mic_validation(v)$gold_standard,
               expect_g)

  expect_g <- AMR::as.mic(
    c("<=0.25", "0.25", "0.5", "1", ">1", "1", "0.5"))
  expect_t <- AMR::as.mic(
    c("<=0.25", "0.06", "<=0.25", "0.5", "1", "0.5", "0.5"))
  # same but flip MICs
  v <- compare_mic(t, g)
  expect_equal(droplevels.mic_validation(v)$gold_standard,
               expect_g)
  expect_equal(droplevels.mic_validation(v)$test,
               expect_t)

  # some more tests
  v <- compare_mic(c("0.5", "0.5", "0.5"), c("0.5", "0.5", "0.5"))
  expect_equal(droplevels.mic_validation(v)$gold_standard,
               AMR::as.mic(c("0.5", "0.5", "0.5")))

  t <- AMR::as.mic(c("0.5", "4", "16", "256", "256"))
  g <- AMR::as.mic(c("0.5", "4", ">4", ">4", "2"))

  expect_t <- AMR::as.mic(
    c("0.5", "4", ">4", ">4", "256")
  )

  v <- compare_mic(g, t)
  v_dropped <- droplevels.mic_validation(v)
  expect_equal(v_dropped$test,
               expect_t)
  expect_equal(v_dropped$gold_standard,
               g)
  expect_s3_class(v_dropped$gold_standard, "mic")
  expect_s3_class(v_dropped$test, "mic")
})
