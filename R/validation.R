#' Essential agreement for MIC validation
#'
#' @description
#' Essential agreement calculation for comparing two MIC vectors.
#'
#' @param x AMR::mic or coercible
#' @param y AMR::mic or coercible
#' @param coerce_mic convert to AMR::mic
#' @param tolerate_censoring "strict", "x", "y", or "both" - whether to tolerate
#' censoring in x, y, or both. See details.
#' @param tolerate_matched_censoring "strict", "x", "y", or "both" - how to handle
#' situations where one of the values is censored, but both values match (e.g.,
#' x = ">2", y = "2"). For most situations, this is considered essential agreement.
#' so should be left as "both".
#' @param tolerate_leq whether to tolerate <= in essential agreement, e.g., <=2
#' and 4 will be considered in essential agreement (because <=2 includes 2mg/L,
#' which is within 1 dilution of 4mg/L). This argument respects the
#' tolerate_censoring argument, so if tolerate_censoring is "strict", this will
#' not be applied.
#' @param tolerate_geq whether to tolerate >= in essential agreement, e.g., >=4
#' and 2 will be considered in essential agreement (because >=4 includes 4mg/L,
#' which is within 1 dilution of 2mg/L). This argument respects the
#' tolerate_censoring argument, so if tolerate_censoring is strict, this will
#' not be applied.
#' @param mode "categorical" or "numeric", see details
#' @return logical vector
#' @export
#'
#' @details
#' Essential agreement is a central concept in the comparison of two sets of MIC
#' values. It is most often used when validating a new method against a gold
#' standard. This function reliably performs essential agreement in line with
#' ISO 20776-2:2021. The function can be used in two modes: categorical and
#' numeric. In categorical mode, the function will use traditional MIC
#' concentrations to determine the MIC (therefore it will use force_mic() to
#' convert both x and y to a clean MIC -- see [force_mic]). In numeric mode,
#' the function will compare the ratio of the two MICs, after removing censoring
#' (values that are ">" and "<" are multiplied and divided by 2, respectively ---
#' see [mic_uncensor]).
#' In most cases, categorical mode provides more reliable results.
#' Values within +/- 1 dilutions are considered to be in essential agreement.
#'
#' The tolerate_censoring argument controls how the function handles censored
#' data. If set to "strict", the function will return NA for any pair of
#' values that are both censored (and not equal).
#' If set to "x" or "y", the function will allow one of the values to be censored
#' and will compare the uncensored value to the other value.
#' When set to "both", the function will allow one of the values to be censored.
#' If using "both" and both values are censored, the function will attempt to
#' determine essential agreement based on the ratio of the two values, but a
#' warning will be raised.
#'
#' @references
#' International Organization for Standardization. ISO 20776-2:2021
#' Available from: https://www.iso.org/standard/79377.html
#'
#' @examples
#' x <- AMR::as.mic(c("<0.25", "8", "64", ">64"))
#' y <- AMR::as.mic(c("<0.25", "2", "16", "64"))
#' essential_agreement(x, y)
#' # TRUE FALSE FALSE TRUE
#'
#' # examples using tolerate_censoring
#' x <- AMR::as.mic("<4")
#' y <- AMR::as.mic("0.25")
#'
#' essential_agreement(x, y, tolerate_censoring = "x") # TRUE
#' essential_agreement(x, y, tolerate_censoring = "y") # FALSE
#' essential_agreement(x, y, tolerate_censoring = "both") # TRUE (same as "x")
#'
#' # strict returns FALSE as it wants the censoring cut-offs to be close
#' essential_agreement(x, y, tolerate_censoring = "strict")
essential_agreement <- function(x,
                                y,
                                coerce_mic = TRUE,
                                tolerate_censoring = "strict",
                                tolerate_matched_censoring = "both",
                                tolerate_leq = TRUE,
                                tolerate_geq = TRUE,
                                mode = "categorical") {
  if (any(!AMR::is.mic(c(x, y))) & !coerce_mic) {
    stop("Both MIC inputs to essential_agreement must be AMR::mic.
Convert using AMR::as.mic() with or without MIC::force_mic().")
  }

  stopifnot(
    "tolerate_censoring must be one of 'strict', 'x', 'y', or 'both'" =
      tolerate_censoring %in% c("strict", "x", "y", "both"))

  if (mode == "categorical") {
    .xbak <- suppressWarnings(
      AMR::as.mic(x))
    .ybak <- suppressWarnings(
      AMR::as.mic(y))

    if (any(is.na(.xbak)) | any(is.na(.ybak))) {
      stop("Both x and y must be mic-compatible values.")
    }

    x <- mic_uncensor(.xbak)
    y <- mic_uncensor(.ybak)

    range <- mic_range()
    range <- c(max(range) * 2,
               range,
               min(range) / 2)

    if (any(.xbak < min(range),
            .xbak > max(range),
            .ybak < min(range),
            .ybak > max(range))) {
      stop("Input MICs must be above {min(range)} and below {max(range)}")
    }

    index_diff <- logical(length(x))
    for (i in seq_along(x)) {
      xi <- x[i]
      yi <- y[i]
      if (is.na(xi) || is.na(yi)) {
        index_diff[i] <- NA
        next
      }
      if (xi == yi) {
        index_diff[i] <- TRUE
        next
      }
      # mismatched censoring
      if (grepl("<", .xbak[i]) & grepl(">", .ybak[i]) |
          grepl(">", .xbak[i]) & grepl("<", .ybak[i])) {
        index_diff[i] <- FALSE
        next
      }

      # x == y, but x is censored (e.g., x = ">2", "y = "2")
      if (all(grepl("<|>", .xbak[i]),
              !grepl("<|>", .ybak[i]),
              as.numeric(.xbak[i]) == as.numeric(.ybak[i]))) {

        if (tolerate_matched_censoring == "both" |
            tolerate_matched_censoring == "x") {
          index_diff[i] <- TRUE
          next
        } else {
          index_diff[i] <- FALSE
          next
        }
      }

      # y == x, but y is censored (e.g., x = "2", "y = ">2")
      if (all(!grepl("<|>", .xbak[i]),
              grepl("<|>", .ybak[i]),
              as.numeric(.xbak[i]) == as.numeric(.ybak[i]))) {
        if (tolerate_matched_censoring == "both" |
            tolerate_matched_censoring == "y") {
          index_diff[i] <- TRUE
          next
        } else {
          index_diff[i] <- FALSE
          next
        }
      }

      # x is 2 * y, but y is equal censored (e.g., x = "4", y = "<=2")
      # if tolerate_censoring is "both" or "y", then this is essential agreement
      if (all(!grepl("<|>", .xbak[i]),
              grepl("<=", .ybak[i]),
              as.numeric(xi) == 2 * as.numeric(.ybak[i]))) {
        if ((tolerate_censoring == "both" | tolerate_censoring == "y") & tolerate_leq) {
          index_diff[i] <- TRUE
          next
        } else {
          index_diff[i] <- FALSE
          next
        }
      }

      # x is y / 2, but y is equal censored (e.g., x = "2", y = ">=4")
      # if tolerate_censoring is "both" or "x", then this is essential agreement
      if (all(grepl(">=", .ybak[i]),
              !grepl("<|>", .xbak[i]),
              as.numeric(xi) == as.numeric(.ybak[i]) / 2)) {
        if ((tolerate_censoring == "both" | tolerate_censoring == "x") & tolerate_geq) {
          index_diff[i] <- TRUE
          next
        } else {
          index_diff[i] <- FALSE
          next
        }
      }

      # y is 2 * x, but x is equal censored (e.g., x = "<=2", y = "4")
      # if tolerate_censoring is "both" or "x", then this is essential agreement
      if (all(grepl("<=", .xbak[i]),
              !grepl("<|>", .ybak[i]),
              as.numeric(yi) == 2 * as.numeric(.xbak[i]))) {
        if ((tolerate_censoring == "both" | tolerate_censoring == "x") & tolerate_leq) {
          index_diff[i] <- TRUE
          next
        } else {
          index_diff[i] <- FALSE
          next
        }
      }

      # y is x / 2, but x is equal censored (e.g., x = ">=4", y = "2")
      # if tolerate_censoring is "both" or "y", then this is essential agreement
      if (all(!grepl("<|>", .ybak[i]),
              grepl(">=", .xbak[i]),
              as.numeric(yi) == as.numeric(.xbak[i]) / 2)) {
        if ((tolerate_censoring == "both" | tolerate_censoring == "y") & tolerate_geq) {
          index_diff[i] <- TRUE
          next
        } else {
          index_diff[i] <- FALSE
          next
        }
      }

      if (all(grepl("<|>", .xbak[i]),
              grepl("<|>", .ybak[i]),
              .xbak[i] != .ybak[i])) {
        if (tolerate_censoring == "strict") {
          index_diff[i] <- NA
          message(glue::glue("Unable to determine essential agreement for
                             censored values: {.xbak[i]} vs {.ybak[i]}.",
                             "Set tolerate_censoring to 'x', 'y', or 'both' to allow comparison.",
                             .sep = "\n"))
          next
        } else if (tolerate_censoring == "both") {
          warning(glue::glue("Both x and y are censored: {.xbak[i]} vs {.ybak[i]}.",
                             "Comparison may not be reliable.",
                             .sep = "\n"))
          diff <- abs(which(range == xi) - which(range == yi))
          index_diff[i] <- diff <= 1
          next
        }
      }
      if (tolerate_censoring == "x" | tolerate_censoring == "both") {
        if (grepl("<", .xbak[i])) {
          if (as.numeric(.ybak[i]) <= .xbak[i]) {
            index_diff[i] <- TRUE
            next
          } else {
            index_diff[i] <- FALSE
            next
          }
        }

        if (grepl(">", .xbak[i])) {
          if (as.numeric(.ybak[i]) >= .xbak[i]) {
            index_diff[i] <- TRUE
            next
          } else {
            index_diff[i] <- FALSE
            next
          }
        }
      }

      if (tolerate_censoring == "y" | tolerate_censoring == "both") {
        if (grepl("<", .ybak[i])) {
          if (as.numeric(.xbak[i]) <= .ybak[i]) {
            index_diff[i] <- TRUE
            next
          } else {
            index_diff[i] <- FALSE
            next
          }
        }

        if (grepl(">", .ybak[i])) {
          if (as.numeric(.xbak[i]) >= .ybak[i]) {
            index_diff[i] <- TRUE
            next
          } else {
            index_diff[i] <- FALSE
            next
          }
        }
      }
      diff <- abs(which(range == xi) - which(range == yi))
      index_diff[i] <- diff <= 1
    }
    return(index_diff)

  }
  if (mode == "numerical") {
    x <- as.numeric(mic_uncensor(x))
    y <- as.numeric(mic_uncensor(y))

    frac <- x / y
    return(as.logical(
      dplyr::case_when(
        is.na(frac) ~ NA,
        frac > 2.0 ~ FALSE,
        frac < 0.5 ~ FALSE,
        frac != 1.0 ~ TRUE,
        TRUE ~ TRUE)))
  }
  stop("Mode must be categorical or numerical")
}

#' Compare and validate MIC values
#'
#' @param gold_standard vector of MICs to compare against.
#' @param test vector of MICs that are under investigation
#' @param ab character vector (same length as MIC) of antibiotic names (optional)
#' @param mo character vector (same length as MIC) of microorganism names (optional)
#' @param accept_ecoff if TRUE, ECOFFs will be used when no clinical breakpoints are available
#' @param simplify if TRUE, MIC values will be coerced into the closest halving
#' dilution (e.g., 0.55 will be converted to 0.5)
#' @param ea_mode "categorical" or "numeric", see [essential_agreement]
#' @param tolerate_censoring "strict", "gold_standard", "test", or "both" - how to handle
#' censored data (see [essential_agreement] for details). Generally, this should be
#' left as "gold_standard" since this setting "tolerates" a test that has higher
#' granularity (i.e., less censoring) than the gold standard. Setting to "test"
#' or "both" should be used with caution but may be appropriate in some cases
#' where the test also produces censored results.
#' @param tolerate_matched_censoring "strict", "gold_standard", "test",
#' or "both" - how to handle situations where one of the values is censored,
#' but both values match (e.g., gold_standard = ">2", test = "2"). Generally, this
#' should be left as "both", since these values are considered to be in
#' essential agreement. For more details, see [essential_agreement].
#' @param tolerate_leq whether to tolerate <= in essential agreement, e.g., <=2
#' and 4 will be considered in essential agreement. See [essential_agreement]
#' for details.
#' @param tolerate_geq whether to tolerate >= in essential agreement, e.g.,  >=4
#' and 2 will be considered in essential agreement. See [essential_agreement]
#' for details.
#' @param ... additional arguments to be passed to AMR::as.sir
#'
#' @return S3 mic_validation object
#'
#' @description
#' This function compares an vector of MIC values to another. Generally, this is
#' in the context of a validation experiment -- an investigational assay or
#' method (the "test") is compared to a gold standard. The rules used by this
#' function are in line with "ISO 20776-2:2021 Part 2: Evaluation of performance
#' of antimicrobial susceptibility test devices against reference broth
#' micro-dilution."
#'
#' There are two levels of detail that are provided. If only the MIC values are
#' provided, the function will look for essential agreement between the two sets
#' of MIC. If the organism and antibiotic arguments are provided, the function
#' will also calculate the categorical agreement using EUCAST breakpoints (or,
#' if breakpoint not available and accept_ecoff = TRUE, ECOFFs).
#'
#' The function returns a special dataframe of results, which is also an
#' mic_validation object. This object can be summarised using summary() for
#' summary metrics, plotted using plot() for an essential agreement confusion
#' matrix, and tabulated using table().
#'
#' @export
#'
#' @examples
#' # Just using MIC values only
#' gold_standard <- c("<0.25", "8", "64", ">64")
#' test <- c("<0.25", "2", "16", "64")
#' val <- compare_mic(gold_standard, test)
#' summary(val)
#'
#' # Using MIC values and antibiotic and organism names
#' gold_standard <- c("<0.25", "8", "64", ">64")
#' test <- c("<0.25", "2", "16", "64")
#' ab <- c("AMK", "AMK", "AMK", "AMK")
#' mo <- c("B_ESCHR_COLI", "B_ESCHR_COLI", "B_ESCHR_COLI", "B_ESCHR_COLI")
#' val <- compare_mic(gold_standard, test, ab, mo)
#' "error" %in% names(val)  # val now has categorical agreement
compare_mic <- function(gold_standard,
                        test,
                        ab = NULL,
                        mo = NULL,
                        accept_ecoff = FALSE,
                        simplify = TRUE,
                        ea_mode = "categorical",
                        tolerate_censoring = "gold_standard",
                        tolerate_matched_censoring = "both",
                        tolerate_leq = TRUE,
                        tolerate_geq = TRUE,
                        ...) {
  if (length(gold_standard) != length(test)) {
    stop("Gold standard and test must be the same length")
  }

  call_list <- list(accept_ecoff = accept_ecoff,
                    simplify = simplify,
                    ea_mode = ea_mode,
                    tolerate_censoring = tolerate_censoring,
                    tolerate_matched_censoring = tolerate_matched_censoring,
                    tolerate_leq = tolerate_leq,
                    tolerate_geq = tolerate_geq)

  gold_standard_mod <- gold_standard |>
    force_mic(levels_from_AMR = !simplify) |>
    AMR::as.mic()

  test_mod <- test |>
    force_mic(levels_from_AMR = !simplify) |>
    AMR::as.mic()

  if (tolerate_censoring == "gold_standard") {
    tolerate_censoring <- "x"
  } else if (tolerate_censoring == "test") {
    tolerate_censoring <- "y"
  }
  output <- list(
    gold_standard = gold_standard_mod,
    test = test_mod,
    essential_agreement = factor(essential_agreement(gold_standard_mod,
                                                     test_mod,
                                                     mode = ea_mode,
                                                     tolerate_censoring = tolerate_censoring,
                                                     tolerate_matched_censoring = tolerate_matched_censoring,
                                                     tolerate_leq = tolerate_leq,
                                                     tolerate_geq = tolerate_geq),
                                 levels = c(FALSE, TRUE))
  )

  if (!is.null(ab)) {
    if (length(ab) == 1) {
      ab <- rep(ab, length(gold_standard))
    }
    if (length(ab) > 1 & length(ab) != length(gold_standard)) {
      stop("Antibiotic names must be the same length as MIC values, or single value")
    }
    output[["ab"]] <- ab
  }

  if (!is.null(mo)) {
    if (length(mo) == 1) {
      mo <- rep(mo, length(gold_standard))
    }
    if (length(mo) > 1 & length(mo) != length(gold_standard)) {
      stop("Microorganism names must be the same length as MIC values, or single value")
    }
    output[["mo"]] <- mo
  }

  if (!is.null(ab) & !is.null(mo)) {
    gold_standard_sir <- as.sir_vectorised(gold_standard_mod, mo, ab, accept_ecoff, ...)
    test_sir <- as.sir_vectorised(test_mod, mo, ab, accept_ecoff, ...)
    output[["gold_standard_sir"]] <- gold_standard_sir
    output[["test_sir"]] <- test_sir
    output[["error"]] <- compare_sir(gold_standard_sir,
                                   test_sir)
  }

  validation_class <- "single_ab_validation"
  if (!is.null(ab) && length(unique(ab)) > 1) {
    validation_class <- "multi_ab_validation"
  }

  class(output) <- append(class(output), c(validation_class, "mic_validation"), 0)
  attr(output, "call") <- call_list
  output
}

drop_levels_mic_validation <- function(x, target, source,
                                       lower = TRUE,
                                       safe = TRUE) {
  than <- ifelse(lower, "<=", ">")
  than_fun <- ifelse(lower, `<`, `>`)
  bound_fun <- ifelse(lower, min, max)

  source_bound <- bound_fun(x[[source]], na.rm = TRUE)
  source_bound <- force_mic(source_bound)

  indices_to_change <- than_fun(x[[target]], source_bound) & x$essential_agreement == TRUE
  x[[target]][indices_to_change] <-
    AMR::as.mic(paste0(than, as.numeric(source_bound)))


  x
}

#' Droplevels for MIC validation object
#'
#' @param x mic_validation object
#' @param safe ensure that essential agreement is not changed after dropping
#' levels
#' @param ... additional arguments
#'
#' @return mic_validation object
#' @export
#'
#' @description
#' Quite often, MIC values are being compared across methods with different
#' levels of granularity. For example, the true MIC may be measured across a
#' higher range of values than the test method. This means that there may be
#' MIC levels that don't provide much additional information (since they are
#' only present in one of the methods). This function removes these unnecessary
#' levels at both ranges of the MIC values.
#'
#' This function ensure that the changes do not "change" the essential
#' agreement interpretation. This can be suppressed using safe = FALSE,
#' however this is probably not desired behaviour.
#'
#' @examples
#' gold_standard <- c("<0.25", "0.25", "0.5", "1", "2", "1", "0.5")
#' test <- c("0.004", "0.08", "<0.25", "0.5", "1", "0.5", "0.5")
#' val <- compare_mic(gold_standard, test)
#' droplevels(val)
droplevels.mic_validation <- function(x,
                                      safe = TRUE,
                                      ...) {
  x <- drop_levels_mic_validation(x, target = "gold_standard",
                                  source = "test",
                                  lower = TRUE,
                                  safe = safe)
  x <- drop_levels_mic_validation(x, target = "test",
                                  source = "gold_standard",
                                  lower = TRUE,
                                  safe = safe)

  x <- drop_levels_mic_validation(x, target = "gold_standard",
                                  source = "test",
                                  lower = FALSE,
                                  safe = safe)
  x <- drop_levels_mic_validation(x, target = "test",
                                  source = "gold_standard",
                                  lower = FALSE,
                                  safe = safe)

  if (safe) {
    tol_censor_call <- attr(x, "call")$tolerate_censoring
    if (tol_censor_call == "gold_standard") {
      tol_censor_call <- "x"
    } else if (tol_censor_call == "test") {
      tol_censor_call <- "y"
    }

    tol_censor_match_call <- attr(x, "call")$tolerate_matched_censoring
    if (tol_censor_match_call == "gold_standard") {
      tol_censor_match_call <- "x"
    } else if (tol_censor_match_call == "test") {
      tol_censor_match_call <- "y"
    }

    new_ea <- essential_agreement(
      x[["gold_standard"]],
      x[["test"]],
      coerce_mic = FALSE,
      mode = "categorical",
      tolerate_censoring = tol_censor_call,
      tolerate_matched_censoring = tol_censor_match_call
    )

    if (!all(new_ea == x$essential_agreement)) {
      stop(
        glue::glue(
          "Essential agreement does not match after dropping levels.
          You can ignore and force the levels to be dropped using safe = FALSE"))
    }
  }
  x

}

#' Print MIC validation object
#'
#' @param x mic_validation object
#' @param ... additional arguments
#'
#' @return character
#' @export
#'
#' @examples
#' gold_standard <- c("<0.25", "8", "64", ">64")
#' test <- c("<0.25", "2", "16", "64")
#' val <- compare_mic(gold_standard, test)
#' print(val)
print.mic_validation <- function(x, ...) {
  if (!is.null(x$ab) & !is.null(x$mo)) {
    antibiotic_str <- paste(unique(x$ab), collapse = ", ")
    organism_str <- paste(unique(x$mo), collapse = ", ")
    return(
      cat(glue::glue("MIC validation object with {length(x$gold_standard)} observations
                     Agreement type: essential and categorical
                     Antibiotics: {antibiotic_str}
                     Organisms: {organism_str}"))
      )
    }

  return(
    cat(glue::glue("MIC validation object with {length(x$gold_standard)} observations
                   Agreement type: essential"))
  )

}

#' Subset MIC validation object
#'
#' @param x mic_validation object
#' @param subset logical expression to subset by
#' @param ... additional arguments
#'
#' @return mic_validation object
#' @export
#'
#' @examples
#' gold_standard <- c("<0.25", "8", "64", ">64")
#' test <- c("<0.25", "2", "16", "64")
#' ab <- AMR::as.ab(c("AMK", "AMK", "CIP", "CIP"))
#' mo <- AMR::as.mo(c("E. coli", "E. coli", "P. mirabilis", "P. mirabilis"))
#' val <- compare_mic(gold_standard, test, ab, mo)
#' subset(val, ab == AMR::as.ab("AMX"))
#' subset(val, mo == AMR::as.mo("E. coli"))
subset.mic_validation <- function(x, subset, ...) {
  filtered <- x |>
    as.data.frame() |>
    subset(subset, ...)
  class(filtered) <- append(class(filtered), "mic_validation", 0)
  filtered
}

#' Summary of MIC validation results
#'
#' @param object S3 mic_validation object
#' @param ... further optional parameters
#'
#' @export
#'
#' @return S3 mic_validation_summary object
#'
#' @description
#' Summarise the results of an MIC validation generated using compare_mic().
#'
#' @examples
#' gold_standard <- c("<0.25", "8", "64", ">64")
#' test <- c("<0.25", "2", "16", "64")
#' val <- compare_mic(gold_standard, test)
#' summary(val)
#' # or, for more detailed results
#' as.data.frame(summary(val))
summary.mic_validation <- function(object,
                                   ...) {
  if (!"ab" %in% names(object) & !"mo" %in% names(object)) {
    output <- list(EA_n = sum(object$essential_agreement == TRUE),
                   EA_pcent = sum(object$essential_agreement == TRUE) / length(object$essential_agreement),
                   bias = bias(object$gold_standard, object$test),
                   n = length(object$essential_agreement))
    class(output) <- append(class(output), "mic_validation_summary", after = 0)
    return(output)
  }

  if ("ab" %in% names(object) & !"mo" %in% names(object)) {
      output <- object |>
        as.data.frame() |>
        dplyr::group_by(.data[["ab"]]) |>
        dplyr::summarise(EA_n = sum(.data[["essential_agreement"]] == TRUE),
                         EA_pcent = sum(.data[["essential_agreement"]] == TRUE) / length(.data[["essential_agreement"]]),
                         bias = bias(object$gold_standard, object$test),
                         n = dplyr::n())
      class(output) <- append(class(output), "mic_validation_summary", after = 0)
      return(output)
  }

  if ("ab" %in% names(object) & "mo" %in% names(object)) {
    output <- object |>
        as.data.frame() |>
        dplyr::group_by(.data[["ab"]], .data[["mo"]]) |>
        dplyr::summarise(
          EA_pcent = sum(.data[["essential_agreement"]] == TRUE) / length(.data[["essential_agreement"]]),
          bias = bias(.data[["gold_standard"]], .data[["test"]]),
          resistant_pcent = AMR::proportion_R(.data[["gold_standard_sir"]], minimum = 1, as_percent = FALSE) * 100,
          minor_error_pcent = sum(.data[["error"]] == "m", na.rm = TRUE) / length(.data[["error"]]) * 100,
          major_error_pcent = sum(.data[["error"]] == "M", na.rm = TRUE) / length(.data[["error"]]) * 100,
          very_major_error_pcent = sum(.data[["error"]] == "vM", na.rm = TRUE) /length(.data[["error"]]) * 100,
          EA_n = sum(.data[["essential_agreement"]] == TRUE),
          resistant_n = AMR::count_resistant(.data[["gold_standard_sir"]]),
          minor_error_n = sum(.data[["error"]] == "m", na.rm = TRUE),
          major_error_n = sum(.data[["error"]] == "M", na.rm = TRUE),
          very_major_error_n = sum(.data[["error"]] == "vM", na.rm = TRUE),
          n = dplyr::n())
    class(output) <- append(class(output), "mic_validation_summary", after = 0)
    return(output)

  }
}

#' Print MIC validation summary
#'
#' @param x mic_validation_summary object
#' @param ... additional arguments
#'
#' @return character
#' @export
#'
#' @examples
#' gold_standard <- c("<0.25", "8", "64", ">64")
#' test <- c("<0.25", "2", "16", "64")
#' val <- compare_mic(gold_standard, test)
#' print(summary(val))
print.mic_validation_summary <- function(x, ...) {
  if (!"ab" %in% names(x) & !"mo" %in% names(x)) {
    return(
      cat(glue::glue("MIC validation summary
                     Essential agreement: {x$EA_n} ({round(x$EA_pcent * 100, 2)}%)
                     Bias: {x$bias}"))
    )
  }
  if ("ab" %in% names(x) & !"mo" %in% names(x)) {
    antibiotic_str <- paste(unique(x$ab), collapse = ", ")
    return(
      cat(glue::glue("MIC validation summary
                     Number of observations: {sum(x$n)}
                     Antibiotic: {antibiotic_str}
                     Essential agreement: {sum(x$EA_n)} ({round(sum(x$EA_n) / sum(x$n) * 100, 2)}%)
                     Mean bias: {mean(x$bias)}
                     *Use as.data.frame() to see full summary*"))
    )
  }
  if (length(x) > 1 & "ab" %in% names(x) & "mo" %in% names(x)) {
    antibiotic_str <- paste(unique(x$ab), collapse = ", ")
    organism_str <- paste(unique(x$mo), collapse = ", ")

    return(
      cat(glue::glue("MIC validation summary
                     Antibiotic: {antibiotic_str}
                     Organism: {organism_str}
                     Essential agreement: {sum(x$EA_n)} ({round(sum(x$EA_n) / sum(x$n) * 100, 2)}%)
                     Resistant: {sum(x$resistant_n)} ({round(sum(x$resistant_n) / sum(x$n) * 100, 2)}%)
                     Minor errors: {sum(x$minor_error_n)} ({round(sum(x$minor_error_n) / sum(x$n) * 100, 2)}%)
                     Major errors: {sum(x$major_error_n)} ({round(sum(x$major_error_n) / sum(x$n) * 100, 2)}%)
                     Very major errors: {sum(x$very_major_error_n)} ({round(sum(x$very_major_error_n) / sum(x$n) * 100, 2)}%)
                     Mean bias: {mean(x$bias)}
                     N: {sum(x$n)}
                     *Use as.data.frame() to see full summary*"))
    )
  }
}

#' Compare SIR results and generate categorical agreement
#'
#' @param gold_standard Susceptibility results in AMR::sir format
#' @param test Susceptibility results in AMR::sir format
#'
#' @return factor vector with the following levels: M, vM, m.
#' @export
#'
#' @description
#' Compare two AMR::sir vectors and generate a categorical agreement vector with
#' the following levels: M (major error), vM (very major error), m (minor error).
#' The error definitions are:
#'
#' 1. Major error (M): The test result is resistant (R) when the gold standard
#' is susceptible (S).
#' 2. vM (very major error): The test result is susceptible (S) when the gold
#' standard is resistant (R).
#' 3. Minor error (m): The test result is intermediate (I) when the gold standard
#' is susceptible (S) or resistant (R), or vice versa.
#'
#' @examples
#' gold_standard <- c("S", "R", "I", "I")
#' gold_standard <- AMR::as.sir(gold_standard)
#' test <- c("S", "I", "R", "R")
#' test <- AMR::as.sir(test)
#' compare_sir(gold_standard, test)
compare_sir <- function(gold_standard, test) {
  if (length(gold_standard) != length(test)) {
    stop("Inputs to compare_sir must be same length")
  }
  if (!all(AMR::is.sir(gold_standard)) | !all(AMR::is.sir(test))) {
    stop("Inputs to compare_sir must be AMR::sir objects")
  }
  return(
    factor(
      dplyr::case_when(
        gold_standard == "S" & test == "I" ~ "m",
        gold_standard == "R" & test == "I" ~ "m",
        gold_standard == "I" & test == "S" ~ "m",
        gold_standard == "I" & test == "R" ~ "m",
        gold_standard == "S" & test == "R" ~ "M",
        gold_standard == "R" & test == "S" ~ "vM",
        TRUE ~ NA
      ),
      levels = c("M", "vM", "m")
    )
  )
}

#' Calculate MIC bias
#'
#' @param gold_standard AMR::mic vector
#' @param test AMR::mic vector
#'
#' @return numeric value
#' @export
#'
#' @description
#' Calculate the bias between two AMR::mic vectors. The bias is calculated as
#' the percentage of test MICs that are above the gold standard MICs minus the
#' percentage of test MICs that are below the gold standard MICs.
#'
#' @references
#' International Organization for Standardization. ISO 20776-2:2021
#' Available from: https://www.iso.org/standard/79377.html
#'
#' @examples
#' gold_standard <- c("<0.25", "8", "64", ">64")
#' test <- c("<0.25", "2", "16", "64")
#' bias(gold_standard, test)
bias <- function(gold_standard, test) {
  if (length(gold_standard) != length(test)) {
    stop("Inputs to bias must be same length")
  }
  bias_above <- sum(test > gold_standard)
  bias_below <- sum(test < gold_standard)

  n <- length(gold_standard)
  return(bias_above / n * 100 - bias_below / n * 100)
}

#' Table
#'
#' @param x an mic_validation object to be tabulated into an essential agreement
#' frequency table, or object/s to be passed to base::table
#' @param ... further arguments
#'
#' @rdname table
#' @export
#'
table <- function(x, ...) {
  UseMethod("table")
}

#' @rdname table
#' @export
table.default <- function(x, ...) {
  return(base::table(x, ...))
}

tabulate_flex <- function(t, ea, bold, ea_color, gold_standard_name, test_name) {
  t_flex <- t |>
    stats::addmargins(FUN = list(Total = sum),
                      quiet = TRUE) |>
    as.data.frame.matrix() |>
    tibble::rownames_to_column(var = "test_mics") |>
    dplyr::mutate("test" = test_name, .before = "test_mics") |>
    flextable::flextable()

  for (i in 1:nrow(t)) {
    for (j in 1:ncol(t)) {
      if (ea[i, j] == TRUE) {
        if (bold) {
          t_flex <- flextable::bold(t_flex, i = i, j = j + 2)
        }
        if (!is.null(ea_color)) {
          t_flex <- flextable::bg(t_flex, i = i, j = j + 2, bg = ea_color)
        }
      }
    }
  }
  t_flex <- t_flex |>
    flextable::set_header_labels("test_mics" = "", test = "") |>
    flextable::add_header_row(values = c("", gold_standard_name, ""),
                              colwidths = c(1, length(t_flex[["col_keys"]]) - 2, 1)) |>
    flextable::align(align = "center", i=1,j=2,part = "header") |>
    flextable::align(align = "right", j=2) |>
    flextable::merge_v(j = 1) |>
    flextable::border_remove() |>
    flextable::align(align = "right", part = "body")
  t_flex
}

#' @rdname table
#'
#' @param x mic_validation S3 object
#' @param format simple or flextable
#' @param fill_dilutions Fill dilutions that are not present in the data in
#' order to match the y- and x- axes
#' @param bold Bold cells where essential agreement is TRUE
#' @param ea_color Background color for essential agreement cells
#' @param gold_standard_name Name of the gold standard to display in output
#' @param test_name Name of the test to display in output
#' @param ... further arguments
#'
#' @return table or flextable object
#'
#' @export
#'
#' @examples
#' gold_standard <- c("<0.25", "8", "64", ">64")
#' test <- c("<0.25", "2", "16", "64")
#' val <- compare_mic(gold_standard, test)
#' table(val)
table.mic_validation <- function(x,
                                 format = 'flextable',
                                 fill_dilutions = TRUE,
                                 bold = TRUE,
                                 ea_color = NULL,
                                 gold_standard_name = "Gold Standard",
                                 test_name = "Test",
                                 ...) {
  test <- match_levels(x[["test"]], match_to = x[["gold_standard"]])
  gold_standard <- match_levels(x[["gold_standard"]], match_to = x[["test"]])

  if (fill_dilutions) {
    test <- fill_dilution_levels(test)
    gold_standard <- fill_dilution_levels(gold_standard)
  }

  t <- table(test,
             gold_standard)

  rnames <- rownames(t)
  cnames <- colnames(t)
  grid <- expand.grid(rnames, cnames)
  ea <- essential_agreement(as.character(grid$Var1), as.character(grid$Var2))
  ea <- matrix(ea, nrow = length(rnames))
  ea <- as.data.frame(ea)

  if (format == "simple") {
    return(t)
  }
  if (format == "flextable") {
    return(tabulate_flex(t, ea, bold, ea_color, gold_standard_name, test_name))
  }

  stop("Format must be simple or flextable")
}
