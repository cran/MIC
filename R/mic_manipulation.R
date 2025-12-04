#' Get MIC meta-data from feature database
#'
#' @param x dataframe containing meta-data
#' @param ids vector of IDs to get meta-data for
#' @param ab_col column name containing MIC results
#' @param id_col column name containing IDs
#' @param as_mic return as AMR::as.mic
#' @param prefer_high_mic where multiple MIC results per ID, prefer the higher MIC
#' @param simplify return as vector of MICs (vs dataframe)
#'
#' @return vector containing MICs, or dataframe of IDs and MICs
#' @export
#'
#' @description
#' This function helps extract MICs from a database of results. It is compatible
#' with the PATRIC meta data format when used on a tidy_patric_db object,
#' created using tidy_patric_db().
#'
#' If more than one MIC is present for a particular observation, the function
#' can return the higher MIC by setting prefer_high_mic = TRUE. If
#' prefer_high_mic = FALSE, the lower MIC will be returned.
#'
#' @examples
#' df <- data.frame(genome_id = c("a_12", "b_42", "x_21", "x_21", "r_75"),
#'                  gentamicin = c(0.25, 0.125, 32.0, 16.0, "<0.0125"))
#' get_mic(df,
#'         ids = c("b_42", "x_21"),
#'         ab_col = "gentamicin",
#'         id_col = "genome_id",
#'         as_mic = FALSE,
#'         prefer_high_mic = TRUE,
#'         simplify = TRUE)
get_mic <- function(x,
                    ids,
                    ab_col,
                    id_col = NULL,
                    as_mic = TRUE,
                    prefer_high_mic = TRUE,
                    simplify = TRUE) {
  if ("save_order" %in% names(x)) {
    stop("Unable to work with x that contains column name 'save_order', please
         rename.")
  }

  if (inherits(x, "tidy_patric_db")) {
    id_col = "genome_id"
  }

  if (!inherits(x, "tidy_patric_db") & is.null(id_col)) {
    stop("Provide id_col or pre-process meta data using tidy_patric_db()")
  }

  x <- x[order(x[[id_col]], x[[ab_col]], decreasing = prefer_high_mic),]
  x <- x[!duplicated(x[[id_col]]),]

  output <- data.frame(ids)
  names(output) <- id_col
  output$save_order <- 1:nrow(output)
  output <- merge(output, x,
                  by = id_col,
                  sort = FALSE,
                  all.x = TRUE)
  if (as_mic) {
    rlang::check_installed("AMR", "To return as MIC class, AMR package must be installed")
    output[[ab_col]] <- AMR::as.mic(output[[ab_col]])
  }
  output <- output[order(output$save_order),]
  output$save_order <- NULL
  if (simplify) {
    return(output[[ab_col]])
  } else {
    return(output)
  }
}

#' Clean up raw MIC for use as a feature
#'
#' @param mic character containing MIC/s
#'
#' @return character of clean MIC/s
#' @description
#' Removes leading "=" which can sometimes be present in raw MIC results. Also converts co-trimoxazole to trimethprim component only.
#'
#' @export
#'
#' @examples
#' clean_raw_mic(c("==>64","0.25/8.0"))
clean_raw_mic <- function(mic) {
  mic <- stringr::str_remove(mic, pattern = "^=+")
  stringr::str_remove(mic, "/(.)+")
}

#' Uncensor MICs
#'
#'
#' @param mic vector of MICs to uncensor; will be coerced to MIC using AMR::as.mic
#' @param method method to uncensor MICs (scale, simple, or bootstrap)
#' @param scale scalar to multiply or divide MIC by (for method = scale)
#' @param ab antibiotic name (for method = bootstrap)
#' @param mo microorganism name (for method = bootstrap)
#' @param distros dataframe of epidemiological distributions (only used,
#' optionally, for method = bootstrap)
#'
#' @return vector of MICs in AMR::mic format
#'
#' @details
#' Censored MIC data is generally unsuitable for modelling without some
#' conversion of censored data. The default behaviour (method = scale) is to
#' halve MICs under the limit of detection (<=) and double MICs above the limit
#' of detection (>). When used with method = simple, this function effectively
#' just removes the censoring symbols, e.g., <=2 becomes 2, and >64 becomes 64.
#'
#' The bootstrap method is the more complex of the three available methods. It
#' attempts to use a second (uncensored) MIC distribution to sample values in
#' the censored range. These values are then used to populate and uncensor
#' the MIC data provided as input (mic). The second (uncensored) MIC
#' distribution is ideally provided from similar experimental conditions.
#' Alternatively, epidemiological distributions can be used. These distributions
#' should be provided as a dataframe to the distros argument. The format for
#' this dataframe is inspired by the EUCAST epidemiological distributions, see:
#' https://www.eucast.org/bacteria/mic-and-zone-distributions-ecoffs/. The dataframe
#' should contain columns for antimicrobial (converted using AMR::as.ab),
#' organism (converted using AMR::as.mo), and MIC concentrations. An example
#' is provided in the 'ecoffs' dataset available with this pacakge. Currently,
#' only Escherichia coli is available in this dataset. Each observation (row)
#' consists of the frequency a particular MIC concentration is observed in the
#' distribution. If such a dataframe is not provided to distros, the function
#' will attempt to use 'ecoffs', but remains limited to E. coli.
#'
#' @references https://www.eucast.org/bacteria/mic-and-zone-distributions-ecoffs/
#'
#' @export
#'
#' @examples
#' mic_uncensor(c(">64.0", "<0.25", "8.0"), method = "scale", scale = 2)
mic_uncensor <- function(mic,
                         method = "scale",
                         scale = 2,
                         ab = NULL,
                         mo = NULL,
                         distros = NULL) {
  if (method == "scale") {
    return(mic_uncensor_scale(mic, scale))
  }
  if (method == "simple") {
    return(mic_uncensor_simple(mic))
  }
  if (method == "bootstrap") {
    if (any(is.null(ab),
            is.null(mo))) {
      stop("For method = bootstrap, both ab and mo must be provided")
    }
    return(AMR::as.mic(mic_uncensor_bootstrap(mic = mic, ab = ab, mo = mo,
                                              distros = distros)))
  }
  stop("Method must be scale, simple or bootstrap")
}

mic_uncensor_simple <- function(mic) {
  AMR::as.mic(as.numeric(AMR::as.mic(mic)))
}

mic_uncensor_scale <- function(mic, scale) {
  suppressWarnings(
    dplyr::case_when(
      stringr::str_detect(mic, ">") ~ AMR::as.mic(AMR::as.mic(stringr::str_remove(force_mic(mic), ">")) * scale),
      stringr::str_detect(mic, "<") ~ AMR::as.mic(AMR::as.mic(stringr::str_remove(force_mic(mic), "<")) / scale),
      .default = AMR::as.mic(force_mic(mic)))
  )
}

mic_uncensor_bootstrap <- Vectorize(function(mic, ab, mo, distros) {
  if (is.na(mic)) {
    return(AMR::NA_mic_)
  }
  mic <- gsub("=", "", mic)
  if (!startsWith(mic, "<") & !startsWith(mic, ">")) {
    return(AMR::as.mic(mic))
  }

  ab <- AMR::as.ab(ab)
  mic <- AMR::as.mic(mic)
  if (!is.null(mo)) {
    mo <- AMR::as.mo(mo)
  }

  if (is.null(distros)) {
    message(
    "No custom distributions provided to mic_uncensor using bootstrap method.
Using built-in epidemiological distributions from EUCAST:
https://www.eucast.org/bacteria/mic-and-zone-distributions-ecoffs/.
Note: Only Escherichia coli is currently supported.")
    distros <- ecoffs
  }

  ecoff_sub <- distros[distros['antibiotic'] == ab & distros['organism'] == mo,]
  ecoff_mics <- ecoff_sub[as.character(mic_range())]

  if (nrow(ecoff_sub) == 0) {
    warning("No ECOFFs available for ", ab)
    return(mic)
  }

  if (startsWith(as.character(mic), "<")) {
    ecoff_mics[AMR::as.mic(names(ecoff_mics)) > mic] <- 0
  }
  if (startsWith(as.character(mic), ">")) {
    ecoff_mics[AMR::as.mic(names(ecoff_mics)) <= mic] <- 0
  }

  ecoff_mics <- sapply(ecoff_mics, as.numeric)

  if (sum(ecoff_mics) == 0) {
    warning(ab, "appears to be below or above epidemiological distribution,
performing a simple mic uncensor.")
    return(mic_uncensor_simple(mic))
  }

  return(AMR::as.mic(sample(names(ecoff_mics), size = 1, replace=T, prob = ecoff_mics)))
},
USE.NAMES = F,
vectorize.args = c("mic", "ab", "mo"))

#' Censor MIC values
#'
#' @param mic MIC (coercible to AMR::as.mic)
#' @param ab antibiotic name (coercible to AMR::as.ab)
#' @param mo microorganism name (coercible to AMR::as.mo)
#' @param rules censor rules - named list of pathogen (in AMR::as.mo code) to
#' antibiotic (in AMR::as.ab code) to censoring rules. The censoring rules
#' should provide a min or max value to censor MICs to. See example for more.
#' @param max maximum concentration to censor to (default = Inf), will override
#' any rules provided
#' @param min minimum concentration to censor to (default = -Inf), will override
#' any rules provided
#'
#' @return censored MIC values (S3 mic class)
#'
#' @description
#' MIC datasets often arise from different laboratories or experimental
#' conditions. In practice, this means that there can be different levels of
#' censoring (<= and >) within the data. This function can be used to harmonise
#' the dataset to a single level of censoring. The function requires a set of
#' rules that specify the censoring levels (see example).
#'
#' @export
#'
#' @examples
#' example_rules <- list("B_ESCHR_COLI" = list(
#'   "AMK" = list(min = 2, max = 32),
#'   "CHL" = list(min = 4, max = 64),
#'   "GEN" = list(min = 1, max = 16),
#'   "CIP" = list(min = 0.015, max = 4),
#'   "MEM" = list(min = 0.016, max = 16),
#'   "AMX" = list(min = 2, max = 64),
#'   "AMC" = list(min = 2, max = 64),
#'   "FEP" = list(min = 0.5, max = 64),
#'   "CAZ" = list(min = 1, max = 128),
#'   "TGC" = list(min = 0.25, max = 1)
#'   ))
#'
#' mic_censor(AMR::as.mic(512),
#'            "AMK",
#'            "B_ESCHR_COLI",
#'            example_rules) == AMR::as.mic(">32")
mic_censor <- function(mic, ab = NULL,
                       mo = NULL,
                       rules = NULL,
                       max = Inf,
                       min = -Inf) {
  mic_censor_vectorize <- Vectorize(
    function(mic, ab, mo, rules) {
      mic <- AMR::as.mic(mic)
      min_thresh <- rules[[mo]][[ab]][["min"]]
      min_thresh <- ifelse(is.null(min_thresh),
                           -Inf,
                           min_thresh)
      max_thresh <- rules[[mo]][[ab]][["max"]]
      max_thresh <- ifelse(is.null(max_thresh),
                           Inf,
                           max_thresh)
      if (mic > max_thresh) {
        return(paste0(">", max_thresh))
      }
      if (mic < min_thresh) {
        return(paste0("<=", min_thresh))
      }
      return(mic)
    },
    USE.NAMES = FALSE,
    vectorize.args = c("mic", "ab", "mo")
  )

  rules_args <- c(ab, mo, rules)
  if (any(sapply(rules_args, is.null)) && any(!sapply(rules_args, is.null))) {
    stop("All of ab, mo, and rules must be provided or all must be NULL.")
  }

  censored <- AMR::as.mic(mic)

  if (all(!is.null(rules_args))) {
    censored <- AMR::as.mic(mic_censor_vectorize(mic, ab, mo, rules))
  }

  # apply overall min and max
  if (max != Inf) {
    max <- AMR::as.mic(max)
    censored[censored > max] <- paste0(">", max)
  }
  if (min != -Inf) {
    min <- AMR::as.mic(min)
    censored[censored < min] <- paste0("<=", min)
  }
  censored
}

#' Generate dilution series
#'
#' @param start starting (highest) concentration
#' @param dilutions number of dilutions
#' @param min minimum (lowest) concentration
#' @param precise force range to be high precision (not usually desired
#' behaviour)
#'
#' @return Vector of numeric concentrations
#' @export
#'
#' @examples
#' mic_range(128)
#' mic_range(128, dilutions = 21) # same results
mic_range <- function(start = 512,
                      dilutions = Inf,
                      min = 0.002,
                      precise = FALSE) {
  recursive_inner <- function(start,
                              dilutions,
                              min) {
    if (start[length(start)] < min) {
      return(utils::head(start, -1))
    }
    if (dilutions == 0) {
      return (start)
    } else {
      return(recursive_inner(c(start, start[length(start)] / 2),
                             dilutions - 1,
                             min))
    }
  }

  if (precise) {
    return(recursive_inner(start = start,
                           dilutions = dilutions,
                           min = min))
  }

  eucast_range <- c(512,
                    256,
                    128,
                    64,
                    32,
                    16,
                    8,
                    4,
                    2,
                    1,
                    0.5,
                    0.25,
                    0.125,
                    0.06,
                    0.03,
                    0.016,
                    0.008,
                    0.004,
                    0.002)
  filtered_range <- eucast_range[eucast_range >= min & eucast_range <= start]
  end_index <- ifelse(is.infinite(dilutions), length(filtered_range), dilutions)
  return(filtered_range[1:end_index])
}

#' Force MIC-like into MIC-compatible format
#'
#' @description
#' Convert a value that is "almost" an MIC into a valid MIC value.
#'
#' @param value vector of MIC-like values (numeric or character)
#' @param levels_from_AMR conform to AMR::as.mic levels
#' @param max_conc maximum concentration to force to
#' @param min_conc minimum concentration to force to
#' @param method method to use when forcing MICs (closest or round_up)
#' @param prefer where value is in between MIC (e.g., 24mg/L) chose the higher
#' MIC ("max") or lower MIC ("min"); only applies to method = "closest"
#' @param leq whether to force <= for lower censored values (i.e., <). If TRUE,
#' then all values below the limit of detection are converted to <=. If FALSE,
#' then they are converted to <. If NULL, they are not changed.
#' @param geq whether to force >= for higher censored values (i.e., >). If TRUE,
#' then all values above the limit of detection are converted to >=. If FALSE,
#' then they are converted to >. If NULL, they are not changed.
#'
#' @details
#' Some experimental or analytical conditions measure MIC (or surrogate) in a
#' way that does not fully conform to traditional MIC levels
#' (i.e., concentrations). This function allows these values to be coerced into
#' an MIC value that is compatible with the AMR::mic class. When using method =
#' "closest", the function will choose the closest MIC value to the input value
#' (e.g., 2.45 will be coerced to 2). When using method = "round up", the
#' function will round up to the next highest MIC value (e.g., 2.45 will be
#' coerced to 4). "Round up" is technically the correct approach if the input
#' value was generated from an experiment that censored between concentrations
#' (e.g., broth or agar dilution). However, "closest" may be more appropriate in
#' some cases.
#'
#' Please note that this function will not make any changes to censored values
#' (beyond some simple cleaning, e.g., <==2 is converted to <=2). This is because
#' it is not possible to make assumptions about censored data.
#'
#' The `leq` and `geq` arguments convert censored values to <= or >=. When MIC
#' is measured using a an inhibitory dilution method, the lower limit should be
#' reported as <= (since the lowest dilution could be inhibitory itself), and
#' the upper limit should be reported as > (growth in the highest dilution means
#' that it is not an inhibitory concentration). The default values for `leq`
#' and `geq` enforce this.
#'
#' @return AMR::as.mic compatible character
#' @export
#'
#' @examples
#' force_mic(c("2.32", "<4.12", ">1.01"))
force_mic <- function(value,
                      levels_from_AMR = FALSE,
                      max_conc = 512,
                      min_conc = 0.002,
                      method = "closest",
                      prefer = 'max',
                      leq = TRUE,
                      geq = NULL) {

  if (levels_from_AMR) {
    mic_levels <- levels(AMR::as.mic(NA))
  } else {
    mic_levels <- mic_range(start = max_conc, min = min_conc)
    mic_levels <- c(max(mic_levels) * 2,
                    mic_levels,
                    min(mic_levels) / 2)
    mic_levels <- c(mic_levels,
                    paste0(">", mic_levels),
                    paste0("<", mic_levels))
  }

  output <- rep(NA_character_, length(value))
  for (i in seq_along(value)) {
    prefix <- NULL
    inner_val <- value[i]

    if (is.na(inner_val)) {
      next
    }

    stopifnot(any(c(AMR::is.mic(inner_val), is.numeric(inner_val), is.character(inner_val), is.na(inner_val))))

    if (AMR::is.mic(inner_val)) {
      inner_val <- as.character(inner_val)
    }
    if (is.numeric(inner_val)) {
      appropriate_levels <- subset(mic_levels, !stringr::str_detect(mic_levels, "[^0-9.]"))
      censored <- FALSE
    }
    if (is.character(inner_val)){
      if (stringr::str_detect(inner_val, "<") &
          stringr::str_detect(inner_val, "=")) {
        appropriate_levels <- subset(mic_levels, stringr::str_detect(mic_levels, "<"))
        appropriate_levels <- stringr::str_remove_all(appropriate_levels, "[^0-9.]")
        censored <- TRUE
        prefix <- "<="
      } else if (stringr::str_detect(inner_val, ">") &
                 stringr::str_detect(inner_val, "=")) {
        appropriate_levels <- subset(mic_levels, stringr::str_detect(mic_levels, ">"))
        appropriate_levels <- stringr::str_remove_all(appropriate_levels, "[^0-9.]")
        censored <- TRUE
        prefix <- ">="
      } else if (stringr::str_detect(inner_val, "<")) {
        appropriate_levels <- subset(mic_levels, stringr::str_detect(mic_levels, "<"))
        appropriate_levels <- stringr::str_remove_all(appropriate_levels, "[^0-9.]")
        censored <- TRUE
        prefix <- "<"
      } else if (stringr::str_detect(inner_val, ">")) {
        appropriate_levels <- subset(mic_levels, stringr::str_detect(mic_levels, ">"))
        appropriate_levels <- stringr::str_remove_all(appropriate_levels, "[^0-9.]")
        censored <- TRUE
        prefix <- ">"
      } else {
        appropriate_levels <- subset(mic_levels, !stringr::str_detect(mic_levels, "[^0-9.]"))
        censored <- FALSE
      }
      inner_val <- stringr::str_remove_all(inner_val, "[^0-9.]")
      inner_val <- as.numeric(inner_val)
    }
    mic_vector <- sort(as.numeric(appropriate_levels))

    if (!method %in% c("closest", "round up")) {
      stop("Method must be closest or round up")
    }

    if (method == "closest") {
      positions <- which(abs(mic_vector - inner_val) == min(abs(mic_vector - inner_val)))
    }

    if (method == "round up") {
      if (inner_val %in% mic_vector) {
        positions <- which(mic_vector == inner_val)
      } else {
        positions <- min(which(inner_val < mic_vector))
      }
    }

    if (censored) {
      # no changes made to censored values
      output[i] <- paste0(prefix, inner_val)
      if (!is.null(leq)) {
        if (leq) {
          if (stringr::str_detect(prefix, "<")) {
            output[i] <- paste0("<=", inner_val)
          }
        } else {
          if (stringr::str_detect(prefix, "<")) {
            output[i] <- paste0("<", inner_val)
          }
        }
      }
      if (!is.null(geq)) {
        if (geq) {
          if (stringr::str_detect(prefix, ">")) {
            output[i] <- paste0(">=", inner_val)
          }
        } else {
          if (stringr::str_detect(prefix, ">")) {
            output[i] <- paste0(">", inner_val)
          }
        }
      }
      next
    }

    if (length(positions) == 1) {
      output[i] <- paste0(prefix, mic_vector[positions])
    }

    if (prefer == 'min') {
      output[i] <- paste0(prefix, mic_vector[min(positions)])
    }

    if (inner_val > max_conc * 2) {
      warning(paste("Value greater than supported MIC:", inner_val))
      output[i] <- NA_character_
    } else if (inner_val < min_conc / 2) {
      warning(paste("Value less than supported MIC:", inner_val))
      output[i] <- NA_character_
    } else {
      output[i] <- paste0(prefix, mic_vector[max(positions)])
    }
  }
  return(output)
}

#' Convert MIC or Disk Diffusion to SIR, vectorised over antimicrobials
#'
#' @param mic vector of MIC values
#' @param mo vector of microorganism names
#' @param ab vector of antibiotic names
#' @param accept_ecoff if TRUE, ECOFFs will be used when no clinical breakpoints are available
#' @param ... additional arguments that are passed to AMR::as.sir
#'
#' @return S3 sir values
#' @export
#' @description
#' The AMR::as.sir function is not vectorised over antimicrobials. This function
#' provides vectorisation over antimicrobials. Due to the overhead of running
#' AMR::as.sir, this function tries to be efficient by only running AMR::as.sir
#' as little as necessary.
#'
#' @examples
#' mic <- c("<0.25", "8", "64", ">64")
#' mo <- c("B_ESCHR_COLI", "B_ESCHR_COLI", "B_ESCHR_COLI", "B_ESCHR_COLI")
#' ab <- c("AMK", "AMK", "AMK", "AMK")
#' as.sir_vectorised(mic, mo, ab)

#' # using different microorganisms and antibiotics
#' mic <- c("<0.25", "8", "64", ">64")
#' mo <- c("B_ESCHR_COLI", "B_ESCHR_COLI", "B_PROTS_MRBL", "B_PROTS_MRBL")
#' ab <- c("AMK", "AMK", "CIP", "CIP")
#' as.sir_vectorised(mic, mo, ab)
as.sir_vectorised <- function(mic, mo, ab, accept_ecoff = FALSE, ...) {
  if (length(unique(ab)) == 1) {
    ab <- ab[[1]]
    output <- AMR::as.sir(mic, mo = mo, ab = ab, ...)
    if (accept_ecoff & any(is.na(output))) {
      output[is.na(output)] <- AMR::as.sir(mic[is.na(output)],
                                           mo = mo[is.na(output)],
                                           ab = ab,
                                           breakpoint_type = "ECOFF",
                                           ...)
    }
    return(output)
  }

  output <- purrr::pmap_vec(list(mic,
                       mo,
                       ab), \(x, mo, ab) AMR::as.sir(x, mo = mo, ab = ab))
  output <- purrr::pmap_vec(list(mic,
                                 mo,
                                 ab,
                                 output), \(x, mo, ab, sir) {
                                   if (is.na(sir)) {
                                     if (accept_ecoff) {
                                       message(paste("No clinical breakpoint for:", ab, mo, "using ECOFF"))
                                       AMR::as.sir(x, mo = mo, ab = ab,
                                                   breakpoint_type = "ECOFF")
                                     } else {
                                       message(paste("No clinical breakpoint for:", ab, mo, "(use accept_ecoff to use ECOFF)"))
                                       sir
                                     }
                                   } else {
                                     sir
                                   }
                                 })
  output
}

#' S breakpoint for MIC
#'
#' @param mo mo name (coerced using AMR::as.mo)
#' @param ab ab name (coerced using AMR::as.ab)
#' @param accept_ecoff if TRUE, ECOFFs will be used when no clinical breakpoints are available
#' @param ... additional arguments to pass to AMR::as.sir, which is used to
#' calculate the S breakpoint
#'
#' @return MIC value
#' @export
#'
#' @examples
#' mic_s_breakpoint("B_ESCHR_COLI", "AMK")
#' mic_s_breakpoint("B_ESCHR_COLI", "CHL", accept_ecoff = TRUE)
mic_s_breakpoint <- function(mo, ab, accept_ecoff = FALSE, ...) {
  mic_range <- mic_range()
  sir_range <- AMR::as.sir(AMR::as.mic(mic_range),
                           mo = mo,
                           ab = ab,
                           ...)
  if (any(is.na(sir_range))) {
    if (accept_ecoff) {
      sir_range <- AMR::as.sir(AMR::as.mic(mic_range),
                               mo = mo,
                               ab = ab,
                               breakpoint_type = "ECOFF",
                               ...)
    } else {
      stop("No clinical breakpoints available for ", ab, mo, "(consider using accept_ecoff if appropriate.")
    }
  }
  for (i in seq_along(sir_range)) {
    if (sir_range[i] == "S") {
      return(AMR::as.mic(mic_range[i]))
    }
  }
}

#' R breakpoint for MIC
#'
#' @param mo mo name (coerced using AMR::as.mo)
#' @param ab ab name (coerced using AMR::as.ab)
#' @param accept_ecoff if TRUE, ECOFFs will be used when no clinical breakpoints are available
#' @param ... additional arguments to pass to AMR::as.sir, which is used to
#' calculate the R breakpoint
#'
#' @return MIC value
#' @export
#'
#' @examples
#' mic_r_breakpoint("B_ESCHR_COLI", "AMK")
#' mic_r_breakpoint("B_ESCHR_COLI", "CHL", accept_ecoff = TRUE)
mic_r_breakpoint <- function(mo, ab, accept_ecoff = FALSE, ...) {
  mic_range <- rev(mic_range())
  sir_range <- AMR::as.sir(AMR::as.mic(mic_range),
                           mo = mo,
                           ab = ab,
                           ...)
  if (any(is.na(sir_range))) {
    if (accept_ecoff) {
      sir_range <- AMR::as.sir(AMR::as.mic(mic_range),
                               mo = mo,
                               ab = ab,
                               breakpoint_type = "ECOFF",
                               ...)
    } else {
      stop("No clinical breakpoints available for ", ab, mo, "(consider using accept_ecoff if appropriate.")
    }
  }
  for (i in seq_along(sir_range)) {
    if (sir_range[i] == "R") {
      return(AMR::as.mic(mic_range[i-1]))
    }
  }
}

match_levels <- function(x, match_to) {
  if (!AMR::is.mic(x) | !AMR::is.mic(match_to)) {
    stop("Both x and match_to must be AMR::mic objects")
  }
  all_levels <- levels(AMR::as.mic(NA))

  keep_these <- levels(droplevels(match_to))
  drop_these <- all_levels[!all_levels %in% keep_these]

  x <- forcats::fct_drop(x, only = drop_these)
  class(x) <- append(class(x), "mic", after = 0)
  x
}

#' Fill MIC dilution levels
#'
#' @param x MIC vector
#' @param cap_upper If True, will the top level will be the highest MIC dilution in x
#' @param cap_lower If True, will the bottom level will be the lowest MIC dilution in x
#' @param as.mic By default, returns an ordered factor. Set as.mic = TRUE to return as AMR::mic
#'
#' @return ordered factor (or AMR::mic if as.mic = TRUE)
#' @export
#'
#' @examples
#' # use in combination with droplevels to clean up levels:
#' x <- AMR::as.mic(c("<0.25", "8", "64", ">64"))
#' x <- droplevels(x)
#' fill_dilution_levels(x)
fill_dilution_levels <- function(x,
                                 cap_upper = TRUE,
                                 cap_lower = TRUE,
                                 as.mic = TRUE) {
  mic_conc_range <- AMR::as.mic(mic_range())
  if (cap_lower) {
    mic_conc_range <- mic_conc_range[mic_conc_range >= min(AMR::as.mic(x))]
  }
  if (cap_upper) {
    mic_conc_range <- mic_conc_range[mic_conc_range <= max(AMR::as.mic(x))]
  }

  x <- forcats::fct_expand(x, as.character(mic_conc_range))

  new_levels <- levels(x)[order(match(levels(x),
                                      as.character(levels(AMR::as.mic(NA)))))]

  x <- ordered(x, levels = new_levels)
  if (as.mic) {
    class(x) <- append(class(x), "mic", after = 0)
  }
  x
}
