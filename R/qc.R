#' Check that MIC is within QC range
#'
#' @param measurement measured QC MIC
#' @param strain control strain identifier (usually ATCC)
#' @param ab antibiotic name (will be coerced to AMR::as.ab)
#' @param ignore_na ignores NA (returns TRUE)
#' @param guideline Guideline to use (EUCAST or CLSI)
#' @param year Guideline year (version)
#' @return logical vector
#'
#' @description
#' Check whether MIC values are within acceptable range for
#' quality control (QC). Every MIC experiment should include a control strain
#' with a known MIC. The results of the experiment are only valid if the control
#' strain MIC falls within the acceptable range. This function checks whether
#' an MIC result is within the acceptable range given: 1) a control strain
#' (usually identified as an ATCC or NCTC number), 2) an antibiotic name, and 3)
#' a guideline (EUCAST or CLSI). The acceptable range is defined by 'QC_table',
#' which is a dataset which is loaded with this package.
#'
#' The source of the QC values is the WHONET QC Ranges and Targets available from
#' the 'Antimicrobial Resistance Test Interpretation Engine' (AMRIE) repository:
#' https://github.com/AClark-WHONET/AMRIE
#'
#' @export
#' @references
#' O’Brien TF, Stelling JM. WHONET: An Information System for Monitoring
#' Antimicrobial Resistance. Emerg Infect Dis. 1995 Jun;1(2):66–66.

#' @examples
#' qc_in_range(AMR::as.mic(0.5), 25922, "GEN") == TRUE
#' qc_in_range(AMR::as.mic(8.0), 25922, "GEN") == FALSE
qc_in_range <- Vectorize(
  function(measurement,
           strain,
           ab,
           ignore_na = TRUE,
           guideline = "EUCAST",
           year = "2023") {
    if (is.na(measurement) | is.na(strain) | is.na(ab)) {
      if (ignore_na) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
    if (!inherits(measurement, "mic") & !inherits(measurement, "disk")) {
      stop("Measurement must be AMR::mic or AMR::disk")
    }
    ab <- as.character(AMR::as.ab(ab))

    strain <- tolower(strain)
    if (!startsWith(strain, "atcc")) {
      strain <- paste0("atcc", strain)
    }
    qc_subset <- QC_table[QC_table$STRAIN == strain &
                            QC_table$ANTIBIOTIC == ab &
                            QC_table$GUIDELINE == guideline &
                            QC_table$YEAR == year &
                            QC_table$METHOD == "MIC", ]

    if (nrow(qc_subset) == 0) {
      # no QC info in table
      return(NA)
    }

    stopifnot(nrow(qc_subset) == 1)

    if (inherits(measurement, "mic")) {
      measurement <- mic_uncensor(measurement)

      if (measurement < AMR::as.mic(qc_subset[, "MINIMUM"])) {
        return(FALSE)
      }
      if (measurement > AMR::as.mic(qc_subset[, "MAXIMUM"])) {
        return(FALSE)
      }
      return(TRUE)
    }
  },
  vectorize.args = c("measurement", "strain", "ab")
)

#' Check that QC measurement is at the required target
#' `r lifecycle::badge('experimental')`
#'
#' @param measurement measured QC MIC
#' @param strain control strain identifier (usually ATCC)
#' @param ab antibiotic name (will be coerced to AMR::as.ab)
#' @param ignore_na ignores NA (returns TRUE)
#' @param guideline Guideline to use (EUCAST or CLSI)
#' @param year Guideline year (version)
#' @description
#' MIC experiments should include a control strain with a known MIC.
#' The MIC result for the control strain should be a particular target MIC. This
#' function checks whether the target MIC was achieved given: 1) a control
#' strain (usually identified as an ATCC or NCTC number), 2) an antibiotic name,
#' and 3) a guideline (EUCAST or CLSI).
#'
#' Since QC target values are currently not publicly available in an easy to
#' use format, this function takes a pragmatic approach -- for most antibiotics
#' and QC strains, the target is assumed to be the midpoint of the acceptable
#' range. This approximation is not necessarily equal to the QC target reported
#' by guideline setting bodies such as EUCAST. Therefore, this function is
#' considered experimental and should be used with caution.
#'
#' This function can be used alongnside qc_in_range(), which checks whether the
#' MIC is within the acceptable range.
#'
#' The source of the QC values is the WHONET QC Ranges and Targets available from
#' the 'Antimicrobial Resistance Test Interpretation Engine' (AMRIE) repository:
#' https://github.com/AClark-WHONET/AMRIE
#'
#' @return logical vector
#' @export
#'
#' @references
#' O’Brien TF, Stelling JM. WHONET: An Information System for Monitoring
#' Antimicrobial Resistance. Emerg Infect Dis. 1995 Jun;1(2):66–66.
#'
#' @examples
#' qc_on_target(AMR::as.mic(0.5), 25922, "GEN") == TRUE
qc_on_target <- Vectorize(
  function(measurement,
           strain,
           ab,
           ignore_na = TRUE,
           guideline = "EUCAST",
           year = "2023") {
    if (is.na(measurement) | is.na(strain) | is.na(ab)) {
      if (ignore_na) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
    if (!inherits(measurement, "mic") & !inherits(measurement, "disk")) {
      stop("Measurement must be AMR::mic or AMR::disk")
    }
    ab <- as.character(AMR::as.ab(ab))

    strain <- tolower(strain)
    if (!startsWith(strain, "atcc")) {
      strain <- paste0("atcc", strain)
    }

    qc_subset <- QC_table[QC_table$STRAIN == strain &
                            QC_table$ANTIBIOTIC == ab &
                            QC_table$GUIDELINE == guideline &
                            QC_table$YEAR == year & QC_table$METHOD == "MIC", ]

    if (nrow(qc_subset) == 0) {
      # no QC info in table
      return(NA)
    }

    stopifnot(nrow(qc_subset) == 1)

    if (inherits(measurement, "mic")) {
      measurement <- mic_uncensor(measurement)

      if (measurement < qc_subset[, "MINIMUM_TARGET"]) {
        return(FALSE)
      }
      if (measurement > qc_subset[, "MAXIMUM_TARGET"]) {
        return(FALSE)
      }
      return(TRUE)
    }
  },
  vectorize.args = c("measurement", "strain", "ab")
)

#' Standardise MIC to control strain
#' `r lifecycle::badge('experimental')`
#'
#' @param test_measurement Measured MIC to standardise
#' @param qc_measurement Measured QC MIC to standardise to
#' @param strain control strain identifier (usually ATCC)
#' @param ab antibiotic name (will be coerced to AMR::as.ab)
#' @param prefer_upper Where the target MIC is a range, prefer the upper value
#' in the range
#' @param ignore_na Ignore NA (returns AMR::NA_mic_)
#' @param guideline Guideline to use (EUCAST or CLSI)
#' @param year Guideline year (version)
#' @param force Force into MIC-compatible format after standardisation
#'
#' @return AMR::mic vector
#' @export
#'
#' @description
#' MIC experiments are generally quality-controlled by including a control strain
#' with a known MIC. The MIC result for the control strain should be a particular
#' target MIC, or at least within an acceptable range. This function standardises
#' a measured MIC to the target MIC given: 1) a control strain (usually identified
#' as an ATCC or NCTC number), 2) an antibiotic name, and 3) a guideline (EUCAST
#' or CLSI). The definition of standardisation in this context is to adjust the
#' measured MIC based on the QC MIC. This is based on the following principles
#' and assumption:
#'
#' 1. A measured MIC is composed of two components: the true MIC and a
#' measurement error. The measurement error is considered to be inevitable when
#' measuring MICs, and is likely to be further composed of variability in
#' laboratory conditions and operator interpretation.
#' 2. It is assumed that the MIC of the control strain in the experiment has
#' also been affected by this error.
#'
#' The standardisation applied by this function uses the measured QC strain
#' MIC as a reference point, and scales the rest of the MICs to this reference.
#' In general, this means that the MICs are doubled or halved, depending on the
#' result of the QC MIC. A worked example is provided below and illustrates
#' the transformation that this function applies.
#'
#' There is no current evidence base for this approach, therefore, this
#' function is considered experimental and should be used with caution.
#'
#' @examples
#' # Ref strain QC MIC for GEN is 0.5
#' standardise_mic(
#'   test_measurement = c(AMR::as.mic(">8.0"),  # QC = 1, censored MIC remains censored
#'                        AMR::as.mic(4.0),  # QC = 0.5 which is on target, so stays same
#'                        AMR::as.mic(2),  # QC = 1, so scaled down to 1
#'                        AMR::as.mic(2)),  # QC = 0.25, so scaled up to 8
#'   qc_measurement = c(AMR::as.mic(1),
#'                      AMR::as.mic(0.5),
#'                      AMR::as.mic(1),
#'                      AMR::as.mic(0.25)),
#'   strain = 25922,
#'   ab = AMR::as.ab("GEN"))
standardise_mic <- function(test_measurement,
                            qc_measurement,
                            strain,
                            ab,
                            prefer_upper = FALSE,
                            ignore_na = TRUE,
                            guideline = "EUCAST",
                            year = "2023",
                            force = TRUE) {
  standardise_mic_vectorised <- Vectorize(
    function(test_measurement,
             qc_measurement,
             strain,
             ab,
             prefer_upper,
             ignore_na,
             guideline,
             year,
             force) {
      if (is.na(test_measurement) | is.na(qc_measurement) | is.na(strain) | is.na(ab)) {
        if (ignore_na) {
          return(NA)
        } else {
          return(NA)
        }
      }
      if (!inherits(test_measurement, "mic") | !inherits(qc_measurement, "mic")) {
        stop("Measurements must be AMR::mic")
      }
      ab <- as.character(AMR::as.ab(ab))

      mic_char <- as.character(test_measurement)
      if (grepl("<|>", mic_char)) {
        return(as.character(test_measurement))
      }

      if (qc_on_target(qc_measurement, strain, ab)) {
        return(as.character(test_measurement))
      }

      if (!qc_in_range(qc_measurement, strain, ab)) {
        warning("QC not in range, standardise_mic may not be appropriate")
      }

      strain <- tolower(strain)
      if (!startsWith(strain, "atcc")) {
        strain <- paste0("atcc", strain)
      }
      qc_subset <- QC_table[QC_table$STRAIN == strain & QC_table$ANTIBIOTIC == ab & QC_table$GUIDELINE == guideline & QC_table$YEAR == year & QC_table$METHOD == "MIC", ]

      if (nrow(qc_subset) == 0) {
        # no QC info in table
        return(NA)
      }

      stopifnot(nrow(qc_subset) == 1)

      if (qc_subset[1, "TARGET_ESTIMATED"]) {
        warning(paste("Target for", ab, "is estimated as mid-point between upper
and lower range. EUCAST/CLSI guidance may be different."))
      }

      scale_factor_lower <- qc_subset[, "MINIMUM_TARGET"]
      scale_factor_upper <- qc_subset[, "MAXIMUM_TARGET"]
      if (scale_factor_lower == scale_factor_upper) {
        scale_factor <- scale_factor_lower
      } else{
        if (prefer_upper) {
          scale_factor <- scale_factor_upper
        } else {
          scale_factor <- scale_factor_lower
        }
      }

      as.character(test_measurement * (scale_factor / qc_measurement))
    },
    vectorize.args = c("test_measurement",
                       "qc_measurement",
                       "strain",
                       "ab")
  )
  if (force) {
    return(
      AMR::as.mic(force_mic(standardise_mic_vectorised(test_measurement,
                                                       qc_measurement,
                                                       strain,
                                                       ab,
                                                       prefer_upper,
                                                       ignore_na,
                                                       guideline,
                                                       year)
      )
      ))
  }
  AMR::as.mic(standardise_mic_vectorised(test_measurement,
                                         qc_measurement,
                                         strain,
                                         ab,
                                         prefer_upper,
                                         ignore_na,
                                         guideline,
                                         year))
}
