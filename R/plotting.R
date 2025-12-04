prepare_mic_validation_plotting_data <- function(x, match_axes, add_missing_dilutions) {
  x <- as.data.frame(x)
  # keep only columns needed for plotting
  if (!"ab" %in% colnames(x)) {
    x <- x[,c("gold_standard", "test", "essential_agreement")]
  } else {
    x <- x[,c("gold_standard", "test", "essential_agreement", "ab")]
  }

  if (match_axes) {
    x[["gold_standard"]] <- match_levels(x[["gold_standard"]], match_to = x[["test"]])
    x[["test"]] <- match_levels(x[["test"]], match_to = x[["gold_standard"]])

    if (add_missing_dilutions) {
      x[["gold_standard"]] <- fill_dilution_levels(x[["gold_standard"]],
                                                   cap_lower = TRUE,
                                                   cap_upper = TRUE)
      x[["test"]] <- fill_dilution_levels(x[["test"]],
                                          cap_lower = TRUE,
                                          cap_upper = TRUE)
    }

    if (length(levels(x[["gold_standard"]])) > length(levels(x[["test"]]))) {
      # after dilution filling, levels may not yet match, force another match
      x[["test"]] <- forcats::fct_expand(x[["test"]],
                                          as.character(levels(x[["gold_standard"]])))
      x[["test"]] <- forcats::fct_relevel(x[["test"]],
                                          levels(x[["gold_standard"]]))
    }

    if (length(levels(x[["test"]])) > length(levels(x[["gold_standard"]]))) {
      x[["gold_standard"]] <- forcats::fct_expand(x[["gold_standard"]],
                                                  as.character(levels(x[["test"]])))
      x[["gold_standard"]] <- forcats::fct_relevel(x[["gold_standard"]],
                                                  levels(x[["test"]]))
    }
  }

  # temp fix - drop use of mic class as a patch to allow AMR v3 compatibility
  x[["gold_standard"]] <- factor(x[["gold_standard"]])
  x[["test"]] <- factor(x[["test"]])

  x
}

#' @export
plot.single_ab_validation <- function(x,
                                      match_axes = TRUE,
                                      add_missing_dilutions = TRUE,
                                      ...) {
  x_df <- prepare_mic_validation_plotting_data(x, match_axes, add_missing_dilutions)

  p <- x_df |>
    dplyr::group_by(.data[["gold_standard"]],
                    .data[["test"]],
                    .data[["essential_agreement"]]) |>
    dplyr::summarise(n = dplyr::n()) |>
    dplyr::rename(`EA` = .data[["essential_agreement"]]) |>
    ggplot2::ggplot(ggplot2::aes(x = .data[["gold_standard"]],
                                 y = .data[["test"]],
                                 fill = .data[["n"]],
                                 color = .data[["EA"]])) +
    ggplot2::geom_tile(alpha=1, show.legend = TRUE) +
    ggplot2::geom_text(ggplot2::aes(label=.data[["n"]]), show.legend = TRUE) +
    ggplot2::scale_fill_gradient(low="white", high="#009194") +
    ggplot2::scale_fill_manual(values=c("red", "black"), aesthetics = "color", drop = FALSE) +
    ggplot2::guides(color=ggplot2::guide_legend(override.aes=list(fill=NA))) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::xlab("Gold standard MIC (mg/L)") +
    ggplot2::ylab("Test (mg/L)")

  if (match_axes) {
    p <- p + ggplot2::scale_x_discrete(drop = FALSE)
    p <- p + ggplot2::scale_y_discrete(drop = FALSE)
  }

  # if ("ab" %in% names(x) & "mo" %in% names(x)) {
  #     bpoints <- AMR::clinical_breakpoints
  #     p <- p + ggplot2::geom_hline(yintercept = AMR::as.mic(bpoints[]))
  # }
  p
}

#' @export
plot.multi_ab_validation <- function(x,
                                     match_axes = TRUE,
                                     add_missing_dilutions = TRUE,
                                     facet_wrap_ncol = NULL,
                                     facet_wrap_nrow = NULL,
                                     ...) {
  if (is.null(facet_wrap_ncol) && is.null(facet_wrap_nrow)) {
    return(plot.single_ab_validation(x, match_axes, add_missing_dilutions, ...))
  }

  x_df <- prepare_mic_validation_plotting_data(x, match_axes, add_missing_dilutions)

  p <- x_df |>
    dplyr::group_by(.data[["gold_standard"]],
                    .data[["test"]],
                    .data[["essential_agreement"]],
                    .data[["ab"]]) |>
    dplyr::mutate(ab = AMR::ab_name(AMR::as.ab(as.character(.data[["ab"]])))) |>
    dplyr::mutate(ab = dplyr::if_else(is.na(.data[["ab"]]), "unknown", .data[["ab"]])) |>
    dplyr::summarise(n = dplyr::n()) |>
    dplyr::rename(`EA` = .data[["essential_agreement"]]) |>
    ggplot2::ggplot(ggplot2::aes(x = .data[["gold_standard"]],
                                 y = .data[["test"]],
                                 fill = .data[["n"]],
                                 color = .data[["EA"]])) +
    ggplot2::geom_tile(alpha=1) +
    ggplot2::geom_text(ggplot2::aes(label=.data[["n"]])) +
    ggplot2::scale_fill_gradient(low="white", high="#009194") +
    ggplot2::scale_fill_manual(values=c("red", "black"), aesthetics = "color", drop = FALSE)

    if (any(!is.null(c(facet_wrap_ncol, facet_wrap_nrow)))) {
      p <- p + ggh4x::facet_wrap2(~ .data[["ab"]],
                                     nrow = facet_wrap_nrow,
                                     ncol = facet_wrap_ncol,
                                     axes = "all")
    }

    p <- p +
    ggplot2::guides(color=ggplot2::guide_legend(override.aes=list(fill=NA))) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::xlab("Gold standard MIC (mg/L)") +
    ggplot2::ylab("Test (mg/L)")

  if (match_axes) {
    p <- p + ggplot2::scale_x_discrete(drop = FALSE)
    p <- p + ggplot2::scale_y_discrete(drop = FALSE)
  }
  p
}

#' Plot MIC validation results
#'
#' @param x object generated using compare_mic
#' @param match_axes Same x and y axis
#' @param add_missing_dilutions Axes will include dilutions that are not
#' @param facet_wrap_ncol Facet wrap into n columns by antimicrobial (optional,
#' only available when more than one antimicrobial in validation)
#' @param facet_wrap_nrow Facet wrap into n rows by antimicrobial (optional,
#' only available when more than one antimicrobial in validation)
#' represented in the data, based on a series of dilutions generated using mic_range().
#' @param ... additional arguments
#'
#' @return ggplot object
#'
#' @export
#'
#' @examples
#' gold_standard <- c("<0.25", "8", "64", ">64")
#' test <- c("<0.25", "2", "16", "64")
#' val <- compare_mic(gold_standard, test)
#' plot(val)
#'
#' # if the validation contains multiple antibiotics, i.e.,
#' ab <- c("CIP", "CIP", "AMK", "AMK")
#' val <- compare_mic(gold_standard, test, ab)
#' # the following will plot all antibiotics in a single plot (pooled results)
#' plot(val)
#' # use the faceting arguments to split the plot by antibiotic
#' plot(val, facet_wrap_ncol = 2)
plot.mic_validation <- function(x,
                                match_axes = TRUE,
                                add_missing_dilutions = TRUE,
                                facet_wrap_ncol = NULL,
                                facet_wrap_nrow = NULL,
                                ...) {
  # Fallback for objects without specific class
  if (!is.null(x$ab) && length(unique(x$ab)) > 1) {
    plot.multi_ab_validation(x,
                             match_axes = match_axes,
                             add_missing_dilutions = add_missing_dilutions,
                             facet_wrap_ncol = facet_wrap_ncol,
                             facet_wrap_nrow = facet_wrap_nrow,
                             ...)
  } else {
    plot.single_ab_validation(x,
                              match_axes = match_axes,
                              add_missing_dilutions = add_missing_dilutions,
                              ...)
  }
}
