#' Convert genomes to kmers in libsvm format
#'
#' @param source_dir directory containing genomes
#' @param target_dir target directory to store kmers in libsvm format
#' @param k k-mer length
#' @param canonical only count canonical kmers
#' @param squeeze remove non-canonical kmers
#' @param ext file extension to filter
#'
#' @description
#' Raw genome data (pre- or post-assembly) is usually transformed by k-mer
#' counting prior to machine learning (ML). XGBoost is a popular ML algorithm
#' for this problem, due to its scalability to high dimensional data. This
#' function converts genomes to k-mer counts stored in XGBoost's preferred
#' format, libsvm. Further information on the libsvm format is available at
#' \url{https://xgboost.readthedocs.io/en/stable/tutorials/input_format.html}.
#' Briefly, libsvm is effectively a text file that stores data points as
#' x:y pairs, where x is the feature index, and y is the feature value. Each
#' observation is stored on its own line, with the first column reserved for
#' labels. Labels can be provided later, during data import.
#'
#' This function converts each individual genome to an individual libsvm
#' text file of k-mer counts (therefore, each .txt file will be 1 line long).
#' This function supports parallel processing using the by setting an appropriate
#' \code{future::plan()} (usually \code{future::multisession}) ---
#' each genome is processed in parallel. To monitor progress, use the
#' \code{progressr} package by wrapping the function in
#' \code{\link[progressr]{with_progress}}.
#'
#' Although XGBoost can load a multiple .txt (libsvm) files by providing the
#' directory as an input, this is generally not recommended as order of
#' import cannot be guaranteed and probably depends on filesystem. Instead,
#' it is recommended that this function is combined with
#' [split_and_combine_files()] which generates a single .txt file (with the
#' order of observations guaranteed and stored in a .csv file).
#'
#' @return TRUE if successful
#'
#' @examples
#' set.seed(123)
#' # create 10 random DNA files
#' tmp_dir <- tempdir()
#' # remove any existing .fna files
#' file.remove(
#'  list.files(tmp_dir, pattern = "*.fna", full.names = TRUE)
#' )
#' for (i in 1:10) {
#' writeLines(paste0(">", i, "\n", paste0(sample(c("A", "T", "C", "G"),
#'  100, replace = TRUE), collapse = "")), file.path(tmp_dir, paste0(i, ".fna")))
#' }
#'
#' tmp_target_dir <- file.path(tmp_dir, "kmers")
#' unlink(tmp_target_dir, recursive = TRUE)
#'
#' # convert genomes to k-mers
#' future::plan(future::sequential)  # use multisession for parallel processing
#' progressr::with_progress(
#'   genomes_to_kmer_libsvm(tmp_dir, tmp_target_dir, k = 3)
#' )
#'
#' # check the output
#' list.files(tmp_target_dir)
#' readLines(list.files(tmp_target_dir, full.names = TRUE)[1])
#'
#' @seealso to convert a single genome, use [genome_to_libsvm()]
#' @export
genomes_to_kmer_libsvm <- function(source_dir,
                                   target_dir,
                                   k = 3,
                                   canonical = TRUE,
                                   squeeze = FALSE,
                                   ext = ".fna") {
  lifecycle::deprecate_warn(
    when = "1.2.0",
    what = "genomes_to_kmer_libsvm()",
    with = "faLearn::genomes_to_kmer_libsvm()",
    details = "This function has been moved to the faLearn package."
  )
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
  }
  ext <- gsub("^\\.", "", ext)
  genome_paths <- list.files(source_dir,
                             pattern = paste0("*.", ext),
                             full.names = TRUE,
                             ignore.case = TRUE)
  p <- progressr::progressor(along = genome_paths)
  future.apply::future_lapply(genome_paths, \(x) {
    genome_to_libsvm(as.character(Biostrings::readDNAStringSet(x)),
                    file.path(normalizePath(target_dir),
                              paste0(strip_filename(x), ".txt")),
                    k = k,
                    canonical = canonical,
                    squeeze = squeeze)
    p(glue::glue("Completed: {basename(x)}"))
  }, future.seed = TRUE)
  return(TRUE)
}

#' @rdname reverse_complement
#' @export
reverse_complement <- function(dna) {
  lifecycle::deprecate_warn(
    when = "1.2.0",
    what = "reverse_complement()",
    with = "faLearn::reverse_complement()",
    details = "This function has been moved to the faLearn package."
  )
  .Call(`_MIC_reverse_complement`, dna)
}

#' @rdname kmers
#' @export
kmers <- function(x, k = 3L, simplify = FALSE, canonical = TRUE, squeeze = FALSE, anchor = TRUE, clean_up = TRUE, key_as_int = FALSE, starting_index = 1L) {
  lifecycle::deprecate_warn(
    when = "1.2.0",
    what = "kmers()",
    with = "faLearn::kmers()",
    details = "This function has been moved to the faLearn package."
  )
  .Call(`_MIC_kmers`, x, k, simplify, canonical, squeeze, anchor, clean_up, key_as_int, starting_index)
}

#' @rdname genome_to_libsvm
#' @export
genome_to_libsvm <- function(x, target_path, label = "0", k = 3L, canonical = TRUE, squeeze = FALSE, overwrite = FALSE) {
  lifecycle::deprecate_warn(
    when = "1.2.0",
    what = "genome_to_libsvm()",
    with = "faLearn::genome_to_libsvm()",
    details = "This function has been moved to the faLearn package."
  )
  .Call(`_MIC_genome_to_libsvm`, x, target_path, label, k, canonical, squeeze, overwrite)
}

#' @rdname squeezed_mers
#' @export
squeezed_mers <- function(k = 3L) {
  lifecycle::deprecate_warn(
    when = "1.2.0",
    what = "squeezed_mers()",
    with = "faLearn::squeezed_mers()",
    details = "This function has been moved to the faLearn package."
  )
  .Call(`_MIC_squeezed_mers`, k)
}

#' @rdname unsqueezed_mers
#' @export
unsqueezed_mers <- function(k = 3L) {
  lifecycle::deprecate_warn(
    when = "1.2.0",
    what = "unsqueezed_mers()",
    with = "faLearn::unsqueezed_mers()",
    details = "This function has been moved to the faLearn package."
  )
  .Call(`_MIC_unsqueezed_mers`, k)
}

#' @rdname squeezed_index_to_str
#' @export
squeezed_index_to_str <- function(x, k, starting_index = 1L) {
  lifecycle::deprecate_warn(
    when = "1.2.0",
    what = "squeezed_index_to_str()",
    with = "faLearn::squeezed_index_to_str()",
    details = "This function has been moved to the faLearn package."
  )
  .Call(`_MIC_squeezed_index_to_str`, x, k, starting_index)
}

#' @rdname unsqueezed_index_to_str
#' @export
unsqueezed_index_to_str <- function(x, k, starting_index = 1L) {
  lifecycle::deprecate_warn(
    when = "1.2.0",
    what = "unsqueezed_index_to_str()",
    with = "faLearn::unsqueezed_index_to_str()",
    details = "This function has been moved to the faLearn package."
  )
  .Call(`_MIC_unsqueezed_index_to_str`, x, k, starting_index)
}
