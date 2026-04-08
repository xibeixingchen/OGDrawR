#' One-Step Mitochondrial Genome Map
#'
#' Parse a GenBank file and draw the circular map in a single call.
#'
#' @param gb_file Path to GenBank file.
#' @param genome_name Character label for the centre. If `NULL`, derived from
#'   the file name.
#' @param output_file Optional output path (pdf/png/tiff/svg).
#' @param gc_window GC window size in bp. Default `200`.
#' @param ... Additional arguments passed to [draw_mito_map()].
#'
#' @return Invisible `mito_genome` object.
#'
#' @examples
#' \dontrun{
#' plot_mito_genome("genome.gb", "ZM4", output_file = "ZM4_mito.pdf")
#' }
#' @export
plot_mito_genome <- function(gb_file, genome_name = NULL,
                             output_file = NULL, gc_window = 200, ...) {
  parsed <- parse_genbank(gb_file, gc_window = gc_window)
  if (is.null(genome_name)) {
    genome_name <- tools::file_path_sans_ext(basename(gb_file))
  }
  draw_mito_map(parsed, genome_name = genome_name,
                output_file = output_file, ...)
  invisible(parsed)
}
