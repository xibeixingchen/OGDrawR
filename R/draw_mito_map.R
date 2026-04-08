#' Draw an OGDraw-Style Circular Mitochondrial Genome Map
#'
#' Renders a circular genome map with gene blocks (intron/exon aware),
#' GC content ring, automatic label placement, and a 17-item legend.
#'
#' @param parsed A `mito_genome` object returned by [parse_genbank()].
#' @param genome_name Character string shown in the centre of the map
#'   (rendered in italic).
#' @param output_file Optional output path. Extension determines format:
#'   `.pdf`, `.png`, `.tiff`, or `.svg`. If `NULL`, draws to the current
#'   device.
#' @param w,h Width and height of the output in inches. Default `12`.
#' @param res Resolution in dpi for raster formats. Default `600`.
#' @param colors Named character vector of gene-category colours.
#'   Default [gene_colors].
#' @param min_label_gap Minimum gap in bp between labels. Default `4500`.
#' @param label_cex Numeric vector of length 3 giving label sizes for
#'   tRNA, rRNA, and other genes. Default `c(0.65, 0.70, 0.73)`.
#'
#' @return Invisible `NULL`. Called for its side effect of producing a plot.
#'
#' @examples
#' \dontrun{
#' parsed <- parse_genbank("genome.gb")
#' draw_mito_map(parsed, "ZM4", output_file = "map.pdf")
#' }
#' @export
#' @importFrom circlize circos.clear circos.par circos.initialize circos.track
#'   circos.rect circos.lines circos.text
#' @importFrom graphics par text legend
#' @importFrom grDevices pdf png tiff dev.off
#' @importFrom tools file_ext
draw_mito_map <- function(parsed, genome_name = "Genome",
                          output_file = NULL, w = 12, h = 12, res = 600,
                          colors = gene_colors,
                          min_label_gap = 4500,
                          label_cex = c(0.65, 0.70, 0.73)) {

  genes        <- parsed$genes
  exons_list   <- parsed$exons
  gc_data      <- parsed$gc
  genome_length <- parsed$genome_length

  # -- open device -----------------------------------------------------------
  if (!is.null(output_file)) {
    ext <- file_ext(output_file)
    switch(ext,
      pdf  = pdf(output_file, width = w, height = h),
      png  = png(output_file, width = w, height = h, units = "in", res = res),
      tiff = tiff(output_file, width = w, height = h, units = "in",
                  res = res, compression = "lzw"),
      svg  = svg(output_file, width = w, height = h),
      stop("Unsupported format: ", ext, call. = FALSE)
    )
    on.exit(dev.off(), add = TRUE)
  }

  par(mar = c(3, 0, 0, 0))
  circos.clear()
  circos.par(
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0, 0),
    start.degree = 90, gap.degree = 0,
    canvas.xlim  = c(-1.4, 1.4),
    canvas.ylim  = c(-1.5, 1.3)
  )
  circos.initialize(factors = "genome", xlim = c(0, genome_length))

  plus_genes  <- genes[genes$strand ==  1, ]
  minus_genes <- genes[genes$strand == -1, ]

  plus_spread  <- spread_labels(plus_genes$mid,  genome_length, min_label_gap)
  minus_spread <- spread_labels(minus_genes$mid, genome_length, min_label_gap)

  cex_tRNA  <- label_cex[1]
  cex_rRNA  <- label_cex[2]
  cex_other <- label_cex[3]

  # -- Track 1: plus-strand labels -------------------------------------------
  circos.track(
    factors = "genome", ylim = c(0, 1),
    track.height = 0.20, bg.border = NA,
    panel.fun = function(x, y) {
      for (i in seq_len(nrow(plus_genes))) {
        g       <- plus_genes[i, ]
        gene_x  <- g$mid
        label_x <- plus_spread$pos[i]
        moved   <- plus_spread$moved[i]
        cv <- switch(g$category, tRNA = cex_tRNA, rRNA = cex_rRNA, cex_other)
        if (moved) {
          circos.lines(c(gene_x, gene_x), c(0.0, 0.15),
                       col = "grey40", lwd = 0.35)
          circos.lines(c(gene_x, label_x), c(0.15, 0.25),
                       col = "grey40", lwd = 0.35, straight = TRUE)
        } else {
          circos.lines(c(gene_x, gene_x), c(0.0, 0.25),
                       col = "grey40", lwd = 0.35)
          label_x <- gene_x
        }
        circos.text(label_x, 0.28, g$name,
                    facing = "clockwise", niceFacing = TRUE,
                    adj = c(0, 0.5), cex = cv, font = 4)
      }
    }
  )

  # -- Track 2: plus-strand gene blocks (exon/intron) -------------------------
  circos.track(
    factors = "genome", ylim = c(0, 1),
    track.height = 0.06, bg.border = NA,
    panel.fun = function(x, y) {
      for (i in seq_len(nrow(plus_genes))) {
        g <- plus_genes[i, ]
        col <- colors[g$category]
        if (is.na(col)) col <- colors["other"]
        .draw_gene_blocks(g, exons_list[[g$name]], col, colors["intron"])
      }
    }
  )

  # -- Track 3: backbone line ------------------------------------------------
  circos.track(
    factors = "genome", ylim = c(0, 1),
    track.height = 0.006, bg.border = NA, bg.col = "black",
    panel.fun = function(x, y) {}
  )

  # -- Track 4: minus-strand gene blocks (exon/intron) -----------------------
  circos.track(
    factors = "genome", ylim = c(0, 1),
    track.height = 0.06, bg.border = NA,
    panel.fun = function(x, y) {
      for (i in seq_len(nrow(minus_genes))) {
        g <- minus_genes[i, ]
        col <- colors[g$category]
        if (is.na(col)) col <- colors["other"]
        .draw_gene_blocks(g, exons_list[[g$name]], col, colors["intron"])
      }
    }
  )

  # -- Track 5: minus-strand labels ------------------------------------------
  circos.track(
    factors = "genome", ylim = c(0, 1),
    track.height = 0.15, bg.border = NA,
    panel.fun = function(x, y) {
      for (i in seq_len(nrow(minus_genes))) {
        g       <- minus_genes[i, ]
        gene_x  <- g$mid
        label_x <- minus_spread$pos[i]
        moved   <- minus_spread$moved[i]
        cv <- switch(g$category, tRNA = cex_tRNA, rRNA = cex_rRNA, cex_other)
        if (moved) {
          circos.lines(c(gene_x, gene_x), c(1.0, 0.85),
                       col = "grey40", lwd = 0.35)
          circos.lines(c(gene_x, label_x), c(0.85, 0.75),
                       col = "grey40", lwd = 0.35, straight = TRUE)
        } else {
          circos.lines(c(gene_x, gene_x), c(1.0, 0.75),
                       col = "grey40", lwd = 0.35)
          label_x <- gene_x
        }
        circos.text(label_x, 0.70, g$name,
                    facing = "clockwise", niceFacing = TRUE,
                    adj = c(1, 0.5), cex = cv, font = 4)
      }
    }
  )

  # -- Spacer ----------------------------------------------------------------
  circos.track(
    factors = "genome", ylim = c(0, 1),
    track.height = 0.06, bg.border = NA,
    panel.fun = function(x, y) {}
  )

  # -- Track 6: GC content ring ----------------------------------------------
  gc_mean <- mean(gc_data$gc)
  gc_min  <- min(gc_data$gc) - 0.01
  gc_max  <- max(gc_data$gc) + 0.01

  circos.track(
    factors = "genome", ylim = c(gc_min, gc_max),
    track.height = 0.08, bg.border = NA,
    panel.fun = function(x, y) {
      circos.rect(0, gc_min, genome_length, gc_max,
                  col = "#E0E0E0", border = NA)
      for (i in seq_len(nrow(gc_data))) {
        s <- gc_data$start[i]; e <- gc_data$end[i]; gv <- gc_data$gc[i]
        bottom <- gc_max - (gv - gc_min)
        circos.rect(s, bottom, e, gc_max, col = "#9A9A9A", border = NA)
      }
      mean_y <- gc_max - (gc_mean - gc_min)
      circos.lines(c(0, genome_length), c(mean_y, mean_y),
                   col = "#444444", lwd = 0.3)
    }
  )

  # -- Centre text -----------------------------------------------------------
  text(0,  0.04, bquote(italic(.(genome_name))), cex = 1.6, font = 1)
  text(0, -0.03, "mitochondrial genome", cex = 1.1)
  text(0, -0.09, paste0(format(genome_length, big.mark = ","), " bp"), cex = 1.1)

  # -- Legend ----------------------------------------------------------------
  legend_labels <- c(
    "complex I (NADH dehydrogenase)",
    "complex II (succinate dehydrogenase)",
    "complex III (ubichinol cytochrome c reductase)",
    "complex IV (cytochrome c oxidase)",
    "ATP synthase",
    "cytochrome c biogenesis",
    "RNA polymerase",
    "ribosomal proteins (SSU)",
    "ribosomal proteins (LSU)",
    "maturases",
    "other genes",
    "ORFs",
    "transfer RNAs",
    "ribosomal RNAs",
    "origin of replication",
    "polycistronic transcripts",
    "introns"
  )
  legend_cols <- c(
    colors["complex_I"], colors["complex_II"], colors["complex_III"],
    colors["complex_IV"], colors["atp_synthase"], colors["cytochrome_c"],
    colors["rna_pol"], colors["ribo_SSU"], colors["ribo_LSU"],
    colors["maturase"], colors["other"], colors["orf"],
    colors["tRNA"], colors["rRNA"], colors["ori_rep"],
    colors["poly_trans"], colors["intron"]
  )
  legend(x = -1.4, y = -0.72,
         legend = legend_labels, fill = legend_cols,
         border = "black", cex = 0.9, bty = "n",
         y.intersp = 1.05, x.intersp = 0.5)

  circos.clear()
  invisible(NULL)
}


# -- internal: draw exon/intron blocks for a single gene ----------------------
#' @noRd
.draw_gene_blocks <- function(g, exon_df, gene_col, intron_col,
                              exon_y0 = 0.0, exon_y1 = 1.0,
                              intron_y0 = 0.3, intron_y1 = 0.7) {
  if (!g$has_intron || is.null(exon_df) || nrow(exon_df) < 2) {
    circos.rect(g$start, exon_y0, g$end, exon_y1,
                col = gene_col, border = "black", lwd = 0.3)
  } else {
    circos.rect(g$start, intron_y0, g$end, intron_y1,
                col = intron_col, border = "black", lwd = 0.2)
    for (j in seq_len(nrow(exon_df))) {
      circos.rect(exon_df$start[j], exon_y0, exon_df$end[j], exon_y1,
                  col = gene_col, border = "black", lwd = 0.3)
    }
  }
}
