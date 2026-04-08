#' OGDraw Standard Color Palette for Mitochondrial Gene Categories
#'
#' A named character vector of 17 colors following the OGDraw convention.
#'
#' @format Named character vector with hex color codes. Names correspond to
#'   gene categories: `complex_I`, `complex_II`, `complex_III`, `complex_IV`,
#'   `atp_synthase`, `cytochrome_c`, `rna_pol`, `ribo_SSU`, `ribo_LSU`,
#'   `maturase`, `other`, `orf`, `tRNA`, `rRNA`, `ori_rep`, `poly_trans`,
#'   `intron`.
#' @export
gene_colors <- c(
  complex_I    = "#FFD700",
  complex_II   = "#32CD32",
  complex_III  = "#FFDEAD",
  complex_IV   = "#FFB6C1",
  atp_synthase = "#9ACD32",
  cytochrome_c = "#228B22",
  rna_pol      = "#B22222",
  ribo_SSU     = "#F5DEB3",
  ribo_LSU     = "#D2B48C",
  maturase     = "#FF8C00",
  other        = "#DA70D6",
  orf          = "#AFEEEE",
  tRNA         = "#00008B",
  rRNA         = "#FF0000",
  ori_rep      = "#FFB6C1",
  poly_trans   = "#BA55D3",
  intron       = "#dad7d7ff"
)
