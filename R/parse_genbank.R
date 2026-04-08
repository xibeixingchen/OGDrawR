#' Parse a GenBank File for Mitochondrial Genome Annotation
#'
#' Extracts gene coordinates (with exon/intron boundaries), strand information,
#' functional categories, and sliding-window GC content from a GenBank flat file.
#'
#' @param gb_file Path to a GenBank format file (`.gb` or `.gbk`).
#' @param gc_window Integer. Window size in bp for GC content calculation.
#'   Default `200`.
#'
#' @return An object of class `mito_genome`, a list containing:
#' \describe{
#'   \item{genes}{A `data.frame` with columns `name`, `start`, `end`, `strand`,
#'     `has_intron`, `category`, `mid`, `len`.}
#'   \item{exons}{A named list of `data.frame`s (one per gene) each with
#'     columns `start` and `end` giving exon coordinates.}
#'   \item{gc}{A `data.frame` with columns `start`, `end`, `mid`, `gc`.}
#'   \item{genome_length}{Integer. Total genome length in bp.}
#' }
#'
#' @examples
#' \dontrun{
#' parsed <- parse_genbank("path/to/genome.gb")
#' print(parsed)
#' }
#' @export
parse_genbank <- function(gb_file, gc_window = 200) {

  if (!file.exists(gb_file)) {
    stop("GenBank file not found: ", gb_file, call. = FALSE)
  }

  raw <- readLines(gb_file)
  content <- paste(raw, collapse = "\n")

 # -- gene annotations -------------------------------------------------------
  blocks <- strsplit(content, "\n     gene ")[[1]]
  genes <- list()
  for (block in blocks[-1]) {
    lines <- strsplit(block, "\n")[[1]]
    loc   <- trimws(lines[1])
    strand <- 1L
    if (grepl("complement", loc)) {
      strand <- -1L
      loc <- gsub("complement\\(|\\)$", "", loc)
    }
    has_intron <- grepl("join", loc)

    exons <- list()
    if (has_intron) {
      loc_inner <- gsub("join\\(|\\)$", "", loc)
      parts <- strsplit(loc_inner, ",")[[1]]
      for (p in parts) {
        nums <- as.numeric(regmatches(p, gregexpr("\\d+", p))[[1]])
        if (length(nums) >= 2) {
          exons[[length(exons) + 1]] <- c(start = nums[1], end = nums[2])
        }
      }
      start <- exons[[1]]["start"]
      end   <- exons[[length(exons)]]["end"]
    } else {
      nums <- as.numeric(regmatches(loc, gregexpr("\\d+", loc))[[1]])
      if (length(nums) < 2) next
      start <- nums[1]; end <- nums[2]
      exons[[1]] <- c(start = start, end = end)
    }

    name <- NA
    for (l in lines[2:min(8, length(lines))]) {
      m <- regmatches(l, regexec('/gene="([^"]+)"', l))[[1]]
      if (length(m) == 2) { name <- m[2]; break }
    }
    if (is.na(name)) next

    exon_df <- do.call(rbind, lapply(exons, function(e) {
      data.frame(start = e["start"], end = e["end"])
    }))

    genes[[length(genes) + 1]] <- list(
      name = name, start = as.numeric(start), end = as.numeric(end),
      strand = strand, has_intron = has_intron, exons = exon_df
    )
  }

  exons_list <- lapply(genes, function(g) g$exons)
  df <- do.call(rbind, lapply(genes, function(g) {
    data.frame(name = g$name, start = g$start, end = g$end,
               strand = g$strand, has_intron = g$has_intron,
               stringsAsFactors = FALSE)
  }))
  dup <- duplicated(df$name)
  df <- df[!dup, ]
  exons_list <- exons_list[!dup]
  names(exons_list) <- df$name

  df$category <- .classify_gene(df$name)
  df$mid <- (df$start + df$end) / 2
  df$len <- abs(df$end - df$start)

  # -- sequence / GC ----------------------------------------------------------
  in_seq <- FALSE; seq_parts <- c()
  for (line in raw) {
    line <- trimws(line)
    if (grepl("^ORIGIN", line)) { in_seq <- TRUE; next }
    if (grepl("^//", line))     break
    if (in_seq) seq_parts <- c(seq_parts, gsub("[^a-zA-Z]", "", line))
  }
  sequence <- toupper(paste(seq_parts, collapse = ""))
  genome_length <- nchar(sequence)

  n_win <- floor(genome_length / gc_window)
  gc_df <- data.frame(start = integer(n_win), end = integer(n_win),
                      mid = numeric(n_win), gc = numeric(n_win))
  for (i in seq_len(n_win)) {
    s <- (i - 1) * gc_window + 1
    e <- min(i * gc_window, genome_length)
    w <- substr(sequence, s, e); wl <- nchar(w)
    gc_df$start[i] <- s; gc_df$end[i] <- e; gc_df$mid[i] <- (s + e) / 2
    gc_df$gc[i] <- (nchar(gsub("[^G]", "", w)) + nchar(gsub("[^C]", "", w))) / wl
  }

  structure(
    list(genes = df, exons = exons_list, gc = gc_df,
         genome_length = genome_length),
    class = "mito_genome"
  )
}

# -- S3 methods ---------------------------------------------------------------

#' @export
print.mito_genome <- function(x, ...) {
  cat("Mitochondrial genome:",
      format(x$genome_length, big.mark = ","), "bp\n")
  cat("Genes:", nrow(x$genes), "\n")
  n_intron <- sum(x$genes$has_intron)
  if (n_intron > 0) cat("Genes with introns:", n_intron, "\n")
  cat("Mean GC:", round(mean(x$gc$gc), 4), "\n")
  invisible(x)
}

#' @export
summary.mito_genome <- function(object, ...) {
  cat("Mitochondrial genome:",
      format(object$genome_length, big.mark = ","), "bp\n")
  cat("Genes:", nrow(object$genes), "\n\n")
  cat("Category breakdown:\n")
  tbl <- table(object$genes$category)
  for (nm in sort(names(tbl))) {
    cat("  ", nm, ": ", tbl[nm], "\n", sep = "")
  }
  intron_genes <- object$genes[object$genes$has_intron, ]
  if (nrow(intron_genes) > 0) {
    cat("\nGenes with introns:\n")
    for (nm in intron_genes$name) {
      ex <- object$exons[[nm]]
      cat("  ", nm, " (", nrow(ex), " exons)\n", sep = "")
    }
  }
  invisible(object)
}

# -- internal helper ----------------------------------------------------------

#' Classify gene name into functional category
#' @noRd
.classify_gene <- function(names) {
  sapply(names, function(n) {
    nl <- tolower(n)
    if (grepl("^trn",  nl)) return("tRNA")
    if (grepl("^rrn",  nl)) return("rRNA")
    if (grepl("^nad",  nl)) return("complex_I")
    if (grepl("^sdh",  nl)) return("complex_II")
    if (grepl("^cob",  nl)) return("complex_III")
    if (grepl("^cox",  nl)) return("complex_IV")
    if (grepl("^atp",  nl)) return("atp_synthase")
    if (grepl("^ccm",  nl)) return("cytochrome_c")
    if (grepl("^rps",  nl)) return("ribo_SSU")
    if (grepl("^rpl",  nl)) return("ribo_LSU")
    if (grepl("^mat",  nl)) return("maturase")
    if (grepl("^orf",  nl)) return("orf")
    if (grepl("^rpo",  nl)) return("rna_pol")
    "other"
  }, USE.NAMES = FALSE)
}
