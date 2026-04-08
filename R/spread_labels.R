#' Spread Labels to Avoid Overlapping on a Circular Layout
#'
#' Iteratively pushes apart label positions that are closer than `min_gap_bp`,
#' wrapping around the circular genome boundary.
#'
#' @param gene_mids Numeric vector of midpoint positions (bp) for each gene.
#' @param genome_length Total genome length in bp.
#' @param min_gap_bp Minimum gap in bp between adjacent labels. Default `4500`.
#'
#' @return A list with:
#' \describe{
#'   \item{pos}{Adjusted label positions (same order as input).}
#'   \item{moved}{Logical vector indicating which labels were displaced
#'     significantly (> 800 bp).}
#' }
#' @export
spread_labels <- function(gene_mids, genome_length, min_gap_bp = 4500) {
  n <- length(gene_mids)
  if (n <= 1) return(list(pos = gene_mids, moved = rep(FALSE, n)))
  ord  <- order(gene_mids)
  adj  <- gene_mids[ord]
  orig <- adj
  for (iter in 1:20) {
    moved <- FALSE
    for (i in 2:n) {
      gap <- adj[i] - adj[i - 1]
      if (gap < min_gap_bp) {
        shift <- (min_gap_bp - gap) / 2
        adj[i - 1] <- adj[i - 1] - shift
        adj[i]     <- adj[i] + shift
        moved <- TRUE
      }
    }
    gap_wrap <- (adj[1] + genome_length) - adj[n]
    if (gap_wrap < min_gap_bp) {
      shift <- (min_gap_bp - gap_wrap) / 2
      adj[n] <- adj[n] - shift
      adj[1] <- adj[1] + shift
    }
    if (!moved) break
  }
  adj <- adj %% genome_length
  was_moved <- abs(adj - orig) > 800
  result_pos   <- numeric(n); result_pos[ord]   <- adj
  result_moved <- logical(n); result_moved[ord] <- was_moved
  list(pos = result_pos, moved = result_moved)
}
