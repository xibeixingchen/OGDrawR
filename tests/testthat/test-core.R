test_that("spread_labels returns correct structure", {
  pos <- c(1000, 2000, 50000)
  res <- spread_labels(pos, genome_length = 100000, min_gap_bp = 4500)
  expect_type(res, "list")
  expect_length(res$pos, 3)
  expect_length(res$moved, 3)
})

test_that("spread_labels handles single gene", {
  res <- spread_labels(5000, genome_length = 100000)
  expect_equal(res$pos, 5000)
  expect_false(res$moved)
})

test_that("gene_colors has 17 entries", {
  expect_length(gene_colors, 17)
  expect_true(all(grepl("^#", gene_colors)))
})

test_that(".classify_gene works", {
  cls <- OGDrawR:::.classify_gene(c("nad1", "trnA", "rrn18", "orf100", "atp6"))
  expect_equal(cls, c("complex_I", "tRNA", "rRNA", "orf", "atp_synthase"))
})
