# OGDrawR <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
<!-- badges: end -->

[**中文文档 / Chinese README**](README_zh.md)

**OGDrawR** generates [OGDraw](https://chlorobox.mpimp-golm.mpg.de/OGDraw.html)-style circular genome maps for plant mitochondria entirely in R, powered by [circlize](https://github.com/jokergoo/circlize).

## Features

- **GenBank parsing** — extracts gene coordinates, exon/intron boundaries, strand, and GC content from `.gb` / `.gbk` files
- **Intron/exon visualization** — exons rendered as full-height blocks, introns as narrow connecting bars
- **Automatic label layout** — iterative de-overlapping algorithm on circular coordinates
- **OGDraw 17-color palette** — standard functional-category colours, fully customisable
- **Italic rendering** — genome name and gene names displayed in italic (publication-ready)
- **Multiple output formats** — PDF, PNG, TIFF, SVG

## Installation

```r
# install.packages("remotes")
remotes::install_github("zhichengjia/OGDrawR")
```

## Quick Start

```r
library(OGDrawR)

# One step
plot_mito_genome("genome.gb", "ZM4", output_file = "ZM4_mito.pdf")

# Or step by step
parsed <- parse_genbank("genome.gb", gc_window = 200)
parsed
summary(parsed)
draw_mito_map(parsed, "ZM4", output_file = "ZM4_mito.png")
```

## Customisation

```r
# Adjust label spacing and font size
draw_mito_map(parsed, "ZM4",
              min_label_gap = 5000,
              label_cex = c(0.6, 0.65, 0.7))

# Custom colours
my_colors <- gene_colors
my_colors["complex_I"] <- "#FF4500"
draw_mito_map(parsed, "ZM4", colors = my_colors)
```

## Exported Functions

| Function | Description |
|---|---|
| `parse_genbank()` | Parse a GenBank file → `mito_genome` object |
| `draw_mito_map()` | Draw circular map from a `mito_genome` object |
| `plot_mito_genome()` | One-step parse + draw |
| `spread_labels()` | Label de-overlapping utility |
| `gene_colors` | OGDraw 17-colour palette (named vector) |

## Citation

If you use OGDrawR in your research, please cite:

> Jia Z. (2026). OGDrawR: OGDraw-Style Circular Genome Maps for Plant Mitochondria. R package version 0.1.0. https://github.com/zhichengjia/OGDrawR

## License

MIT © [Zhicheng Jia](mailto:zhicheng.jia@cau.edu.cn)
