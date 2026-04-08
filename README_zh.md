# OGDrawR <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
<!-- badges: end -->

[**English README**](README.md)

**OGDrawR** 是一个纯 R 实现的植物线粒体基因组环形图谱绘制工具，生成 [OGDraw](https://chlorobox.mpimp-golm.mpg.de/OGDraw.html) 风格的图谱，基于 [circlize](https://github.com/jokergoo/circlize) 包。

## 功能特性

- **GenBank 解析** — 从 `.gb` / `.gbk` 文件中提取基因坐标、exon/intron 边界、链向和 GC 含量
- **Intron/exon 可视化** — exon 以全高色块绘制，intron 以窄条连接
- **标签自动防重叠** — 环形坐标下迭代推开算法
- **OGDraw 标准 17 色** — 按功能分类着色，完全可自定义
- **斜体渲染** — 基因组名称和基因名以斜体显示（可直接用于发表）
- **多格式输出** — PDF、PNG、TIFF、SVG

## 安装

```r
# install.packages("remotes")
remotes::install_github("zhichengjia/OGDrawR")
```

## 快速开始

```r
library(OGDrawR)

# 一步出图
plot_mito_genome("genome.gb", "ZM4", output_file = "ZM4_mito.pdf")

# 分步使用
parsed <- parse_genbank("genome.gb", gc_window = 200)
parsed
summary(parsed)
draw_mito_map(parsed, "ZM4", output_file = "ZM4_mito.png")
```

## 自定义

```r
# 调整标签间距和字号
draw_mito_map(parsed, "ZM4",
              min_label_gap = 5000,
              label_cex = c(0.6, 0.65, 0.7))

# 自定义配色
my_colors <- gene_colors
my_colors["complex_I"] <- "#FF4500"
draw_mito_map(parsed, "ZM4", colors = my_colors)
```

## 导出函数

| 函数 | 说明 |
|---|---|
| `parse_genbank()` | 解析 GenBank 文件 → `mito_genome` 对象 |
| `draw_mito_map()` | 从 `mito_genome` 对象绘制环形图 |
| `plot_mito_genome()` | 一步完成解析 + 绘图 |
| `spread_labels()` | 标签防重叠工具函数 |
| `gene_colors` | OGDraw 17 色命名向量 |

## 引用

如果 OGDrawR 对您的研究有帮助，请引用：

> Jia Z. (2026). OGDrawR: OGDraw-Style Circular Genome Maps for Plant Mitochondria. R package version 0.1.0. https://github.com/zhichengjia/OGDrawR

## 许可证

MIT © [贾志诚](mailto:zhicheng.jia@cau.edu.cn)
