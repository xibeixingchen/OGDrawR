# OGDrawR 示例脚本 - 绘制拟南芥线粒体基因组图

library(OGDrawR)

# 获取示例数据路径
gb_file <- system.file("extdata", "Arabidopsis_thaliana.gb", package = "OGDrawR")

# 绘制基因组图
plot_mito_genome(
  input_file = gb_file,
  output_file = "Arabidopsis_mito_map.png",
  title = "Arabidopsis thaliana"
)

cat("基因组图已保存到: Arabidopsis_mito_map.png\n")
