check_data <- function(geno, pheno, geno_id_col = 1, pheno_id_col = 1,
                       na_strings = c("NA", "", ".", "-")) {
  # 初始化结果列表
  result <- list()

  ## 1. 转换可能的缺失值表示为真正的NA
  convert_na <- function(x, na_strings) {
    if(is.character(x) || is.factor(x)) {
      x[x %in% na_strings] <- NA
    }
    x
  }

  if(is.data.frame(geno)) {
    geno[-geno_id_col] <- lapply(geno[-geno_id_col], convert_na, na_strings)
  } else {
    geno[geno %in% na_strings] <- NA
  }

  pheno[-pheno_id_col] <- lapply(pheno[-pheno_id_col], convert_na, na_strings)

  ## 2. 检查ID匹配情况（同前）
  geno_ids <- if(is.null(dim(geno))) rownames(geno) else geno[[geno_id_col]]
  pheno_ids <- pheno[[pheno_id_col]]

  common_ids <- intersect(geno_ids, pheno_ids)
  only_geno <- setdiff(geno_ids, pheno_ids)
  only_pheno <- setdiff(pheno_ids, geno_ids)

  result$ID_match <- list(
    common_samples = length(common_ids),
    only_in_geno = length(only_geno),
    only_in_pheno = length(only_pheno)
  )

  ## 3. 检查geno文件结构并计算缺失值
  if(is.null(dim(geno))) {
    # 矩阵格式
    geno_na <- sum(is.na(geno))
    marker_na <- colSums(is.na(geno))
    sample_na <- rowSums(is.na(geno))

    result$geno <- list(
      format = "matrix",
      samples = nrow(geno),
      markers = ncol(geno),
      total_na = geno_na,
      na_rate = geno_na / length(geno),
      marker_names = colnames(geno),
      sample_ids = rownames(geno)
    )
  } else {
    # 数据框格式
    geno_data <- geno[-geno_id_col]
    geno_na <- sum(is.na(geno_data))
    marker_na <- sapply(geno_data, function(x) sum(is.na(x)))
    sample_na <- apply(geno_data, 1, function(x) sum(is.na(x)))

    result$geno <- list(
      format = "data.frame",
      samples = nrow(geno),
      markers = ncol(geno) - 1,
      total_na = geno_na,
      na_rate = geno_na / length(unlist(geno_data)),
      marker_names = colnames(geno)[-geno_id_col],
      sample_ids = geno[[geno_id_col]]
    )
  }

  ## 4. 检查pheno文件结构并计算缺失值
  pheno_data <- pheno[-pheno_id_col]
  pheno_na <- sum(is.na(pheno_data))
  trait_na <- sapply(pheno_data, function(x) sum(is.na(x)))
  sample_na_pheno <- apply(pheno_data, 1, function(x) sum(is.na(x)))

  result$pheno <- list(
    samples = nrow(pheno),
    traits = ncol(pheno) - 1,
    total_na = pheno_na,
    na_rate = pheno_na / length(unlist(pheno_data)),
    trait_names = colnames(pheno)[-pheno_id_col],
    sample_ids = pheno[[pheno_id_col]]
  )

  ## 5. 输出总结信息
  cat("=== 数据结构检查结果（含缺失值统计） ===\n\n")

  # 基因型数据信息
  cat("基因型数据(geno):\n")
  cat("- 格式:", result$geno$format, "\n")
  cat("- 样本量:", result$geno$samples, "\n")
  cat("- 分子标记数量:", result$geno$markers, "\n")
  cat("- 总缺失值数:", result$geno$total_na, "\n")
  cat("- 缺失率:", round(result$geno$na_rate * 100, 2), "%\n")
  cat("\n")

  # 表型数据信息
  cat("表型数据(pheno):\n")
  cat("- 样本量:", result$pheno$samples, "\n")
  cat("- 表型数量:", result$pheno$traits, "\n")
  cat("- 表型名称:", result$pheno$trait_names, "\n")
  cat("- 总缺失值数:", result$pheno$total_na, "\n")
  cat("- 缺失率:", round(result$pheno$na_rate * 100, 2), "%\n")
  cat("\n")

  # ID匹配情况
  cat("ID匹配情况:\n")
  cat("- 共同样本数:", result$ID_match$common_samples, "\n")
  cat("- 仅存在于geno中的样本数:", result$ID_match$only_in_geno, "\n")
  cat("- 仅存在于pheno中的样本数:", result$ID_match$only_in_pheno, "\n\n")

  ## 7. 返回完整结果
  invisible(result)
}
