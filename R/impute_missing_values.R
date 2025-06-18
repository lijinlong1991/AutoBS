impute_missing_values <- function(geno, pheno, geno_id_col = 1, pheno_id_col = 1) {

  ## 1. 补全基因型数据缺失值
  if(is.data.frame(geno)) {
    # 数据框格式
    geno_imputed <- geno
    for(col in setdiff(1:ncol(geno), geno_id_col)) {
      if(any(is.na(geno[[col]]))) {
        col_mean <- mean(geno[[col]], na.rm = TRUE)
        geno_imputed[[col]][is.na(geno[[col]])] <- round(col_mean)
        message(paste("基因型列", colnames(geno)[col],
                      "中的缺失值已用平均值", round(col_mean, 2), "补全"))
      }
    }
  } else {
    # 矩阵格式
    geno_imputed <- geno
    for(col in 1:ncol(geno)) {
      if(any(is.na(geno[, col]))) {
        col_mean <- mean(geno[, col], na.rm = TRUE)
        geno_imputed[is.na(geno[, col]), col] <- round(col_mean)
        message(paste("基因型标记", colnames(geno)[col],
                      "中的缺失值已用平均值", round(col_mean, 2), "补全"))
      }
    }
  }

  ## 2. 补全表型数据缺失值
  pheno_imputed <- pheno
  for(col in setdiff(1:ncol(pheno), pheno_id_col)) {
    if(any(is.na(pheno[[col]]))) {
      col_mean <- mean(pheno[[col]], na.rm = TRUE)
      pheno_imputed[[col]][is.na(pheno[[col]])] <- col_mean
      message(paste("表型列", colnames(pheno)[col],
                    "中的缺失值已用平均值", round(col_mean, 2), "补全"))
    }
  }

  ## 3. 返回补全后的数据
  return(list(
    geno_imputed = geno_imputed,
    pheno_imputed = pheno_imputed
  ))
}
