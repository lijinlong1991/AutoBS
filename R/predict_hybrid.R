predict_hybrid <- function(genotype, pheno,parent_pairs) {
  # 检查输入数据
  if (ncol(parent_pairs) != 2) {
    stop("parent_pairs should have exactly two columns (P1 and P2)")
  }
  if (nrow(genotype) != length(pheno)) {
    stop("Number of genotypes does not match number of phenotypes")
  }
  if (!all(c(parent_pairs[,1], parent_pairs[,2]) %in% rownames(genotype))) {
    stop("Some parent IDs in parent_pairs are not found in genotype data")
  }

  # 加载必要包
  if (!requireNamespace("rrBLUP", quietly = TRUE)) {
    stop("Package 'rrBLUP' needed for this function to work. Please install it.")
  }
  #parent_pairs=hybrid_list;genotype=geno_new[,-1];pheno=pheno$HEI
  # 计算原始G矩阵
  #X <- scale(genotype, center = TRUE, scale = TRUE)
  #G <- tcrossprod(X) / ncol(X)
  G <- rrBLUP::A.mat(genotype - 1)  # 输入需为等位基因计数（0/1/2）减去1
  rownames(G) <- colnames(G) <- row.names(genotype)
  n_individuals <- nrow(genotype)
  n_hybrids <- nrow(parent_pairs)

  # 扩展G矩阵（原始个体 + 杂交后代）
  G_hybrid <- matrix(0, nrow = n_individuals + n_hybrids, ncol = n_individuals + n_hybrids)
  G_hybrid[1:n_individuals, 1:n_individuals] <- G

  # 获取亲本在原始G矩阵中的索引
  parent_indices <- matrix(match(as.matrix(parent_pairs), rownames(genotype)),
                           ncol = 2, byrow = FALSE)

  for (k in 1:n_hybrids) {
    p1 <- parent_indices[k, 1]
    p2 <- parent_indices[k, 2]

    # 后代与原始个体的关系
    G_hybrid[n_individuals + k, 1:n_individuals] <- 0.5 * (G[p1, ] + G[p2, ])
    G_hybrid[1:n_individuals, n_individuals + k] <- G_hybrid[n_individuals + k, 1:n_individuals]

    # 后代自身的近交系数
    G_hybrid[n_individuals + k, n_individuals + k] <- 1 + 0.5 * G[p1, p2]

    # 后代之间的关系
    for (m in 1:k) {
      p1_m <- parent_indices[m, 1]
      p2_m <- parent_indices[m, 2]
      G_hybrid[n_individuals + k, n_individuals + m] <- 0.25 * (G[p1, p1_m] + G[p1, p2_m] + G[p2, p1_m] + G[p2, p2_m])
      G_hybrid[n_individuals + m, n_individuals + k] <- G_hybrid[n_individuals + k, n_individuals + m]
    }
  }

  # 准备数据（杂交后代表型设为NA）
  pheno_combined <- c(pheno, rep(NA, n_hybrids))
  data <- data.frame(
    ID = c(rownames(genotype),
           paste0(parent_pairs[,1], "_x_", parent_pairs[,2])),
    y = pheno_combined
  )

  # 将G矩阵转换为rrBLUP所需的kinship矩阵
  kinship <- G_hybrid

  # 拟合GBLUP模型
  fit <- rrBLUP::mixed.solve(
    y = data$y,
    K = kinship,      # 亲缘关系矩阵
    method = "REML"
  )

  # 提取育种值
  breeding_values <- fit$u +  as.vector(fit$beta)
  hybrid_pred <- tail(breeding_values, n_hybrids)
  names(hybrid_pred) <- paste0(parent_pairs[,1], "_x_", parent_pairs[,2])

  # 返回结果
  list(
    predicted_values = hybrid_pred,
    G_matrix = G_hybrid,
    model_fit = fit
  )
}
