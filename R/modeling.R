modeling <- function(geno, pheno, k = 5, plot = TRUE, seed = 123, verbose = TRUE) {
  # 检查并安装必要包
  required_pkgs <- c("rrBLUP", "caret")
  if (plot) required_pkgs <- c(required_pkgs, "ggplot2")

  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }

  # 加载包
  suppressPackageStartupMessages({
    library(rrBLUP)
    library(caret)
    if (plot) library(ggplot2)
  })

  # 数据准备和验证
  if (is.data.frame(pheno) || is.matrix(pheno)) {
    if (ncol(pheno) != 1) stop("pheno应为单列数据框或向量")
    y <- pheno[, 1]
    names(y) <- rownames(pheno)
  } else {
    y <- pheno
    if (is.null(names(y))) names(y) <- 1:length(y)
  }

  # 检查geno格式并提取ID
  if (is.data.frame(geno)) {
    geno_ids <- geno[, 1]
    geno_for_G <- as.matrix(geno[, -1, drop = FALSE])
    rownames(geno_for_G) <- geno_ids
  } else if (is.matrix(geno)) {
    geno_ids <- rownames(geno)
    if (is.null(geno_ids)) {
      geno_ids <- 1:nrow(geno)
      rownames(geno) <- geno_ids
    }
    geno_for_G <- geno
  } else {
    stop("geno应为数据框或矩阵")
  }

  # 检查geno矩阵是否有效
  if (any(is.na(geno_for_G))) {
    warning("基因型数据包含缺失值，将在计算G矩阵前用均值填补")
    geno_for_G <- apply(geno_for_G, 2, function(x) {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      return(x)
    })
  }

  # 检查数据范围并转换为0-2格式
  if (max(geno_for_G) > 2 | min(geno_for_G) < 0) {
    stop("基因型数据应为等位基因计数（0/1/2格式）")
  }

  # 按pheno的顺序排列geno数据
  if (!all(names(y) %in% rownames(geno_for_G))) {
    stop("表型数据与基因型数据的样本ID不匹配")
  }
  geno_for_G <- geno_for_G[match(names(y), rownames(geno_for_G)), , drop = FALSE]

  # 计算G矩阵（使用rrBLUP的A.mat函数）
  if (verbose) message("计算基因组关系矩阵(G矩阵)...")
  tryCatch({
    # 确保输入是数值矩阵
    geno_numeric <- apply(geno_for_G, 2, as.numeric)
    rownames(geno_numeric) <- rownames(geno_for_G)
    G <- A.mat(geno_numeric - 1)  # 输入需为等位基因计数（0/1/2）减去1
    rownames(G) <- colnames(G) <- names(y)
  }, error = function(e) {
    stop(paste("计算G矩阵出错:", e$message,
               "\n请检查基因型数据是否为数值矩阵，样本在行，标记在列"))
  })

  # 设置随机种子
  set.seed(seed)

  # 移除缺失值
  complete_cases <- !is.na(y)
  G <- G[complete_cases, complete_cases, drop = FALSE]
  y <- y[complete_cases]

  if (verbose) {
    message(paste("开始", k, "折交叉验证，有效样本量:", length(y)))
  }

  # 创建交叉验证分组
  folds <- createFolds(y, k = k, list = TRUE, returnTrain = FALSE)

  # 初始化存储
  predictions <- numeric(length(y))
  names(predictions) <- names(y)

  # 存储评估指标
  eval_metrics <- data.frame(
    Fold = 1:k,
    Cor = numeric(k),
    MAE = numeric(k),
    MSE = numeric(k),
    RMSE = numeric(k),
    stringsAsFactors = FALSE
  )

  # 交叉验证循环
  for (i in 1:k) {
    if (verbose) message(paste("处理第", i, "折..."))

    test_indices <- folds[[i]]
    train_indices <- setdiff(1:length(y), test_indices)

    # 训练模型
    fit <- mixed.solve(
      y = y[train_indices],
      K = G[train_indices, train_indices, drop = FALSE]
    )

    # 计算测试集的预测值
    G_train_test <- G[train_indices, test_indices, drop = FALSE]
    pred_test <- t(G_train_test) %*% fit$u

    # 存储预测结果
    predictions[test_indices] <- pred_test

    # 计算评估指标
    actual <- y[test_indices]
    eval_metrics$Cor[i] <- cor(pred_test, actual, use = "complete.obs")
    eval_metrics$MAE[i] <- mean(abs(pred_test - actual))
    eval_metrics$MSE[i] <- mean((pred_test - actual)^2)
    eval_metrics$RMSE[i] <- sqrt(eval_metrics$MSE[i])
  }

  # 计算整体评估指标
  overall_metrics <- data.frame(
    Metric = c("Correlation", "MAE", "MSE", "RMSE"),
    Value = c(
      cor(predictions, y, use = "complete.obs"),
      mean(abs(predictions - y)),
      mean((predictions - y)^2),
      sqrt(mean((predictions - y)^2))
    ),
    stringsAsFactors = FALSE
  )

  # 创建图形
  plot_obj <- NULL
  if (plot) {
    plot_df <- data.frame(
      Observed = y,
      Predicted = predictions,
      ID = names(y)
    )

    plot_obj <- ggplot(plot_df, aes(x = Observed, y = Predicted)) +
      geom_point(alpha = 0.6, color = "blue") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
      labs(title = paste("GBLUP", k, "折交叉验证结果"),
           subtitle = sprintf("整体相关性: %.3f, RMSE: %.3f",
                              overall_metrics$Value[1],
                              overall_metrics$Value[4]),
           x = "观测值",
           y = "预测值") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  }

  # 返回结果
  result <- list(
    predictions = predictions,
    fold_metrics = eval_metrics,
    overall_metrics = overall_metrics,
    plot = plot_obj,
    G_matrix = G,
    parameters = list(
      k = k,
      seed = seed,
      n_samples = length(y),
      na_removed = sum(!complete_cases),
      geno_dim = dim(geno_for_G),
      pheno_length = length(pheno)
    )
  )

  class(result) <- "gblup_cv_result"
  return(result)
}
