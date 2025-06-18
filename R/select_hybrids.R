select_hybrids <- function(predicted_values,
                           sort_decreasing = TRUE,
                           method = c("top_percent", "top_n", "threshold"),
                           value,
                           direction = c("greater", "less")) {

  # 参数验证
  method <- match.arg(method)
  direction <- match.arg(direction)

  if (missing(value)) {
    stop("必须提供value参数")
  }

  if (!is.numeric(predicted_values) || is.null(names(predicted_values))) {
    stop("predicted_values必须是有名称的数值向量")
  }

  # 创建包含名称和值的数据框
  data_df <- data.frame(
    Name = names(predicted_values),
    GEBV = as.numeric(predicted_values),
    stringsAsFactors = FALSE
  )

  # 按指定方向排序
  data_df <- data_df[order(data_df$GEBV, decreasing = sort_decreasing), ]

  # 根据选择方法筛选
  selected_idx <- switch(
    method,
    "top_percent" = {
      if (value <= 0 || value > 1) {
        stop("top_percent方法的value必须在0到1之间")
      }
      n_select <- ceiling(nrow(data_df) * value)
      1:n_select
    },
    "top_n" = {
      if (value <= 0 || value > nrow(data_df)) {
        stop("top_n方法的value必须在1到总杂交组合数之间")
      }
      1:value
    },
    "threshold" = {
      if (direction == "greater") {
        which(data_df$GEBV > value)
      } else {
        which(data_df$GEBV < value)
      }
    }
  )

  # 提取选中的行
  selected_data <- data_df[selected_idx, ]

  # 返回结果
  return(selected_data
  )
}
