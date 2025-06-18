\name{select_hybrids}
\alias{select_hybrids}
\title{Selecting based on the predicted values of the hybrid offspring}
\usage{
Selecting based on the predicted values of the hybrid offspring
}
\description{
Selecting based on the predicted values of the hybrid offspring
}
\examples{
#假设已有预测结果
pre <- predict_hybrid_performance(hybrid_list, geno_new[,-1], pheno$HEI)
hybrid_gebv <- pre$predicted_values

# 选择前20%（按GEBV从高到低排序）
result <- select_hybrids(hybrid_gebv, sort_decreasing = TRUE,
                       method = "top_percent", value = 0.2)
 print(result$selected_values)

# 选择前10个（按GEBV从低到高排序）
result <- select_hybrids(hybrid_gebv, sort_decreasing = FALSE,
                      method = "top_n", value = 10)

#选择大于5的（按GEBV从高到低排序）
result <- select_hybrids(hybrid_gebv, sort_decreasing = TRUE,

}
