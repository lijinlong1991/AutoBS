\name{predict_hybrid }
\alias{predict_hybrid }
\title{check genotype and phenotype data}
\usage{
This function predicts the performance of hybrid offspring based on the genotype and phenotype data of the parents, using the GBLUP model.
}
\description{
This function predicts the performance of hybrid offspring based on the genotype and phenotype data of the parents, using the GBLUP model.}
\examples{
set.seed(123)
geno <- matrix(sample(0:2, 100*20, replace=TRUE), nrow=100, ncol=20)
rownames(geno) <- paste0("Genotype_", 1:100)
pheno <- rnorm(100)
hybrids <- data.frame(P1 = c("Genotype_1", "Genotype_2"),
                    P2 = c("Genotype_3", "Genotype_4"))

#' # 使用函数
result <- predict_hybrid_performance(hybrids, geno, pheno)
print(result$predicted_values)
}
