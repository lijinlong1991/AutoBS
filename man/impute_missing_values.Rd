\name{impute_missing_values}
\alias{impute_missing_values}
\title{impute genotype and phenotype data}
\usage{
Used to impute genotype and phenotype data.
}
\description{
It is used to fill in the missing values in genotype and phenotype data, mainly by using the average value method.
}
\examples{
impute_missing_values(geno, pheno, geno_id_col = 1, pheno_id_col = 1)
}
