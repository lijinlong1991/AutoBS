\name{modeling}
\alias{modeling}
\title{modeling genotype and phenotype data}
\usage{
Used to create model for genotype and phenotype data.
}
\description{
Building a genomic prediction model based on GBLUP and conduct cross-validation.
}
\examples{
modeling(geno = my_geno_data, pheno = my_pheno_data$trait, k = 5,plot = TRUE, seed = 123, verbose = TRUE)
}
