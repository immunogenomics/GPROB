
# GPROB <img src="man/figures/gprob.svg" width="181px" align="right" />

Multiple diseases can present with similar initial symptoms, making it
difficult to clinically differentiate between these conditions. GPROB
uses patients’ genetic information to help prioritize a diagnosis. This
genetic diagnostic tool can be applied to any situation with
phenotypically similar diseases with different underlying genetics.

<!-- badges: start -->

[![R build
status](https://github.com/immunogenomics/GPROB/workflows/R-CMD-check/badge.svg)](https://github.com/immunogenomics/GPROB/actions)
<!-- badges: end -->

## Citation

Please cite:

-   Knevel, R. et al. [Using genetics to prioritize diagnoses for
    rheumatology outpatients with inflammatory
    arthritis.](http://dx.doi.org/10.1126/scitranslmed.aay1548) Sci.
    Transl. Med. 12, (2020)

## License

Please see the [LICENSE](LICENSE) file for details. [Contact
us](mailto:soumya@broadinstitute.org) for other licensing options.

## Installation

Install and load the GPROB R package.

``` r
devtools::install_github("immunogenomics/GPROB")
library(GPROB)
```

## Synopsis

GPROB estimates the probability that each individual has a given
phenotype.

We need three inputs:

-   Population prevalences of the phenotypes of interest.

-   Odds ratios for SNP associations with the phenotypes.

-   SNP genotypes (0, 1, 2) for each individual.

### Example

Let’s use a small example with artificial data to learn how to use
GPROB.

Suppose we have 10 patients, and we know of 7 single nucleotide
polymorphisms (SNPs) associated with rheumatoid arthritis (RA) or
systemic lupus erythematosus (SLE).

#### Prevalence

First, we should find out the prevalence of RA and SLE in the population
that is representative of our patients.

``` r
prevalence <- c("RA" = 0.001, "SLE" = 0.001)
```

#### Odds Ratios

Next, we need to obtain the odds ratios (ORs) from published genome-wide
association studies (GWAS). We should be careful to note which alleles
are associated with the phenotype to compute the risk in the correct
direction.

``` r
or <- read.delim(
  sep = "",
  row.names = 1,
  text = "
snp  RA SLE
SNP1 1.0 0.4
SNP2 1.0 0.9
SNP3 1.0 1.3
SNP4 0.4 1.6
SNP5 0.9 1.0
SNP6 1.3 1.0
SNP7 1.6 1.0
")
or <- as.matrix(or)
```

#### Genotypes

Finally, we need the genotype data for each of our 10 patients. Here,
the data is coded in the form (0, 1, 2) to indicate the number of copies
of the risk allele.

``` r
geno <- read.delim(
  sep = "",
  row.names = 1,
  text = "
id SNP1 SNP2 SNP3 SNP4 SNP5 SNP6
 1    0    1    0    2    1    0
 2    0    0    1    0    2    2
 3    1    0    1    1    0    2
 4    1    1    0    2    0    0
 5    0    1    1    1    1    0
 6    0    0    1    3    0    2
 7    2    2    2    2    2    2
 8    1    2    0    2    1    1
 9    0    2    1   NA    1    2
10    1    0    2    2    2    0
")
geno <- as.matrix(geno)
```

#### Dealing with missing or invalid data

Before we run the `GPROB()` function, we need to deal with invalid and
missing data.

We remove individuals who have `NA` for any SNP:

``` r
ix <- apply(geno, 1, function(x) !any(is.na(x)))
geno <- geno[ix,]
```

We remove individuals who have invalid allele counts:

``` r
ix <- apply(geno, 1, function(x) !any(x < 0 | x > 2))
geno <- geno[ix,]
```

And we make sure that we use the same SNPs in the `or` and `geno`
matrices:

``` r
or <- or[colnames(geno),]
```

#### Run GPROB

Then we can run the GPROB function to estimate probabilities:

``` r
library(GPROB)
res <- GPROB(prevalence, or, geno)
res
#> $pop_prob
#>              RA          SLE
#> 1  0.0003556116 0.0017797758
#> 2  0.0033703376 0.0010049932
#> 3  0.0016672084 0.0006434285
#> 4  0.0003951084 0.0007126714
#> 5  0.0008885550 0.0014465506
#> 7  0.0005407850 0.0004337103
#> 8  0.0004622457 0.0006414499
#> 10 0.0003200618 0.0013374018
#> 
#> $cond_prob
#>           RA       SLE
#> 1  0.1665326 0.8334674
#> 2  0.7703046 0.2296954
#> 3  0.7215363 0.2784637
#> 4  0.3566669 0.6433331
#> 5  0.3805203 0.6194797
#> 7  0.5549386 0.4450614
#> 8  0.4188163 0.5811837
#> 10 0.1931034 0.8068966
```

In this example, we might interpret the numbers as follows:

-   Individual 2 has RA with probability 0.003, given individual genetic
    risk factors, disease prevalence, and the number of patients used in
    genetic risk score calculations.

-   Individual 2 has RA with probability 0.77, conditional on the
    additional assumption that individual 2 has either RA or SLE.

## Calculations, step by step

Let’s go through each step of GPROB to understand how how it works.

The genetic risk score <i>S<sub>ki</sub></i> of individual *i* for
disease *k* is defined as:

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\Large&space;S_{ki}=\sum_{j}{\beta_{kj}x_{ij}}"/>
</p>

where:

-   <i>x<sub>ij</sub></i> is the number of risk alleles of SNP *j* in
    individual *i*

-   <i>β<sub>kj</sub></i> is the log odds ratio for SNP *j* reported in
    a genome-wide association study (GWAS) for disease *k*

<table>
<tr>
<td>
<b>Note:</b> We might want to consider shrinking the risk by some factor
(e.g. 0.5) to correct for possible overestimation of the effect sizes
due to publication bias. In other words, consider running <code>geno \<-
0.5 \* geno</code>.
</td>
</tr>
</table>

``` r
risk <- geno %*% log(or)
risk
#>            RA         SLE
#> 1  -1.9379420  0.83464674
#> 2   0.3140075  0.26236426
#> 3  -0.3915622 -0.18392284
#> 4  -1.8325815 -0.08164399
#> 5  -1.0216512  0.62700738
#> 7  -1.5185740 -0.57856671
#> 8  -1.6755777 -0.18700450
#> 10 -2.0433025  0.54844506
```

The known prevalence <i>V<sub>k</sub></i> of each disease in the general
population:

``` r
prevalence
#>    RA   SLE 
#> 0.001 0.001
```

We can calculate the population level probability <i>P<sub>ki</sub></i>
that each individual has the disease.

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\Large&space;P_{ki}=\frac{1}{1+\exp{(S_{ki}-\alpha_k)}}"/>
</p>

We find <i>α<sub>k</sub></i> for each disease *k* by minimizing
<i>(P̅<sub>k</sub> - V<sub>k</sub>)<sup>2</sup></i>. This ensures that
the mean probability <i>P̅<sub>k</sub></i> across individuals is equal to
the known prevalence <i>V<sub>k</sub></i> of the disease in the
population.

``` r
# @param alpha A constant that we choose manually.
# @param risk A vector of risk scores for individuals.
# @returns A vector of probabilities for each individual.
prob <- function(alpha, risk) {
  1 / (
    1 + exp(alpha - risk)
  )
}
alpha <- sapply(seq(ncol(risk)), function(i) {
  o <- optimize(
    f        = function(alpha, risk, prevalence) {
      ( mean(prob(alpha, risk)) - prevalence ) ^ 2
    },
    interval = c(-100, 100),
    risk = risk[,i],
    prevalence = prevalence[i]
  )
  o$minimum
})
alpha
#> [1] 6.003374 7.164133
```

Now that we have computed alpha, we can compute the population-level
probabilities of disease for each individual.

``` r
# population-level disease probability
p <- sapply(seq_along(alpha), function(i) prob(alpha[i], risk[,i]))
p
#>            [,1]         [,2]
#> 1  0.0003556116 0.0017797758
#> 2  0.0033703376 0.0010049932
#> 3  0.0016672084 0.0006434285
#> 4  0.0003951084 0.0007126714
#> 5  0.0008885550 0.0014465506
#> 7  0.0005407850 0.0004337103
#> 8  0.0004622457 0.0006414499
#> 10 0.0003200618 0.0013374018
```

Next we assume that each individual has one of the diseases:

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\Large&space;\text{Pr}(Y_k=1|(\textstyle\sum_k{Y_k})=1)"/>
</p>

Then, we calculate the conditional probability <i>C<sub>ki</sub></i> of
each disease *k*:

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\Large&space;C_{ki}=\frac{P_{ki}}{\sum_k{P_{ki}}}"/>
</p>

``` r
# patient-level disease probability
cp <- p / rowSums(p)
cp
#>         [,1]      [,2]
#> 1  0.1665326 0.8334674
#> 2  0.7703046 0.2296954
#> 3  0.7215363 0.2784637
#> 4  0.3566669 0.6433331
#> 5  0.3805203 0.6194797
#> 7  0.5549386 0.4450614
#> 8  0.4188163 0.5811837
#> 10 0.1931034 0.8068966
```
