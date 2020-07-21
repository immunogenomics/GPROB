# GPROB calculates genetic probability of multiple diseases
# Copyright (C) 2020  Rachel Knevel
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' GPROB
#'
#' Many diseases may have similar presentations in the clinic. By
#' using genetic information, it may sometimes be possible to prioritize some
#' diagnoses over others. GPROB implements calculations of genetic probability
#' from genotypes and effect sizes for genetic variants associated with
#' diseases.
#'
#' @param prevalence A numeric vector with prevalence of each phenotype in the population.
#' @param or A matrix with odds ratios for each phenotype (column) and genetic variant (row).
#' @param geno A matrix with risk allele counts for each person (row) and genetic variant (column).
#' @return A list with the population-level probabilities for each person to
#' have each phenotype, and the conditional probabilities, assuming that each
#' person has one of the phenotypes.
#' @export
#' @importFrom Rdpack reprompt
#' @importFrom stats optimize
#' @references{
#' \insertRef{Knevel2020}{GPROB}
#' }
#' @examples
#' prevalence <- c("RA" = 0.001, "SLE" = 0.001)
#' or <- read.delim(
#'   sep = "",
#'   row.names = 1,
#'   text = "
#' snp  RA SLE
#' SNP1 1.0 0.4
#' SNP2 1.0 0.9
#' SNP3 1.0 1.3
#' SNP4 0.4 1.6
#' SNP5 0.9 1.0
#' SNP6 1.3 1.0
#' SNP7 1.6 1.0
#' ")
#' or <- as.matrix(or)
#' geno <- read.delim(
#'   sep = "",
#'   row.names = 1,
#'   text = "
#' id SNP1 SNP2 SNP3 SNP4 SNP5 SNP6
#'  1    0    1    0    2    1    0
#'  2    0    0    1    0    2    2
#'  3    1    0    1    1    0    2
#'  4    1    1    0    2    0    0
#'  5    0    1    1    1    1    0
#'  6    0    0    1    3    0    2
#'  7    2    2    2    2    2    2
#'  8    1    2    0    2    1    1
#'  9    0    2    1   NA    1    2
#' 10    1    0    2    2    2    0
#' ")
#' geno <- as.matrix(geno)
#' ix <- apply(geno, 1, function(x) !any(is.na(x)))
#' geno <- geno[ix,]
#' ix <- apply(geno, 1, function(x) !any(x < 0 | x > 2))
#' geno <- geno[ix,]
#' or <- or[colnames(geno),]
#' GPROB(prevalence, or, geno)
GPROB <- function(prevalence, or, geno) {
  # Risk scores for each individual.
  risk <- geno %*% log(or)
  # Population-level probability for a risk score.
  prob <- function(alpha, risk) {
    1 / ( 1 + exp(alpha - risk) )
  }
  # We find alpha by optimization to ensure that the mean
  # population-level probability is equal to the prevalence.
  alpha <- sapply(seq(ncol(risk)), function(i) {
    o <- optimize(
      f = function(alpha, risk, prevalence) {
        ( mean(prob(alpha, risk)) - prevalence ) ^ 2
      },
      interval = c(-100, 100),
      risk = risk[,i],
      prevalence = prevalence[i]
    )
    o$minimum
  })
  # Population-level disease probability.
  p <- sapply(seq_along(alpha), function(i) prob(alpha[i], risk[,i]))
  colnames(p) <- names(prevalence)
  # Conditional disease probability, assuming each person has one of the
  # phenotypes.
  cp <- p / rowSums(p)
  list(pop_prob = p, cond_prob = cp)
}
