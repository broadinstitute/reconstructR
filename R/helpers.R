# MIT License
#
# Copyright (c) 2023 Ivan Specht
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

### Helper functions used in the outbreak reconstruction algorithm.

## Helper functions for MCMC

# Max possible value of infection time for a host
get_max_t_I <- function(i, l){
  min(c(l$t_E[which(l$anc == i)], l$t_test[i]))
}

# Gamma density function, given mean and variance
gamma_dens <- function(x, mean, var){
  if(all(x > 0) & mean > 0 & var > 0){
    dgamma(x, shape = mean^2/var, rate = mean/var, log = T)
  }else{
    -Inf
  }
}


# Survival function of Gamma distribution given mean and variance, evaluated at x - y
gamma_surv <- function(x, y, mean, var){
  if(mean > 0 & var > 0 & all(x > y, na.rm = T)){
    # In order for this to integrate to 1, must divide by the integral of survival function from 0 to infinity
    # But this is simply equal to the mean! So really, just need to divide by the mean
    out <- pgamma(x - y, shape = mean^2/var, rate = mean/var, lower.tail = F) / mean
    out[is.na(out)] <- 1
    out
  }else{
    0
  }
}

# Normalization of gamma density assuming t_I must exceed time to reach 1/sqrt(p) virions
t_I_norm <- function(l){
  # time to reach 1/sqrt(p) virions
  mean <- l$tau_E
  var <- l$var_E
  lambda <- l$mu/l$p
  delta_t <- (1/lambda)*log(1/(l$A[2:l$N] * sqrt(l$p)))
  pgamma(delta_t, shape = mean^2/var, rate = mean/var, lower.tail = F, log.p = T)
}

# Jukes-cantor evolution function
# input is a vector of minor allele frequencies, followed by vectors of 1st alleles (1,2,3,4 = A,C,G,T), followed by vector of 2nd alleles.
# Followed by mu (mutation rate) and T (time)
## Note: x is always the fraction of allele 2
evolve <- function(x, a1, a2, mu, e, t){
  # First, evolve by JC96
  out <- matrix(1/4 - (1/4)*exp(-(4/3) * (mu * t)), nrow = length(x), ncol = 4) # for entries that start as 0
  out[cbind(1:length(x), a2)] <- 1/4 + exp(-(4/3) * (mu * t)) * (x - 1/4)
  out[cbind(1:length(x), a1)] <- 1/4 + exp(-(4/3) * (mu * t)) * ((1-x) - 1/4)

  # Then, add error
  out <- out*(1-e) + (1-out)*e/3

  out
}

#microbenchmark({evolve(pdemo, maj, min, 0.1, 1)})

# Special distribution that captures de novo mutations in expo growth phase
dspecial <- function(x, A, p){
  -(((1-x)^(-1+A) * (-1+x-A* x+(1-x)^(1/sqrt(p)) *(1+(-1+A+(1/sqrt(p))) *x)))/((1/sqrt(p)) *x^2))
}

# Probability of no mutations being observed out of k reads, given bottleneck size, mutation rate, delta t, and p, assuming a mutation happens in expo growth phase
# pnomuts <- function(k, A, p, mu, t){
#   -(((1 - mu*t)^k * (-(1/sqrt(p)) - k * digamma(k + A) + k * digamma(k + A + (1/sqrt(p)))))/(1/sqrt(p)))
# }


# Probability density function that's equal to multinomial PMF for nonzero values in the data, and binomial CDF when MAF is reported at 0
# This is designed to account for MAFs below the the limit of detection / when read counts aren't reliable
trunc_dmnom <- function(x, size, prob, filters, log = FALSE){
  # Threshold value (gets reported as 0 if below this value)
  thresh <- pmax(size*filters$af, filters$call)

  # Which observations have just one allele reported?
  below_lod <- Rfast::rowsums(x > 0) == 1

  # Which allele is that?
  which_detected <- Rfast::rowMaxs(x[below_lod, , drop = F])

  # What is the probability of observing an alternate allele in the below_lod category?
  p_alt <- 1 - (prob[below_lod, , drop = F])[cbind(1:length(which_detected), which_detected)]

  # Suppressing warnings in case of extremely low log probabilities which may be reported as -Inf
  suppressWarnings({
    out1 <- pbinom(thresh[below_lod], size[below_lod], p_alt, log.p = log)
  })

  # And, for all other observations...
  out2 <- extraDistr::dmnom(x[!below_lod, ], size[!below_lod], prob[!below_lod, ], log = log)

  out <- rep(0, length(size))
  out[below_lod] <- out1
  out[!below_lod] <- out2
  out
}


## Efficient means of computing the Hypergeometric2F1 function at 1-n, 2-n, 2, z, for n ≥ 1

## LOG
poch <- function(q,n){
  if(n==0){
    0
  }else if(n==1){
    log(abs(q))
  }else{
    sum(
      log(abs(q + (0:(n-1))))
    )
  }
}

## LOG coefficients
get_h2f1_coefs <- function(a,b,c){
  if(a >= 0){
    return(NA)
  }else{
    m <- -a
    binoms <- lchoose(m, 0:m)
    bpoch <- sapply(0:m, poch, q = b)
    cpoch <- sapply(0:m, poch, q = c)

    return(
      binoms + bpoch - cpoch
    )
  }
}

h2f1 <- function(A, z, h2f1_coefs){
  if(A <= 0){
    return(NA)
  }else if(A==1){
    rep(1, length(z))
  }else if(A <= maxs$A){
    if(length(z)==1){
      sum(h2f1_coefs[[A]]*z^(0:(A-1)))
    }else{
      #zs <- sapply(z, function(z){z^(0:(A-1))})
      #colsums(h2f1_coefs[[A]]*zs)

      log_zs <- sapply(z, function(z){(0:(A-1)) * log(z)})
      Rfast::colsums(exp(log(h2f1_coefs[[A]]) + log_zs))
    }
  }else{
    NA
  }
}

## PDF for binomial-beta compound dist
# Note, not quite a PDF: would need to normalize by cases in which bottleneck is homozygous. But it's more useful in this form.
# Really, this is the sum of Beta(k, A-k) with k ~ Bin(A, q), summed over k from 1 to A-1.
dbinbeta <- function(x, A, q, h2f1_coefs){
  if(A <= 1 | any(q==1)){
    0
  }else if(length(x) == 0){
    1
  }else{

    z <- (q*x)/((-1+q) *(-1+x))
    log_fn_coefs <- log(A-1) + log(A) + log(1-q) + log(q) + (A-2)*(log(1-q) + log(1-x))

    if(length(z)==1){
      sum(exp(h2f1_coefs[[A]] + (0:(A-1))*log(z) + log_fn_coefs))
    }else{

      log_zs <- sapply(z, function(z){(0:(A-1)) * log(z)})
      Rfast::colsums(exp(h2f1_coefs[[A]] + log_zs + matrix(log_fn_coefs, ncol = length(log_fn_coefs), nrow = A, byrow = T)))
    }
  }
}

# Does the number of virions at time t_I exceed k = 1/sqrt(p)?
exceed <- function(l){
  lambda <- l$mu/l$p
  delta_t <- l$t_I[2:l$N] - l$t_E[2:l$N]
  n <- l$A[2:l$N] * exp(delta_t * lambda)
  if(any(n < 1/sqrt(l$p))){
    return(-Inf)
  }else{
    return(0)
  }
}

# Probability that by the end of the exponential growth phase, the number of mutated particles is less than filters$af
p_below_af <- function(p, sum_beta_CDF, n_outgoing){
  k <- floor(1/sqrt(p)) # Number of particles at end of expo growth phase

  # n_outgoing is number of outgoing particles.
  # if it exceeds the number of outgoing particles for which we've pre-computed the probabilities,
  # then we set it to the max number of outgoing particles we've accounted for.
  # this should be exceedingly rare (under defaults, means > 300 outgoing particles)

  if(n_outgoing >= length(sum_beta_CDF)){
    n_outgoing <- length(sum_beta_CDF) - 1
  }

  if(k <= length(sum_beta_CDF[[n_outgoing + 1]])){
    return(sum_beta_CDF[[n_outgoing + 1]][k]/k)
  }else{
    return((sum_beta_CDF[[n_outgoing + 1]][length(sum_beta_CDF[[n_outgoing + 1]])] + k - length(sum_beta_CDF[[n_outgoing + 1]]))/k)
  }
}



## Helper functions for data processing

# File for matching case names to their vcf files
match_vcf <- function(names, files){
  out <- c()
  for (i in 1:length(names)) {
    out[i] <- which(grepl(names[i], files))
  }
  out
}


# FASTA files use raw format. Here's a simple guide to how that works:
# a = 88 = as.raw(136)
# c = 28 = as.raw(40)
# g = 48 = as.raw(72)
# t = 18 = as.raw(24)

# Functions to convert between raw format and nucleotide letters:
base_to_raw <- function(b){
  if(b == "A"){
    as.raw(136)
  }else if(b == "C"){
    as.raw(40)
  }else if(b == "G"){
    as.raw(72)
  }else if(b == "T"){
    as.raw(24)
  }
}

raw_to_base <- function(r){
  if(r == as.raw(136)){
    "A"
  }else if(r == as.raw(40)){
    "C"
  }else if(r == as.raw(72)){
    "G"
  }else if(r == as.raw(24)){
    "T"
  }else{
    "N"
  }
}

# Get genomic coverage
get_coverage <- function(x){
  mean(x %in% as.raw(c(24,40,72,136)))
}

# Total alt allele frequency at a given position
get_tot_af <- function(pos, vcf){
  sum(vcf$ALT_COUNT[which(vcf$POS == pos)]) / (sum(vcf$ALT_COUNT[which(vcf$POS == pos)]) + vcf$REF_COUNT[which(vcf$POS == pos)][1])
}

# Read depth at a given position from VCF file
get_tot_count <- function(pos, vcf){
  (sum(vcf$ALT_COUNT[which(vcf$POS == pos)]) + vcf$REF_COUNT[which(vcf$POS == pos)][1])
}

# Helper function to perform fisher's exact test and extract p value
fisher <- function(s){
  s <- as.numeric(strsplit(s, ",")[[1]])
  fisher.test(matrix(s, ncol = 2))$p
}

# Process vcf file
process_vcf <- function(i, filepath, files, filters, cons){

  # Read in the file
  vcf <- try(read.table(paste0(filepath, "vcf/", files[i])), silent = T)

  if(class(vcf) == "try-error"){
    vcf <- matrix(nrow = 0, ncol = 13)
    colnames(vcf) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "REF_COUNT", "ALT_COUNT", "SB", "TOT_COUNT", "TOT_AF")
    vcf <- as.data.frame(vcf)
  }else{

    # Rename the columns (for clarity)
    colnames(vcf) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

    # Convert the INFO column to ref and alt reads, filtering for low values
    info <- vcf$INFO
    info <- gsub(".*;DP4=", "", info)
    counts <- sapply(info, strsplit, split = ",")
    counts <- (matrix(as.numeric(unlist(counts)), ncol = 4, byrow = T))
    ref <- counts[,1] + counts[,2]
    ref[ref < filters$call] <- 0
    alt <- counts[,3] + counts[,4]
    alt[alt < filters$call] <- 0

    # DP4 category
    dp4 <- sub(".*DP4=", "", vcf$INFO)
    dp4_pval <- sapply(dp4, fisher, USE.NAMES = F)

    # Strand bias
    sb <- gsub(".*;SB=", "", vcf$INFO)
    sb <- sub(";.*","",sb)
    sb <- as.numeric(sb)

    # Read depth
    dp <- ref + alt

    # Sites that do not pass our filters
    problem <- which(sb > filters$sb | dp4_pval < filters$dp4 | dp < filters$depth | (ref / (ref + alt)) < filters$af | (alt / (ref + alt)) < filters$af)
    cons_problem <- sapply(cons[[i]][vcf$POS[problem]], raw_to_base)

    ref[problem][cons_problem == vcf$REF[problem]] <- dp[problem][cons_problem == vcf$REF[problem]]
    alt[problem][cons_problem == vcf$REF[problem]] <- 0

    ref[problem][cons_problem == vcf$ALT[problem]] <- 0
    alt[problem][cons_problem == vcf$ALT[problem]] <- dp[problem][cons_problem == vcf$ALT[problem]]

    # If equals neither, set to 0
    ref[problem][cons_problem != vcf$REF[problem] & cons_problem != vcf$ALT[problem]] <- 0
    alt[problem][cons_problem != vcf$REF[problem] & cons_problem != vcf$ALT[problem]] <- 0

    #ref[ref < alt & (sb > filters$sb | dp4_pval < filters$dp4)] <- 0
    #alt[alt <= ref & (sb > filters$sb | dp4_pval < filters$dp4)] <- 0

    # Filter by MAF
    # ref[ref < alt & (ref / (ref + alt)) < filters$af] <- 0
    # alt[alt < ref & (alt / (ref + alt)) < filters$af] <- 0

    # Append new columns
    vcf$REF_COUNT <- ref
    vcf$ALT_COUNT <- alt
    vcf$SB <- sb

    # Check that we have consensus data at all relevant positions
    # vcf <- vcf[cons[[i]][vcf$POS] %in% as.raw(c(136,40,72,24)), ]

    # Delete observations with no iSNV (after filtering)
    vcf <- vcf[vcf$ALT_COUNT > 0, ]
    #vcf <- vcf[vcf$POS > 265 & vcf$POS < 29558, ]

    # Make a new column with total counts
    TOT_COUNT <- sapply(vcf$POS, get_tot_count, vcf = vcf)
    vcf <- cbind(vcf, TOT_COUNT)

    # Delete rows with too low of a read depth
    # vcf <- vcf[vcf$TOT_COUNT >= filters$depth, ]

    # Make a new column with total frequency of an alternate nucleotide
    TOT_AF <- sapply(vcf$POS, get_tot_af, vcf = vcf)
    vcf <- cbind(vcf, TOT_AF)

    # Info column is now alternate frequency
    vcf$INFO <- vcf$ALT_COUNT / vcf$TOT_COUNT

    # Delete observations where we're not actually making a substitution
    vcf <- vcf[vcf$REF %in% c("A", "C", "G", "T"), ]
    vcf <- vcf[vcf$ALT %in% c("A", "C", "G", "T"), ]
  }

  vcf
}

# Convert VCF into matrix format (see process.R)
vcf_to_read <- function(i, vcfs, all_pos, depth, cons, filters){

  read <- vcfs[[i]]

  # New columns for counts of each nucleotide
  read$A <- rep(0, nrow(read))
  read$C <- rep(0, nrow(read))
  read$G <- rep(0, nrow(read))
  read$T <- rep(0, nrow(read))

  # Loop through all iSNVs and add them to "read"
  if(nrow(read) > 0){
    for (j in 1:nrow(read)) {
      base1 <- read$REF[j]
      base2 <- read$ALT[j]
      n_base1 <- read$REF_COUNT[j]
      n_base2 <- read$ALT_COUNT[j]

      row <- min(which(read$POS == read$POS[j]))
      read[row, base1] <- read[row, base1] + n_base1
      read[row, base2] <- read[row, base2] + n_base2
    }
  }

  # Delete additional rows at the same position
  read <- read[!duplicated(read$POS), ]

  # Get relevant columns
  read <- read[, c("POS", "A", "C", "G", "T")]

  # Filter to positions in all_pos
  read <- read[read$POS %in% all_pos, ]

  # If iSNV already exists at a certain position, do nothing. Otherwise, add the read for that position using consensus and depth info.
  for (p in all_pos) {
    if(!(p %in% read[, "POS"])){
      newrow <- rep(0, 4)
      newrow[which(c("A","C","G","T") == raw_to_base(cons[[i]][p]))] <- depth[p, i-1]
      if(max(newrow) < filters$depth){
        newrow <- rep(0, 4)
      }
      newrow <- c(p, newrow)
      read <- rbind(read, newrow)
      colnames(read) <- c("POS", "A", "C", "G", "T")
    }
  }

  # Order by position
  read <- read[order(read[,"POS"]), ]

  if(nrow(read) > 0){
    rownames(read) <- 1:nrow(read)
  }


  # Now that we've reordered, we don't need the POS column anymore
  read <- read[, c("A", "C", "G", "T")]

  colnames(read) <- NULL
  rownames(read) <- NULL

  as.matrix(read)
}

# Person-to-person SNP distance
get_cons_dist <- function(i,j,reads){
  if(i==1){
    NA
  }else if(j==1){
    NA
  }else{
    sum(Rfast::rowMaxs(reads[[i]]) != Rfast::rowMaxs(reads[[j]]) & Rfast::rowsums(reads[[i]]) > 0 & Rfast::rowsums(reads[[j]]) > 0)
  }
}

# Do two people share at least one allele at every site?
get_overlap <- function(i,j,reads){
  if(i==1){
    NA
  }else if(j==1){
    NA
  }else{
    all(Rfast::rowsums(reads[[i]] * reads[[j]]) > 0 | Rfast::rowsums(reads[[i]]) == 0 | Rfast::rowsums(reads[[j]]) == 0)
  }
}

# Are there any sites at which two people share MORE than one allele?
get_shared_isnv <- function(i,j,reads){
  if(i==1){
    NA
  }else if(j==1){
    NA
  }else{
    any(Rfast::rowsums(reads[[i]] * reads[[j]] > 0) >= 2) & get_overlap(i,j,reads)
  }
}

# Summarize relationship
summarize <- function(i,j,reads){
  if(i==1 | j==1){
    "Index Case"
  }else if(get_shared_isnv(i,j,reads)){
    "Split Bottleneck"
  }else if(
    get_overlap(i,j,reads) &
    any(
      Rfast::rowsums(reads[[i]] > 0) >= 2 & # iSNV in donor
      Rfast::rowMaxs(reads[[i]]) != Rfast::rowMaxs(reads[[j]]) & # consensus change
      Rfast::rowsums(reads[[i]]) > 0 & # Data for donor
      Rfast::rowsums(reads[[j]]) > 0 # Data for recipient
    )
  ){
    "Minor Allele to Fixation"
  }else if(
    get_overlap(i,j,reads) &
    any(
      Rfast::rowsums(reads[[j]] > 0) >= 2 & # iSNV in recipient
      Rfast::rowMaxs(reads[[i]]) != Rfast::rowMaxs(reads[[j]]) & # consensus change
      Rfast::rowsums(reads[[i]]) > 0 & # Data for donor
      Rfast::rowsums(reads[[j]]) > 0 # Data for recipient
    )
  ){
    "Fixation to Minor Allele"
  }else if(get_cons_dist(i,j,reads) == 0){
    "Same Consensus"
  }else{
    "Undetected Consensus Change"
  }
}

## Miscellaneous helper functions

# Get the largest value per row of a matrix, or NA if every element of the row is 0.
get_max_or_na <- function(r){
  out <- Rfast::rowMaxs(r)
  # Assign NA if no read data
  out[Rfast::rowsums(r) == 0] <- NA
  out
}

# Return 1 for max in vector, 0 otherwise
get_max_pos <- function(v){
  out <- rep(0, length(v))
  out[which.max(v)] <- 1
  out
}

#### Contact tracing likelihood
get_n_true_pos <- function(anc, contacts){
  sum(mapply(function(x,y){x %in% y}, anc, contacts))
}

get_n_false_neg <- function(anc, contacts){
  sum(mapply(function(x,y){!(x %in% y)}, anc, contacts)) - sum(anc[!is.na(anc)] == 1) - 1 # don't penalize index cases, and subtract 1 because ancestor of "host 1" is NA
}

get_ct_lik <- function(anc, contacts, nu, xi, N_ct){
  n_true_pos <- get_n_true_pos(anc, contacts)
  n_false_pos <- N_ct - n_true_pos

  n_false_neg <- get_n_false_neg(anc, contacts)
  n_true_neg <- (length(anc)*(length(anc) - 1) / 2) - n_false_neg # of ordered pairs - number of pairs who do transmit, but don't have CT

  # Likelihood is:
  # nu^n_true_pos * xi^n_false_positive *
  # (1-nu)^n_false_neg * (1-xi)^n_true_neg
  return(
    n_true_pos * log(nu) + n_false_pos * log(xi) +
      n_false_neg * log(1-nu) + n_true_neg * log(1-xi)
  )

}








