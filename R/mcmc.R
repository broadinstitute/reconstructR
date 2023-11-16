#' Run MCMC Algorithm
#'
#' This function executes the MCMC algorithm for outbreak reconstruction.
#'
#' @param N_iters Number of MCMC iterations.
#' @param filters Filters for input data. NULL = defaults.
#' @param mins Minimum values of model parameters. NULL = defaults.
#' @param maxs Maximum values of model parameters. NULL = defaults.
#' @param vars Variances of MCMC moves for model paramaters. NULL = defaults.
#' @param seed Whether or not to set the seed, for replicable results. Defaults to TRUE.
#' @return A list, where each entry is the current state of the Markov chain (multiples of 100 steps).
#' @export

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

run_mcmc <- function(
    N_iters = 1e5,
    filters = NULL,
    mins = NULL,
    maxs = NULL,
    vars = NULL,
    prior_params = NULL,
    seed = T
){

  ### Establish default values for filters, mins, maxs, and vars

  if(is.null(filters)){
    filters <- list(
      omit = character(0), # any samples we should omit?
      coverage = 0.75, # minimum coverage
      sb = Inf, # filter for strand bias (Inf if not using)
      dp4 = 0.05, # p-value cutoff for fishers exact test on ref-fwd, ref-rev, alt-fwd, alt-rev
      af = 0.03, # filter for minor allele frequency
      call = 10, # minimum reads of the alternate allele to call iSNV
      depth = 100, # minimum total depth
      problem_qtile = 0.05 # if an iSNV is one of the problem_qtile-most prevalent in the Broad dataset, throw it out
    )
  }

  if(is.null(mins)){

    mins <- list()
    mins$mu <- 0 # mutation rate per day
    mins$e <- 0 # sequencing error rate
    mins$gamma <- 0 # mean bottleneck size
    mins$alpha <- 0 # probability of index case
    mins$p <- 0 # probability of mutation per site per replication event
    mins$w <- 2 # number of offspring of a virion per cycle
    mins$tau_T <- 0 # time from becoming infectious to time of test. Assume FOR NOW no positive tests before infectious.
    mins$tau_E <- 0 # time from exposed to infectious
    mins$tau_I <- 0 # time in infectious state
    mins$var_T <- 0.1
    mins$var_E <- 0.1
    mins$var_I <- 0.1

    mins$A <- 1 # minimum bottleneck size
  }

  if(is.null(maxs)){
    maxs <- list()
    maxs$mu <- 0.01 # mutation rate per day
    maxs$e <- 0.01 # sequencing error rate
    maxs$gamma <- 10
    maxs$alpha <- 1 # probability of index case
    maxs$p <- 0.01
    maxs$w <- 10000 # number of offspring of a virion per cycle
    maxs$tau_T <- 10 # time from becoming infectious to time of test. Assume FOR NOW no positive tests before infectious.
    maxs$tau_E <- 10 # time from exposed to infectious
    maxs$tau_I <- 10 # time in infectious state
    maxs$var_T <- 15
    maxs$var_E <- 15
    maxs$var_I <- 15

    maxs$A <- 300 # max bottleneck size
  }

  if(is.null(vars)){
    # For each of the epidemiological parameters, specify the standard deviation of the (normally-distributed) proposal density
    vars <- list()
    vars$mu <- 1e-6
    vars$e <- 1e-12
    vars$gamma <- 0.5
    vars$alpha <- 1e-9
    vars$p <- 1e-6
    vars$w <- 50 # number of offspring of a virion per cycle

    # We also have a proposal density for all of the x's, as well as CT
    vars$x <- 1e-2

    vars$tau_T <- 0.5 # time from becoming infectious to time of test
    vars$tau_E <- 0.5 # time from exposed to infectious
    vars$tau_I <- 0.5
    vars$var_T <- 0.5
    vars$var_E <- 0.5
    vars$var_I <- 0.5
    vars$t_I <- 1
    vars$t_E <- 1
  }

  if(is.null(prior_params)){
    # For each of the epidemiological parameters, specify the standard deviation of the (normally-distributed) proposal density
    prior_params <- list()
    prior_params$alpha <- c(1, 1e10)
  }



  ### First step: process deep-sequencing data in the form of VCF files, a FASTA file, and a read depth CSV into a useful format for outbreak reconstruction.

  message("Processing input data...")

  filepath <- "./input_data/"

  # Read in consensus genomes
  cons <- ape::read.FASTA(paste0(filepath, "aligned.fasta"))

  # Read the reference sequence
  master <- ape::read.FASTA(paste0(filepath, "ref.fasta"))
  names(master) <- "NC_045512.2"

  # Get the coverage for each consensus sequence
  coverage <- sapply(cons, get_coverage)

  # Determine which consensus sequences don't meet the coverage threshold
  low_coverage <- as.numeric(which(coverage < filters$coverage))

  # Omit low coverage sequences, as well as any other sequences the user specifies we should omit
  omit <- union(low_coverage, match(filters$omit, names(cons)))
  omit <- omit[!is.na(omit)]
  if(length(omit) > 0){
    cons <- cons[-omit]
  }

  # Test dates
  t_test <- read.csv("./input_data/date.csv")
  t_test <- t_test[match(names(cons), t_test[,1]), 2]

  # Make sure min of t_test is 0
  t_test <- t_test - min(t_test)

  # Add an entry to the front of t_test to offset for the test date of index case in larger population. Fixed to -Inf.
  t_test <- c(-Inf, t_test) # To account for index case

  ## Read and process the contact tracing file, if present
  ct_yes <- F
  if("ct.csv" %in% list.files(filepath)){

    ct <- read.csv(paste0(filepath, "ct.csv"))

    # Filter to cases in "cons"
    ct <- ct[ct[,1] %in% names(cons) & ct[,2] %in% names(cons), , drop = F]
    ct[,1] <- match(ct[,1], names(cons))
    ct[,2] <- match(ct[,2], names(cons))

    if(nrow(ct) > 0){
      ct_yes <- T

      # Add 1, to match new indexing of cases
      ct <- ct + 1

      # Store number of contact links identified
      N_ct <- nrow(ct)

      ## Reformatting: create a list of contacts for each individual. Will be helpful for computing likelihood.
      contacts <- list()
      contacts[[1]] <- numeric(0)
      for (j in 2:(length(cons) + 1)) {
        sub <- ct[ct[,1] == j | ct[,2] == j, ]
        sub <- unique(as.numeric(as.matrix(sub)))
        sub <- sort(setdiff(sub, j))
        contacts[[j]] <- sub
      }
    }
  }

  # Read and process the read depth file
  depth <- read.csv(paste0(filepath, "depth.csv"), check.names = F)

  # Make sure the column names match the names on the consensus sequences
  depth <- depth[, match(names(cons), colnames(depth))]

  # Number of cases
  N <- length(cons)

  # VCF file paths
  files <- list.files(path = paste0(filepath, "vcf"))

  # Filter to files that match the cases of interest. In order to find a match, the name of the case must be a substring of the file name
  files <- files[match_vcf(names(cons), files)]

  # Process the VCFs (adding a blank file at the front to represent the index case)
  vcfs <- lapply(1:N, process_vcf, filepath = filepath, files = files, filters = filters, cons = cons)
  blank_vcf <- vcfs[[1]][integer(),]
  vcfs <- c(list(blank_vcf), vcfs)

  # Merge wild type and consensus genomes into single VCF file
  cons <- insect::join(master, cons)

  # Adjust N to include a designator for index case, i.e. a transmission 1-->2 we will ultimately take to mean "2 is an index case"
  N <- length(cons)

  # Extract all positions on the viral genome across all hosts that exhibit iSNVs
  all_pos <- c()
  for (i in 1:N) {
    all_pos <- c(all_pos, vcfs[[i]]$POS)
  }
  all_pos <- sort(unique(all_pos))


  # Reformat vcf files into a length(all_pos) x 4 matrix, where entry [i,j] is the number of nucleotide j at site i
  # Take j=1 to mean "A", j=2 to mean "C", j=3 to mean "G", j=4 to mean "T"
  reads <- list()
  reads[[1]] <- matrix(0, nrow = length(all_pos), ncol = 4)
  reads[2:N] <- lapply(2:N, vcf_to_read, vcfs = vcfs, all_pos = all_pos, depth = depth, cons = cons, filters = filters)

  # Sum of all reads of all nucleotides across all people
  all_reads <- Reduce("+", reads)


  ## Obtain a list of problematic sites on the genome, based on those that consistently demonstrate iSNVs
  context <- read.csv(system.file("extdata", "context.csv", package = "reconstructR"))
  tot_isnvs <- data.frame(MUT = context$aa_change_full, POS = context$POS, COUNT = Rfast::rowsums(as.matrix(context[,4:ncol(context)]), na.rm = T))
  N_context <- 172519
  N_isnv_cases <- sum(tot_isnvs$COUNT)

  hgeom_probs <- c()
  for (i in 1:nrow(depth)) {
    sub <- tot_isnvs[tot_isnvs$POS == i, ]
    hgeom_probs[i] <- 1 - phyper(1, sum(sub$COUNT), N_context - sum(sub$COUNT), 1000) # Probability
  }


  # What sites have a probability exceeding 0.05?
  problem_sites <- which(hgeom_probs > filters$problem_qtile)

  ###

  # At which sites among all_pos are the consensus genomes all the same?
  if(length(all_pos) > 0){
    all_row_max <- sapply(reads[2:N], get_max_or_na)
    if(length(all_pos == 1)){
      all_row_max = matrix(all_row_max, nrow = 1)
    }
    cons_change <- apply(all_row_max, 1, function(v){length(unique(v[!is.na(v)]))}) > 1

    problem_sites <- setdiff(problem_sites, all_pos[cons_change])

    # Delete sites with common iSNVs AND no consensus-level change
    to_delete <- which(all_pos %in% problem_sites)
    to_delete <- sort(union(to_delete, which(Rfast::rowsums(all_reads > 0) < 2)))
    to_keep <- setdiff(1:length(all_pos), to_delete)
  }else{
    to_keep <- c()
  }



  # Revise reads and all_pos to only include sites with some mutation across the dataset
  for (i in 1:N) {
    reads[[i]] <- reads[[i]][to_keep, , drop = F]
  }
  all_pos <- all_pos[to_keep]
  L <- length(all_pos)

  all_reads <- all_reads[to_keep, , drop = F]

  # Create a matrix of the most-commonly-read nucleotide at each site, across all hosts.
  common <- matrix(0, nrow = length(all_pos), ncol = 4)
  if(length(all_pos) > 0){
    common[cbind(1:length(all_pos), Rfast::rowMaxs(all_reads))] <- 1
  }

  # Sum of (read depth > filters$depth) for sites not included in all_pos
  not_all_pos <- setdiff(1:nrow(depth), all_pos)
  K <- length(not_all_pos)

  depths <- list()
  sum_depths <- c()
  depths[[1]] <- rep(0, K)
  sum_depths[1] <- 0
  for (i in 2:N) {
    depths[[i]] <- depth[not_all_pos,i-1]
    sum_depths[i] <- sum(depths[[i]] > filters$depth)
  }

  # Extra step: store coefficients for hypergeometric function, which depends on maxs
  h2f1_coefs <- list()
  for (i in 1:maxs$A) {
    h2f1_coefs[[i]] <- get_h2f1_coefs(1-i, 2-i, 2)
  }

  # Also, store list of beta CDFs evaluated at filter$af
  # the ith entry is the CDF evaluated when we have i non-mutated particles and 1 mutated particle
  # N_beta_CDF is how many values we want; afterwards, assumed to be 1
  # N_tot_outgoing is the max number of outgoing particles for which we pre-compute
  N_beta_CDF <- 1e5
  N_tot_outgoing <- maxs$A

  sum_beta_CDF <- list()

  for (s in 0:N_tot_outgoing) {
    beta_CDFs <- pbeta(filters$af, 1, (1+s):(N_beta_CDF+s))

    # helpful in MCMC algo: cumsum of the beta CDFs
    sum_beta_CDF[[s+1]] <- cumsum(beta_CDFs * 1:N_beta_CDF / ((1+s):(N_beta_CDF+s)))

    #print(s)
  }



  ####### End of data processing.



  ######## Initialize the markov chain
  message("Initializing Markov chain...")

  mcmc <- list()
  mcmc$N <- N # number of people (including "case 1" which designates index case)
  mcmc$L <- L # number of sites that exhibit a mutation or iSNV in at least one case
  mcmc$K <- K # number of sites that never exhibit a mutation or iSNV in any case
  if(ct_yes){
    mcmc$N_ct <- N_ct # Number of CT links
    mcmc$nu <- 0.5 # true positive rate for CT
    mcmc$xi <- 0.5 # false positive rate for CT
    mcmc$contacts <- contacts

    # Also set mins, maxs, and proposals
    vars$nu <- 0.1
    vars$xi <- 0.1
    mins$nu <- 0
    mins$xi <- 0
    maxs$nu <- 1
    maxs$xi <- 1
  }else{
    mcmc$N_ct <- -1 # Placeholder for no contact tracing
  }
  mcmc$mu <- 1e-6 # mutation rate parameter after exponential growth phase
  mcmc$e <- 1e-10 # error rate parameter
  mcmc$gamma <- 2 # bottleneck size is Geom(1/(1+gamma))
  mcmc$alpha <- min(1e-10, maxs$alpha / 2) # number of index cases is 1 + Bin(N-1, alpha)
  mcmc$p <- 1e-6 # substitution rate at given nucleotide per RNA replication event
  mcmc$w <- 1000 # number of offspring of a virion per cycle

  mcmc$tau_T <- 2 # Mean duration of sojourn interval
  mcmc$tau_E <- 2 # Mean duration of exposed period
  mcmc$tau_I <- 6 # Mean duration of infectious period
  mcmc$var_T <- 2 # Variance of duration of sojourn interval
  mcmc$var_E <- 2 # Variance of duration of exposed period
  mcmc$var_I <- 6 # Variance of duration of infectious period

  mcmc$t_I <- t_test - mcmc$tau_T # Time of exposed-to-infectious transition for each host
  mcmc$t_E <- mcmc$t_I - 2 # Time of exposure for each host
  mcmc$t_test <- t_test # store time of test (never changes)

  mcmc$A <- rep(1, N) # Bottleneck sizes

  ### Allele fractions
  mcmc$x <- list() # vector of 2nd allele fractions at end of expo growth phase

  mcmc$a1 <- list() # vector of 1st allele (1 thru 4)

  mcmc$a2 <- list() # vector of 2nd allele (1 thru 4)

  # Should we start with all x's for missing data points as 0 or 1? (except for index case)
  start_x_binary <- T

  for (i in 1:N) {

    if(i==1){
      ## Most common allele across all reads
      majors1 <- Rfast::rowMaxs(all_reads)
      ## 2nd most common allele across all reads
      mod_all_reads <- all_reads
      if(length(all_pos) > 0){
        mod_all_reads[cbind(1:length(all_pos), majors1)] <- 0
      }
      minors1 <- Rfast::rowMaxs(mod_all_reads)

      # 1st allele (in order: A,C,G,T = 1,2,3,4)
      a1s <- pmin(majors1, minors1)
      # 2nd allele
      a2s <- pmax(majors1, minors1)

      mcmc$a1[[i]] <- a1s
      mcmc$a2[[i]] <- a2s

      mcmc$x[[i]] <- all_reads[cbind(1:L, a2s)] / Rfast::rowsums(all_reads)
    }else{

      ## Most common allele across all reads
      majors <- Rfast::rowMaxs(reads[[i]])
      ## 2nd most common allele across all reads
      mod_all_reads <- reads[[i]]
      if(length(all_pos) > 0){
        mod_all_reads[cbind(1:length(all_pos), majors)] <- 0
      }
      minors <- Rfast::rowMaxs(mod_all_reads)

      # If 0 reads, default to major/minor for overall outbreak
      majors[Rfast::rowsums(reads[[i]]) == 0] <- majors1[Rfast::rowsums(reads[[i]]) == 0]
      minors[Rfast::rowsums(mod_all_reads) == 0] <- minors1[Rfast::rowsums(mod_all_reads) == 0]

      # And finally, if we've assigned a minor allele that's the same as the major allele, switch it to the other option
      minors[minors == majors] <- majors1[minors == majors]


      # 1st allele (in order: A,C,G,T = 1,2,3,4)
      a1s <- pmin(majors, minors)
      # 2nd allele
      a2s <- pmax(majors, minors)

      mcmc$a1[[i]] <- a1s
      mcmc$a2[[i]] <- a2s
      mcmc$x[[i]] <- reads[[i]][cbind(1:L, mcmc$a2[[i]])] / Rfast::rowsums(reads[[i]])

      # If no read data at this site, default to proportion in overall outbreak
      missing <- which(is.nan(mcmc$x[[i]]))
      mcmc$x[[i]][missing] <- mcmc$x[[1]][missing]

      # If binary, simply round the start
      if(start_x_binary){
        mcmc$x[[i]][missing] <- round(mcmc$x[[1]][missing])
      }
    }
  }


  # Ancestor for each case (initialized to everyone being an index case)
  mcmc$anc <- rep(1,N)
  mcmc$anc[1] <- NA

  # Epidemiological likelihood (for the outbreak at large)
  mcmc$e_lik <- get_e_lik(mcmc, mins, maxs, prior_params)

  # Genomic likelihood (per person)
  if(mcmc$L > 0){
    mcmc$g_lik <- sapply(1:N, get_g_lik, l = mcmc, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs)
    if(mcmc$L == 1){
      mcmc$g_lik <- matrix(mcmc$g_lik, nrow = 1)
    }
  }else{
    mcmc$g_lik <- matrix(nrow = 0, ncol = N)
  }


  # Genomic likelihood for unmutated sites
  mcmc$u_lik <- sapply(1:N, get_u_lik, l = mcmc, sum_depths = sum_depths, filters = filters, sum_beta_CDF = sum_beta_CDF)


  ####### Now we can run the Markov chain



  if(seed){
    set.seed(110) # Homage to Professor Joe Blitzstein's Stat 110 course at Harvard, which first taught me about Markov chains
  }



  new_res <- mcmc

  results <- rep(list(new_res), N_iters / 100)

  for (count in 1:N_iters) {
    # Update epi params

    mcmc <- update(mcmc, "gamma", reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, vars, sum_beta_CDF)
    mcmc <- update(mcmc, "alpha", reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, vars, sum_beta_CDF)
    mcmc <- update(mcmc, "mu", reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, vars, sum_beta_CDF)
    mcmc <- update(mcmc, "e", reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, vars, sum_beta_CDF)
    mcmc <- update(mcmc, "p", reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, vars, sum_beta_CDF)
    mcmc <- update(mcmc, "w", reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, vars, sum_beta_CDF)

    if(ct_yes){
      mcmc <- update(mcmc, "nu", reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, vars, sum_beta_CDF)
      mcmc <- update(mcmc, "xi", reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, vars, sum_beta_CDF)
    }

    mcmc <- update(mcmc, "tau_T", reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, vars, sum_beta_CDF)
    mcmc <- update(mcmc, "tau_E", reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, vars, sum_beta_CDF)

    # Not updating variance in soujourn/exposed period
    #mcmc <- update(mcmc, "var_T")
    #mcmc <- update(mcmc, "var_E")

    mcmc <- update_t_E(mcmc, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF)
    mcmc <- update_t_E_ancestor(mcmc, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF)
    mcmc <- update_t_I(mcmc, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF)
    if(mcmc$L > 0){
      mcmc <- update_ancestor_x(mcmc, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF)
    }
    mcmc <- update_swap(mcmc, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF)
    mcmc <- update_ancestor(mcmc, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF)
    mcmc <- update_neck(mcmc, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF)
    mcmc <- update_swap_kids(mcmc, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF)
    if(mcmc$L > 0){
      mcmc <- update_x(mcmc, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, vars, sum_beta_CDF)
    }
    mcmc <- update_epi_ancestor(mcmc, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF)

    # Code to ensure we're correctly calculating the genomic likelihood
    # if(!identical(mcmc$g_lik, sapply(1:N, get_g_lik, l = mcmc, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs))){
    #   print("OH NO")
    #   print(count)
    #   print(sum(mcmc$g_lik))
    #   print(sum(sapply(1:N, get_g_lik, l = mcmc, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs)))
    # }
    #
    # if(!identical(mcmc$e_lik, get_e_lik(mcmc, mins=mins, maxs=maxs))){
    #   print("OH NO 2")
    #   print(count)
    # }

    if(count %% 100 == 0){

      message(paste(count, "iterations complete. Log-likeihood =", round(sum(mcmc$g_lik) + sum(mcmc$u_lik) + mcmc$e_lik, 2)))
      new_res <- mcmc
      results[[count / 100]] <- new_res
    }
  }

  # Store reads in results[[1]], for visualization purposes
  results[[1]]$reads <- reads

  # Store names in results[[1]], for visualization purposes
  results[[1]]$names <- names(cons)

  # Store positions with iSNVs in results[[1]], for analysis purposes
  results[[1]]$all_pos <- all_pos


  return(results)
}






