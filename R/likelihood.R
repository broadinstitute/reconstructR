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

### Likelihood functions

# Epidemiological likelihood function
get_e_lik <- function(l, mins, maxs){

  has_anc <- which(!is.na(l$anc) & (l$anc != 1))
  # Sojourn intervals
  sum(gamma_dens(l$t_test[2:l$N] - l$t_I[2:l$N], l$tau_T, l$var_T)) +

    # Exposed periods
    sum(gamma_dens(l$t_I[2:l$N] - (l$t_E[2:l$N] + (1/((l$mu/l$p )*log(l$w)))*log(1/(l$A[2:l$N] * sqrt(l$p)))), l$tau_E, l$var_E)) +

    # Does the number of virions at time t_I exceed 1/sqrt(p)?
    #exceed(l) +

    # Lambda can't be crazy large....
    # dnorm(l$mu/l$p , mean = 5, sd = 2, log = T) +

    # Number of virions can't be crazy large, or crazy small...
    #sum(dnorm(log10(l$A[2:l$N] * exp((l$t_I[2:l$N] - l$t_E[2:l$N]) * (l$mu/l$p ))), mean = 10, sd = 2, log = T)) +

    # Infectious at time of transmission
    sum(log(gamma_surv(l$t_E[has_anc], (l$t_I[l$anc[has_anc]]), mean = l$tau_I, var = l$var_I))) +

    # Size of bottleneck
    sum(dpois(l$A[2:l$N] - 1, l$gamma, log = T)) +

    # Number of index cases
    dbinom(sum(l$anc[2:l$N] == 1) - 1, l$N - 1, l$alpha, log = T) +

    # Contact tracing
    ifelse(
      l$N_ct != -1,
      get_ct_lik(l$anc, l$contacts, l$nu, l$xi, l$N_ct),
      0
    ) +

    # log priors
    prior("mu", l$mu, mins, maxs) + prior("e", l$e, mins, maxs) + prior("gamma", l$gamma, mins, maxs) + prior("alpha", l$alpha, mins, maxs) + prior("p", l$p, mins, maxs) + prior("w", l$w, mins, maxs) +
    prior("tau_T", l$tau_T, mins, maxs) + prior("tau_E", l$tau_E, mins, maxs) + prior("var_T", l$var_T, mins, maxs) + prior("var_E", l$var_E, mins, maxs)

}

# Likelihood of a singular transmission event, i to j
get_trans_lik <- function(i, j, l, h2f1_coefs){
  # Time from i's ending expo growth phase stage to transmission
  # lambda is the exponential growth rate of the virus, i.e. we have A*exp(t*l) virions at time t
  lambda <- l$mu/l$p

  # Solving for the difference in time of transmission and time to get 1/sqrt(p) virions...
  # (if index case, just set this to 3 days)
  if(i==1){
    delta_t <- 3
  }else{
    delta_t <- l$t_E[j] - (l$t_E[i] + (1/(lambda * log(l$w)))*log(1/(l$A[i] * sqrt(l$p))))
  }

  if(delta_t > 0){

    ### We want to get the fractions of alleles a1 and a2 in j that are present in i at time of transmission.
    # Call these q1 and q2
    q1 <- rep(0,l$L)
    q2 <- rep(0,l$L)

    # Where does allele a1 for j match allele a1 for i?
    a1_same <- which(l$a1[[j]] == l$a1[[i]])
    # Where does allele a1 for j match allele a2 for i?
    a1_swap <- which(l$a1[[j]] == l$a2[[i]])

    # Where does allele a2 for j match allele a2 for i?
    a2_same <- which(l$a2[[j]] == l$a2[[i]])
    # Where does allele a2 for j match allele a1 for i?
    a2_swap <- which(l$a2[[j]] == l$a1[[i]])

    q1[a1_same] <- 1 - l$x[[i]][a1_same]
    q1[a1_swap] <- l$x[[i]][a1_swap]

    q2[a2_same] <- l$x[[i]][a2_same]
    q2[a2_swap] <- 1 - l$x[[i]][a2_swap]

    # (and otherwise, 0 if no match)

    ## Now we need to evolve q1 and q2 per Jukes-Cantor (except for index case)
    ## Note, in index case, we take q1 and q2 to represent the fractions of the general population that have alleles a1 and a2, respectively
    ### NOTE we will want to come back to this and evolve at a different rate for index cases

    q1 <- (1-l$mu*delta_t)*q1 + l$mu*delta_t*(1-q1)/3
    q2 <- (1-l$mu*delta_t)*q2 + l$mu*delta_t*(1-q2)/3

    ### First, deal with the case that a2 fraction in j is 0 or 1.
    all_a1 <- which(l$x[[j]] == 0)
    all_a2 <- which(l$x[[j]] == 1)

    if(i == 1){
      # If j is an index case, our statistical model always assumes the bottleneck is homozygous.
      # Here, q1 and q2 represent the probability that an infector not in our dataset has allele a1 or a2, at 100 frequency
      # Then, we have to have no mutation in the growth phase of j

      # Log of: P(pick outside host with correct major allele) * P(every bottleneck virion features that allele) * No mutations in growth phase
      p_all_a1 <- log(q1[all_a1]) + log(1-l$mu*delta_t) + (1/sqrt(l$p))*log(1-l$p)
      p_all_a2 <- log(q2[all_a2]) + log(1-l$mu*delta_t) + (1/sqrt(l$p))*log(1-l$p)


    }else{

      # If j isn't an index case, the only way we can have all a1 or a2 is if we draw a1 or a2 (respectively) for every bottleneck virion, plus no mutation in growth phase
      p_all_a1 <- l$A[j]*log(q1[all_a1]) +  (1/sqrt(l$p))*log(1-l$p)
      p_all_a2 <- l$A[j]*log(q2[all_a2]) +  (1/sqrt(l$p))*log(1-l$p)
    }




    ### Now the case when a2 fraction is in (0,1). Here we need to split into cases where bottleneck is hetero/homozygous
    cont <- which(l$x[[j]] != 0 & l$x[[j]] != 1)

    # Let's deal with the homozygous bottleneck case first.
    # Two subcases: bottleneck is all a1, bottleneck is all a2
    if(i==1){

      # first term here: q1 is the probability we pick someone in general population with correct major allele
      # we then assume this person DOESN'T have an iSNV at the given site, so we evolve the 100% fraction by (1 - mu*delta_t).
      p_cont_homo <- (q1[cont] * (1-l$mu*delta_t) * dspecial(l$x[[j]][cont], 1, l$p) +
                        q2[cont] * (1-l$mu*delta_t) * dspecial(1 - l$x[[j]][cont], 1, l$p)) *
        (1 - (1-l$p)^(1/sqrt(l$p)))

    }else{

      p_cont_homo <- (q1[cont]^(l$A[j]) * dspecial(l$x[[j]][cont], l$A[j], l$p) +
                        q2[cont]^(l$A[j]) * dspecial(1 - l$x[[j]][cont], l$A[j], l$p)) *
        (1 - (1-l$p)^(1/sqrt(l$p)))

    }

    # Now for the heteozygous case.
    if(i==1){

      # Here we need to multiply by a factor that every choice of bottleneck virion is EITHER a1 or a2.
      # This will always be equal to (1-l$mu*delta_t) + l$mu*delta_t/3, since we assume external cases have no appreciable iSNVs

      p_cont_hetero <- (1 - l$mu*delta_t + l$mu*delta_t/3) *
        (dbinbeta(l$x[[j]][cont], 1, l$mu*delta_t/3, h2f1_coefs) * q1[cont] +
           dbinbeta(l$x[[j]][cont], 1, 1 - l$mu*delta_t, h2f1_coefs) * q2[cont])


    }else{
      # Here we need to multiply by a factor that every choice of bottleneck virion is EITHER a1 or a2
      p_cont_hetero <- (q1[cont] + q2[cont])^(l$A[j]) * dbinbeta(l$x[[j]][cont], l$A[j], q2[cont], h2f1_coefs)

    }

    p_cont <- log(p_cont_homo + p_cont_hetero)

    out <- rep(0,l$L)
    out[all_a1] <- p_all_a1
    out[all_a2] <- p_all_a2
    out[cont] <- p_cont

    return(out)
  }else{
    return(rep(-Inf, l$L))
  }


}


# Genomic likelihood function. Returns vector of genomic likelihood at position 1:L in host i
get_g_lik <- function(i, l, reads, filters, h2f1_coefs){

  # Who does the virus get passed on to?
  onward <- which(l$anc == i)

  # Log-likelihood from transmissions
  if(length(onward) >= 2){
    if(l$L > 1){
      out <- Rfast::rowsums(sapply(onward, get_trans_lik, i=i, l=l, h2f1_coefs=h2f1_coefs))
    }else{
      out <- sum(sapply(onward, get_trans_lik, i=i, l=l, h2f1_coefs=h2f1_coefs))
    }
  }else if(length(onward) == 1){
    out <- get_trans_lik(i, onward, l, h2f1_coefs=h2f1_coefs)
  }else{
    out <- rep(0,l$L)
  }


  # Log-likelihood from read data
  if(i!=1){
    # Time from i's ending expo growth phase stage to test
    # lambda is the exponential growth rate of the virus, i.e. we have A*exp(t*l) virions at time t
    lambda <- l$mu/l$p

    # Solving for the difference in time of transmission and time to get 1/sqrt(p) virions...
    delta_t <- l$t_test[i] - (l$t_E[i] + (1/(lambda * log(l$w)))*log(1/(l$A[i] * sqrt(l$p))))

    if(delta_t > 0){
      # Proportions of each of four alleles, evolved
      qs <- evolve(l$x[[i]], l$a1[[i]], l$a2[[i]], l$mu, l$e, delta_t)

      read_lik <- trunc_dmnom(reads[[i]], Rfast::rowsums(reads[[i]]), qs, filters, log = T)
      read_lik[is.nan(read_lik)] <- 0

      out <- out + read_lik
    }else{
      out <- out - Inf
    }

  }

  return(out)
}

################


# Likelihood of a singular transmission event, i to j
get_u_trans_lik <- function(i, j, l){
  # Time from i's ending expo growth phase stage to transmission
  # lambda is the exponential growth rate of the virus, i.e. we have A*exp(t*l) virions at time t
  lambda <- l$mu/l$p

  # Solving for the difference in time of transmission and time to get 1/sqrt(p) virions...
  # (if index case, just set this to 3 days)
  if(i==1){
    delta_t <- 3
  }else{
    delta_t <- l$t_E[j] - (l$t_E[i] + (1/(lambda * log(l$w)))*log(1/(l$A[i] * sqrt(l$p))))
  }

  if(delta_t > 0){

    # Per Jukes-Cantor, we compute the fraction of the major allele in i at time of transmission (starting from 100%)
    q <- (1-l$mu*delta_t)

    # The probability that we get 100% major allele in j is (approx) the probability that the bottleneck is all major
    p <- l$A[j]*log(q)

    # And, this has to happen at all K sites for which no mutation ever occurs
    return(l$K*p)
  }else{
    return(-Inf)
  }

}

# Genomic likelihood function for un-mutated positions on the genome, i.e. sites where everyone has the same allele at 100%
get_u_lik <- function(i, l, sum_depths, filters, sum_beta_CDF){

  out <- 0
  # Who does the virus get passed on to?
  onward <- which(l$anc == i)
  for (j in onward) {
    out <- out + get_u_trans_lik(i,j,l)
  }

  # Time from i's ending expo growth phase stage to transmission
  # lambda is the exponential growth rate of the virus, i.e. we have A*exp(t*l) virions at time t
  if(i == 1){
    return(out)
  }else{

    lambda <- l$mu/l$p

    # Solving for the difference in time of transmission and time to get 1/sqrt(p) virions...
    delta_t <- l$t_test[i] - (l$t_E[i] + (1/(lambda * log(l$w)))*log(1/(l$A[i] * sqrt(l$p))))

    if(delta_t > 0){
      out <- out + sum_depths[i] * log(
        p_below_af(l$p, sum_beta_CDF, sum(l$A[onward])) * (1 - (1 - l$p)^(1/sqrt(l$p))) + # probability if we do get a mutation in exp growth phase
          (1 - l$p)^(1/sqrt(l$p))
      ) +

        # ... plus indicator that quiescent term doesn't exceed the filter
        ifelse(l$mu*delta_t + l$e < filters$af, 0, -Inf)

      return(out)
    }else{
      return(-Inf)
    }
  }
}














