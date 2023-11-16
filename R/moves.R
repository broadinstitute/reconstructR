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

### MCMC moves


# General update function for epi parameters
update <- function(l_old, param, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, vars, sum_beta_CDF){
  l_new <- l_old

  l_new[[param]] <- l_old[[param]] + rnorm(1, 0, vars[[param]])

  if(l_new[[param]] < mins[[param]] | l_new[[param]] > maxs[[param]]){
    return(l_old)
  }else{

    l_new$e_lik <- get_e_lik(l_new, mins, maxs, prior_params)
    delta <- l_new$e_lik - l_old$e_lik

    if(param == "mu" | param == "e" | param == "p" | param == "w"){
      if(l_new$L > 0){
        l_new$g_lik <- sapply(1:l_new$N, get_g_lik, l = l_new, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs)
        if(l_new$L == 1){
          l_new$g_lik <- matrix(l_new$g_lik, nrow = 1)
        }
      }
      l_new$u_lik <- sapply(1:l_new$N, get_u_lik, l = l_new, sum_depths = sum_depths, filters = filters, sum_beta_CDF = sum_beta_CDF)
      delta <- delta + sum(l_new$g_lik) - sum(l_old$g_lik) + sum(l_new$u_lik) - sum(l_old$u_lik)
    }

    if(log(runif(1)) > delta){
      return(l_old)
    }else{
      return(l_new)
    }
  }
}

# Update proportions of a2 at end of growth phase, for all sites simultaneously in a given host
update_x <- function(l_old, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, vars, sum_beta_CDF){
  l_new <- l_old
  j <- sample(2:l_new$N, 1) # should we also do this for index case? probably not.

  # h is the ancestor of j
  h <- l_new$anc[j]

  # Propose new x values for host j
  l_new$x[[j]] <- l_new$x[[j]] + rnorm(l_new$L, 0, vars$x)

  # Negative values get re-assigned to 0
  l_new$x[[j]][l_new$x[[j]] < 0] <- 0

  # Values above 1 get re-assigned to 1
  l_new$x[[j]][l_new$x[[j]] > 1] <- 1

  # Log probability of proposing old state, starting in new state
  p_new_old <- log(
    (l_old$x[[j]] == 0) * pnorm(0, l_new$x[[j]], vars$x) +
      (l_old$x[[j]] > 0 & l_old$x[[j]] < 1)*dnorm(l_old$x[[j]], l_new$x[[j]], vars$x) +
      (l_old$x[[j]] == 1) * pnorm(1, l_new$x[[j]], vars$x, lower.tail = F)
  )

  # Log probability of proposing new state, starting in old state
  p_old_new <- log(
    (l_new$x[[j]] == 0) * pnorm(0, l_old$x[[j]], vars$x) +
      (l_new$x[[j]] > 0 & l_new$x[[j]] < 1)*dnorm(l_new$x[[j]], l_old$x[[j]], vars$x) +
      (l_new$x[[j]] == 1) * pnorm(1, l_old$x[[j]], vars$x, lower.tail = F)
  )

  # New genomic log likelihood (only changes for h and j)
  l_new$g_lik[,h] <- get_g_lik(h, l = l_new, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs)
  l_new$g_lik[,j] <- get_g_lik(j, l = l_new, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs)

  # Which proposals get rejected?
  rejections <- which(log(runif(l_new$L)) > l_new$g_lik[,h] + l_new$g_lik[,j] - l_old$g_lik[, h] - l_old$g_lik[, j] + p_new_old - p_old_new)

  #print(L - length(rejections))

  # If rejected, re-assign the old value for x and for the genetic likelihood
  l_new$x[[j]][rejections] <- l_old$x[[j]][rejections]
  l_new$g_lik[rejections, h] <- l_old$g_lik[rejections, h]
  l_new$g_lik[rejections, j] <- l_old$g_lik[rejections, j]

  return(l_new)

}

# Update t_E
update_t_E <- function(l_old, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF){

  l_new <- l_old
  j <- sample(2:l_new$N, 1)

  # h is the ancestor of j
  h <- l_new$anc[j]

  # Get max time for t_E for individual j
  max_t_E <- l_new$t_I[j]

  # Get min time for t_E for individual j
  if(h==1){
    min_t_E <- max_t_E - 5
  }else{
    min_t_E <- l_new$t_I[l_new$anc[j]]
  }

  # Assign uniformly-drawn new value for t_E
  l_new$t_E[j] <- runif(1, min_t_E, max_t_E)

  # New epi likelihood
  l_new$e_lik <- get_e_lik(l_new, mins, maxs, prior_params)

  # New genomic likelihood for h and j
  if(l_new$L > 0){
    l_new$g_lik[,h] <- get_g_lik(h, l = l_new, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs)
    l_new$g_lik[,j] <- get_g_lik(j, l = l_new, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs)
  }
  l_new$u_lik[h] <- get_u_lik(h, l_new, sum_depths = sum_depths, filters = filters, sum_beta_CDF = sum_beta_CDF)
  l_new$u_lik[j] <- get_u_lik(j, l_new, sum_depths = sum_depths, filters = filters, sum_beta_CDF = sum_beta_CDF)

  if(log(runif(1)) > sum(l_new$g_lik) - sum(l_old$g_lik) + sum(l_new$u_lik) - sum(l_old$u_lik) + l_new$e_lik - l_old$e_lik){
    return(l_old)
  }else{
    return(l_new)
  }
}

# Update t_I
update_t_I <- function(l_old, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF){

  l_new <- l_old
  j <- sample(2:l_new$N, 1)

  # Get min time for t_I for individuals 2:N
  min_t_I <- l_new$t_E[j]

  # Get max time for t_I for individuals 2:N
  max_t_I <- get_max_t_I(j, l_new)

  # Assign new values via uniform draws
  l_new$t_I[j] <- runif(1, min_t_I, max_t_I)

  # Get new epi likelihood
  l_new$e_lik <- get_e_lik(l_new, mins, maxs, prior_params)

  # New genomic likelihood for j
  if(l_new$L > 0){
    l_new$g_lik[,j] <- get_g_lik(j, l = l_new, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs)
  }
  l_new$u_lik[j] <- get_u_lik(j, l = l_new, sum_depths = sum_depths, filters = filters, sum_beta_CDF = sum_beta_CDF)

  if(log(runif(1)) > sum(l_new$g_lik) - sum(l_old$g_lik) + sum(l_new$u_lik) - sum(l_old$u_lik) + l_new$e_lik - l_old$e_lik){
    return(l_old)
  }else{
    return(l_new)
  }
}



# Update ancestor
update_ancestor <- function(l_old, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF){

  l_new <- l_old

  # resample the ancestor of j
  j <- sample(2:l_new$N, 1)

  # h_old is the OLD ancestor of j
  h_old <- l_new$anc[j]

  # h_new is the NEW ancestor of j
  choices <- which(l_new$t_I < l_new$t_E[j])
  h_new <- ifelse(length(choices) > 1, sample(choices, 1), choices)

  l_new$anc[j] <- h_new

  ## New idea: if no downstream cases, flip sites at which no read depth in j with prob 1/2
  if(!(j %in% l_new$anc)){
    no_depth <- which(Rfast::rowsums(reads[[j]]) == 0 & abs(l_new$x[[h_new]] - l_new$x[[h_old]]) > 0.5)
    to_flip <- no_depth
    l_new$x[[j]][to_flip] <- 1 - l_new$x[[j]][to_flip]
  }

  if(l_new$L > 0){
    l_new$g_lik[,c(j, h_old, h_new)] <- sapply(c(j, h_old, h_new), get_g_lik, l = l_new, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs)
  }
  l_new$u_lik[c(j, h_old, h_new)] <- sapply(c(j, h_old, h_new), get_u_lik, l_new, sum_depths = sum_depths, filters = filters, sum_beta_CDF = sum_beta_CDF)
  l_new$e_lik <- get_e_lik(l_new, mins, maxs, prior_params)

  if(log(runif(1)) > sum(l_new$g_lik) - sum(l_old$g_lik) + sum(l_new$u_lik) - sum(l_old$u_lik) + l_new$e_lik - l_old$e_lik){
    return(l_old)
  }else{
    return(l_new)
  }
}

# Update ancestor, plus x (in new/old ancestor) at sites with consensus change
update_ancestor_x <- function(l_old, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF){

  l_new <- l_old

  # resample the ancestor of j
  j <- sample(2:l_new$N, 1)

  # h_old is the OLD ancestor of j
  h_old <- l_new$anc[j]

  # h_new is the NEW ancestor of j
  choices <- which(l_new$t_I < l_new$t_E[j])
  h_new <- ifelse(length(choices) > 1, sample(choices, 1), choices)

  l_new$anc[j] <- h_new

  ## Which sites exhibit consensus change?
  # cons_change <- which(
  #   ((rowMaxs(reads[[j]]) != rowMaxs(reads[[h_old]])) & (Rfast::rowsums(reads[[j]]) != 0) & (rowsums(reads[[h_old]]) != 0)) |
  #     ((rowMaxs(reads[[j]]) != rowMaxs(reads[[h_new]])) & (rowsums(reads[[j]]) != 0) & (rowsums(reads[[h_new]]) != 0))
  # )

  cons_change <- which(
    (abs(l_new$x[[j]] - l_new$x[[h_old]]) > 0.5) | (abs(l_new$x[[j]] - l_new$x[[h_new]]) > 0.5)
  )

  # Swap x at these sites
  if(length(cons_change) > 0 & h_old != 1 & h_new != 1){
    l_new$x[[h_new]][cons_change] <- l_old$x[[h_old]][cons_change]
    l_new$x[[h_old]][cons_change] <- l_old$x[[h_new]][cons_change]
  }

  # Ancestors' ancestors, whose genomic likelihoods also need to be updated
  g_old <- l_new$anc[h_old]
  g_new <- l_new$anc[h_new]

  to_update <- c(j, h_old, h_new, g_old, g_new)
  to_update <- to_update[!is.na(to_update)]

  l_new$g_lik[,to_update] <- sapply(to_update, get_g_lik, l = l_new, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs)
  l_new$u_lik[to_update] <- sapply(to_update, get_u_lik, l_new, sum_depths = sum_depths, filters = filters, sum_beta_CDF = sum_beta_CDF)
  l_new$e_lik <- get_e_lik(l_new, mins, maxs, prior_params)

  if(log(runif(1)) > sum(l_new$g_lik) - sum(l_old$g_lik) + sum(l_new$u_lik) - sum(l_old$u_lik) + l_new$e_lik - l_old$e_lik){
    return(l_old)
  }else{
    #print("great success")
    return(l_new)
  }
}

# Update ancestor, plus x at all sites . NEED TO ADJUST G_LIK for ancestors of h_old and h_new if we want to revive this move
# update_ancestor_all_x <- function(l_old){
#
#   l_new <- l_old
#
#   # resample the ancestor of j
#   j <- sample(2:N, 1)
#
#   # h_old is the OLD ancestor of j
#   h_old <- l_new$anc[j]
#
#   # h_new is the NEW ancestor of j
#   choices <- which(l_new$t_I < l_new$t_E[j])
#   h_new <- ifelse(length(choices) > 1, sample(choices, 1), choices)
#
#   l_new$anc[j] <- h_new
#
#   # Swap x at all sites
#   if(h_old != 1 & h_new != 1){
#     l_new$x[[h_new]] <- l_old$x[[h_old]]
#     l_new$x[[h_old]] <- l_old$x[[h_new]]
#   }
#
#
#
#   l_new$g_lik[,c(j, h_old, h_new)] <- sapply(c(j, h_old, h_new), get_g_lik, l_new, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs)
#   l_new$u_lik[c(j, h_old, h_new)] <- sapply(c(j, h_old, h_new), get_u_lik, l_new, sum_depths = sum_depths)
#   l_new$e_lik <- get_e_lik(l_new)
#
#   if(log(runif(1)) > sum(l_new$g_lik) - sum(l_old$g_lik) + sum(l_new$u_lik) - sum(l_old$u_lik) + l_new$e_lik - l_old$e_lik){
#     return(l_old)
#   }else{
#     print("new move accepteddd")
#     return(l_new)
#   }
# }



## Update t_E, t_I, and ancestor
## Also idea: pick the new ancestor first, then pick the new t_E, t_I
### Also also idea: pick new x, based on original read data
update_epi_ancestor <- function(l_old, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF){
  l_new <- l_old
  # resample the ancestor of j
  j <- sample(2:l_new$N, 1)

  # h_old is the OLD ancestor of j
  h_old <- l_new$anc[j]

  # Pick new ancestor to be anyone
  choices <- 1:l_new$N
  h_new <- sample(choices, 1)
  l_new$anc[j] <- h_new

  # Maximum time of becoming infectious for j is min of t_E's for people infected by j and l_new$t_test[j]
  max_t_I <- min(c(
    l_new$t_E[which(l_new$anc == j)],
    l_new$t_test[j]
  ))

  # Old minimum time of exposure for j is the time h_old becomes infectious
  if(h_old==1){
    min_t_E_old <- max_t_I - 10
  }else{
    min_t_E_old <- l_new$t_I[h_old]
  }

  # New minimum time of exposure for j is the time h_new becomes infectious
  if(h_new==1){
    min_t_E_new <- max_t_I - 10
  }else{
    min_t_E_new <- l_new$t_I[h_new]
  }


  if(min_t_E_new >= max_t_I){
    return(l_old)
  }else{
    # Draw new t_E, t_I
    draws <- runif(2, min_t_E_new, max_t_I)
    l_new$t_E[j] <- min(draws)
    l_new$t_I[j] <- max(draws)

    # Draw new x, beta from read data
    #l_new$x[[j]] <- rbeta(L, reads[[j]][cbind(1:L, l_new$a2[[j]])]+1, rowsums(reads[[j]]) - reads[[j]][cbind(1:L, l_new$a2[[j]])]+1)

    # Probability density of going from new state to old state
    p_new_old <- -2*log(max_t_I - min_t_E_old) #+ sum(dbeta(l_old$x[[j]], reads[[j]][cbind(1:L, l_new$a2[[j]])]+1, rowsums(reads[[j]]) - reads[[j]][cbind(1:L, l_new$a2[[j]])]+1, log = T))

    # Probability density of going from old state to new state
    p_old_new <- -2*log(max_t_I - min_t_E_new) #+ sum(dbeta(l_new$x[[j]], reads[[j]][cbind(1:L, l_new$a2[[j]])]+1, rowsums(reads[[j]]) - reads[[j]][cbind(1:L, l_new$a2[[j]])]+1, log = T))

    if(l_new$L > 0){
      l_new$g_lik[,c(j, h_old, h_new)] <- sapply(c(j, h_old, h_new), get_g_lik, l = l_new, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs)
    }
    l_new$u_lik[c(j, h_old, h_new)] <- sapply(c(j, h_old, h_new), get_u_lik, l_new, sum_depths = sum_depths, filters = filters, sum_beta_CDF = sum_beta_CDF)
    l_new$e_lik <- get_e_lik(l_new, mins, maxs, prior_params)

    if(log(runif(1)) >
       sum(l_new$g_lik) - sum(l_old$g_lik) +
       sum(l_new$u_lik) - sum(l_old$u_lik) +
       l_new$e_lik - l_old$e_lik +
       p_new_old - p_old_new
    ){
      return(l_old)
    }else{
      #print("based")
      return(l_new)
    }
  }
}



# Update t_E and ancestor for a host, simultaneously
update_t_E_ancestor <- function(l_old, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF){
  l_new <- l_old
  # resample the ancestor of j
  j <- sample(2:l_new$N, 1)

  # h_old is the OLD ancestor of j
  h_old <- l_new$anc[j]

  # Get max time for t_E for individual j
  max_t_E <- l_new$t_I[j]

  # Get min time for t_E for individual j
  if(h_old==1){
    min_t_E_old <- max_t_E - 5
  }else{
    min_t_E_old <- l_new$t_I[h_old]
  }

  old_length <- max_t_E - min_t_E_old

  # Assign uniformly-drawn new value for t_E
  l_new$t_E[j] <- runif(1, min_t_E_old, max_t_E)



  # h_new is the NEW ancestor of j
  old_choices <- which(l_old$t_I < l_old$t_E[j])
  new_choices <- which(l_new$t_I < l_new$t_E[j])
  h_new <- ifelse(length(new_choices) > 1, sample(new_choices, 1), new_choices)

  l_new$anc[j] <- h_new

  # Get min time for t_E for individual j
  if(h_new==1){
    min_t_E_new <- max_t_E - 5
  }else{
    min_t_E_new <- l_new$t_I[h_new]
  }

  new_length <- max_t_E - min_t_E_new

  ## New idea: if no downstream cases, flip sites at which no read depth in j
  if(!(j %in% l_new$anc)){
    no_depth <- which(Rfast::rowsums(reads[[j]]) == 0 & abs(l_new$x[[h_new]] - l_new$x[[h_old]]) > 0.5)
    to_flip <- no_depth
    l_new$x[[j]][to_flip] <- 1 - l_new$x[[j]][to_flip]
  }

  if(l_new$L > 0){
    l_new$g_lik[,c(j, h_old, h_new)] <- sapply(c(j, h_old, h_new), get_g_lik, l = l_new, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs)
  }
  l_new$u_lik[c(j, h_old, h_new)] <- sapply(c(j, h_old, h_new), get_u_lik, l_new, sum_depths = sum_depths, filters = filters, sum_beta_CDF = sum_beta_CDF)
  l_new$e_lik <- get_e_lik(l_new, mins, maxs, prior_params)

  if(log(runif(1)) >
     sum(l_new$g_lik) - sum(l_old$g_lik) +
     sum(l_new$u_lik) - sum(l_old$u_lik) +
     l_new$e_lik - l_old$e_lik +
     log(length(new_choices)) + log(old_length) - log(length(old_choices)) - log(new_length)
  ){
    return(l_old)
  }else{
    return(l_new)
  }
}


# Swap infector/infectee
update_swap <- function(l_old, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF){
  l_new <- l_old
  # resample the ancestor of j
  if(length(which(l_new$anc != 1)) > 1){
    j <- sample(which(l_new$anc != 1), 1)
  }else{
    j <- which(l_new$anc != 1)
  }

  if(length(j) == 1){
    # h is the ancestor of j
    h <- l_new$anc[j]

    # g is the ancestor of h
    g <- l_new$anc[h]

    l_new$anc[j] <- g
    l_new$anc[h] <- j

    # Swap t_E and t_I as well
    l_new$t_E[c(h,j)] <- l_new$t_E[c(j,h)]
    l_new$t_I[c(h,j)] <- l_new$t_I[c(j,h)]

    # Swap bottleneck sizes
    l_new$A[c(h,j)] <- l_new$A[c(j,h)]

    l_new$e_lik <- get_e_lik(l_new, mins, maxs, prior_params)

    if(l_new$e_lik == -Inf){
      return(l_old)
    }else{
      if(l_new$L > 0){
        l_new$g_lik[,c(g,h,j)] <- sapply(c(g,h,j), get_g_lik, l = l_new, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs)
      }
      l_new$u_lik[c(g,h,j)] <- sapply(c(g,h,j), get_u_lik, l_new, sum_depths = sum_depths, filters = filters, sum_beta_CDF = sum_beta_CDF)

      if(log(runif(1)) > sum(l_new$g_lik) - sum(l_old$g_lik) + sum(l_new$u_lik) - sum(l_old$u_lik) + l_new$e_lik - l_old$e_lik){
        return(l_old)
      }else{
        return(l_new)
      }
    }
  }else{
    return(l_old)
  }
}

# Swap infector/infectee + kids
update_swap_kids <- function(l_old, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF){
  l_new <- l_old
  # resample the ancestor of j
  if(length(which(l_new$anc != 1)) > 1){
    j <- sample(which(l_new$anc != 1), 1)
  }else{
    j <- which(l_new$anc != 1)
  }

  if(length(j) == 1){
    # h is the ancestor of j
    h <- l_new$anc[j]

    # g is the ancestor of h
    g <- l_new$anc[h]

    l_new$anc[j] <- g
    l_new$anc[h] <- j

    # Swap t_E and t_I as well
    l_new$t_E[c(h,j)] <- l_new$t_E[c(j,h)]
    l_new$t_I[c(h,j)] <- l_new$t_I[c(j,h)]

    # Swap bottleneck sizes
    l_new$A[c(h,j)] <- l_new$A[c(j,h)]

    # Swap the kids
    l_new$anc[which(l_old$anc == j)] <- h
    l_new$anc[which(l_new$anc == h)] <- j

    l_new$e_lik <- get_e_lik(l_new, mins, maxs, prior_params)

    if(l_new$e_lik == -Inf){
      return(l_old)
    }else{
      if(l_new$L > 0){
        l_new$g_lik[,c(g,h,j)] <- sapply(c(g,h,j), get_g_lik, l = l_new, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs)
      }
      l_new$u_lik[c(g,h,j)] <- sapply(c(g,h,j), get_u_lik, l_new, sum_depths = sum_depths, filters = filters, sum_beta_CDF = sum_beta_CDF)

      if(log(runif(1)) > sum(l_new$g_lik) - sum(l_old$g_lik) + sum(l_new$u_lik) - sum(l_old$u_lik) + l_new$e_lik - l_old$e_lik){
        return(l_old)
      }else{
        #print("yip")
        return(l_new)
      }
    }
  }else{
    return(l_old)
  }
}


# Add or remove particle from bottleneck
update_neck <- function(l_old, reads, filters, h2f1_coefs, sum_depths, mins, maxs, prior_params, sum_beta_CDF){
  l_new <- l_old
  # change a particle in the bottleneck infecting j
  j <- sample(2:l_new$N, 1)
  h <- l_new$anc[j]

  if(runif(1) < 0.5){
    l_new$A[j] <- l_new$A[j] - 1
  }else{
    l_new$A[j] <- l_new$A[j] + 1
  }

  if(l_new$A[j] < 1){
    return(l_old)
  }else{
    if(l_new$L > 0){
      l_new$g_lik[,h] <- get_g_lik(h, l = l_new, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs)
      l_new$g_lik[,j] <- get_g_lik(j, l = l_new, reads = reads, filters = filters, h2f1_coefs = h2f1_coefs)
    }
    l_new$u_lik[h] <- get_u_lik(h, l_new, sum_depths = sum_depths, filters = filters, sum_beta_CDF = sum_beta_CDF)
    l_new$u_lik[j] <- get_u_lik(j, l_new, sum_depths = sum_depths, filters = filters, sum_beta_CDF = sum_beta_CDF)
    l_new$e_lik <- get_e_lik(l_new, mins, maxs, prior_params)

    if(log(runif(1)) > sum(l_new$g_lik) - sum(l_old$g_lik) + sum(l_new$u_lik) - sum(l_old$u_lik) + l_new$e_lik - l_old$e_lik){
      return(l_old)
    }else{
      return(l_new)
    }
  }
}








