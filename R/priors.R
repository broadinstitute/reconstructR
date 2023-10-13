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

### LOG Priors
prior_mu <- function(x, mins, maxs){
  norm <- pexp(maxs$mu, 0.001) - pexp(mins$mu, 0.001)
  if(x < mins$mu){
    -Inf
  }else if(x < maxs$mu){
    log(dexp(x, 0.001) / norm)
  }else{
    -Inf
  }
}

prior_e <- function(x, mins, maxs){
  norm <- pexp(maxs$e, 0.001) - pexp(mins$e, 0.001)
  if(x < mins$e){
    -Inf
  }else if(x < maxs$e){
    log(dexp(x, 0.001) / norm)
  }else{
    -Inf
  }
}

prior_p <- function(x, mins, maxs){
  norm <- pexp(maxs$p, 0.001) - pexp(mins$p, 0.001)
  if(x < mins$p){
    -Inf
  }else if(x < maxs$p){
    log(dexp(x, 0.001) / norm)
  }else{
    -Inf
  }
}

prior_w <- function(x, mins, maxs){
  norm <- 1/(maxs$w - mins$w)
  if(x < mins$w){
    -Inf
  }else if(x < maxs$w){
    log(norm)
  }else{
    -Inf
  }
}

prior_gamma <- function(x, mins, maxs){
  norm <- pexp(maxs$gamma, 0.001) - pexp(mins$gamma, 0.001)
  if(x < mins$gamma){
    -Inf
  }else if(x < maxs$gamma){
    log(dexp(x, 0.001) / norm)
  }else{
    -Inf
  }
}


#### NOTE THE PRIOR HERE, depends a lot on what we're doing
prior_alpha <- function(x, mins, maxs){
  dbeta(x, 1, 1e10, log = T)
}

prior_tau_T <- function(x, mins, maxs){
  norm <- pexp(maxs$tau_T, 0.001) - pexp(mins$tau_T, 0.001)
  if(x < mins$tau_T){
    -Inf
  }else if(x < maxs$tau_T){
    log(dexp(x, 0.001) / norm)
  }else{
    -Inf
  }
}
prior_tau_E <- function(x, mins, maxs){
  norm <- pexp(maxs$tau_E, 0.001) - pexp(mins$tau_E, 0.001)
  if(x < mins$tau_E){
    -Inf
  }else if(x < maxs$tau_E){
    log(dexp(x, 0.001) / norm)
  }else{
    -Inf
  }
}
prior_tau_I <- function(x, mins, maxs){
  norm <- pexp(maxs$tau_I, 0.001) - pexp(mins$tau_I, 0.001)
  if(x < mins$tau_I){
    -Inf
  }else if(x < maxs$tau_I){
    log(dexp(x, 0.001) / norm)
  }else{
    -Inf
  }
}
prior_var_T <- function(x, mins, maxs){
  norm <- pexp(maxs$var_T, 0.001) - pexp(mins$var_T, 0.001)
  if(x < mins$var_T){
    -Inf
  }else if(x < maxs$var_T){
    log(dexp(x, 0.001) / norm)
  }else{
    -Inf
  }
}
prior_var_E <- function(x, mins, maxs){
  norm <- pexp(maxs$var_E, 0.001) - pexp(mins$var_E, 0.001)
  if(x < mins$var_E){
    -Inf
  }else if(x < maxs$var_E){
    log(dexp(x, 0.001) / norm)
  }else{
    -Inf
  }
}
prior_var_I <- function(x, mins, maxs){
  norm <- pexp(maxs$var_I, 0.001) - pexp(mins$var_I, 0.001)
  if(x < mins$var_I){
    -Inf
  }else if(x < maxs$var_I){
    log(dexp(x, 0.001) / norm)
  }else{
    -Inf
  }
}

prior <- function(param, x, mins, maxs){
  if(param == "mu"){
    prior_mu(x, mins, maxs)
  }else if(param == "e"){
    prior_e(x, mins, maxs)
  }else if(param == "p"){
    prior_p(x, mins, maxs)
  }else if(param == "w"){
    prior_w(x, mins, maxs)
  }else if(param == "gamma"){
    prior_gamma(x, mins, maxs)
  }else if(param == "alpha"){
    prior_alpha(x, mins, maxs)
  }else if(param == "tau_T"){
    prior_tau_T(x, mins, maxs)
  }else if(param == "tau_E"){
    prior_tau_E(x, mins, maxs)
  }else if(param == "tau_I"){
    prior_tau_I(x, mins, maxs)
  }else if(param == "var_T"){
    prior_var_T(x, mins, maxs)
  }else if(param == "var_E"){
    prior_var_E(x, mins, maxs)
  }else if(param == "var_I"){
    prior_var_I(x, mins, maxs)
  }
}



