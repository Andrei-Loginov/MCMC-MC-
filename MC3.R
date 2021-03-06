library(Rcpp)

sourceCpp('MH.cpp')

T <- TRUE
F <- FALSE

const.estimator <- function(str.len, alphabet.len, gamma, N){
  start.str = sample(0:(alphabet.len - 1), str.len, replace = T)
  scores <- MHScoreSample(N + 50000, rep(0, str.len), alphabet.len, gamma, start.str)
  scores <- scores[-(1:50000)]
  weight = vapply(scores, function(x) exp(-x*gamma), numeric(1))
  return(mean(weight))
}

N.batches <- function(N) {
  u <- floor(sqrt(N))
  while (N %% u != 0) u <- u - 1
  return(u)
}


MC3.UC.Rcpp <- function(N, target.str, thr, alphabet.len, first.gamma, last.gamma, nchains,
                                        step, batch = T, ACF = F){
  str.len = length(target.str)
  scores <- numeric(0)
  gammas <- numeric(0)
  if (nchains == 1){
    gammas = first.gamma
    scores <- MC3ScoreSample(N + 50000, target.str, alphabet.len, gammas, step)
  }
  else {
    gammas = if (nchains > 2)vapply(1:nchains, function(i) first.gamma + (i-1)*(last.gamma - first.gamma) / (nchains - 2), numeric(1)) else c(first.gamma, last.gamma)
    scores <- MC3ScoreSample(N + 50000, target.str, alphabet.len, gammas, step)
  }
  scores <- scores[-(1:50000)]
  Z = const.estimator(str.len, alphabet.len, gammas[1], N)
  q = vapply(scores, function(x)  exp(x * gammas[1]) * Z, numeric(1))
  h = as.numeric(scores >= thr)
  if (ACF){
    acf(scores)
  }
  if (batch){
    u = N.batches(N)
    v = N / u
    Y = vapply(1:u, 
               function(n) mean(h[((n - 1) * v + 1) : (n * v)] * weights[((n - 1) * v + 1) : (n * v)]) / mean(weights),
               numeric(1))
    p_est = sum(h*weights) / sum(weights)
    var_est = (v / (u - 1)) * sum((Y - p_est)^2) / N
    return(list(p = p_est, var = var_est))
    
  }
  return(mean(h / q))
}


MC3.UC.mult.Rcpp <- function(N, target.str, thr, alphabet.len, first.gamma, omega, nchains,
                        step, batch = T, ACF = F){
  str.len = length(target.str)
  
  gammas <- vapply(0:(nchains - 1), function(i) first.gamma * omega^i, numeric(1))
  scores <- MC3ScoreSample(N + 50000, target.str, alphabet.len, gammas, step)
  scores <- scores[-(1:50000)]
  Z = const.estimator(str.len, alphabet.len, gammas[1], N)
  q = vapply(scores, function(x)  exp(x * gammas[1]) * Z, numeric(1))
  h = as.numeric(scores >= thr)
  if (ACF){
    acf(scores)
  }
  if (batch){
    u = N.batches(N)
    v = N / u
    Y = vapply(1:u, 
               function(n) mean(h[((n - 1) * v + 1) : (n * v)] * weights[((n - 1) * v + 1) : (n * v)]) / mean(weights),
               numeric(1))
    p_est = sum(h*weights) / sum(weights)
    var_est = (v / (u - 1)) * sum((Y - p_est)^2) / N
    return(list(p = p_est, var = var_est))
    
  }
  return(mean(h / q))
}

