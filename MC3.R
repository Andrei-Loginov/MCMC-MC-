library(Rcpp)

sourceCpp('MH.cpp')

T <- TRUE
F <- FALSE

const.estimator <- function(str.len, alphabet.len, gamma, N){
  start.str = sample(0:(alphabet.len - 1), str.len, replace = T)
  scores <- MHScoreSample(N, rep(0, str.len), alphabet.len, gamma, start.str)
  weight = vapply(scores, function(x) exp(-x*gamma), numeric(1))
  return(mean(weight))
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
  Z = const.estimator(str.len, alphabet.len, gammas[1], 2*10^5)
  q = vapply(scores, function(x)  exp(x * gammas[1]) * Z, numeric(1))
  h = as.numeric(scores >= thr)
  if (ACF){
    acf(scores)
  }
  if (batch){
    u = trunc(sqrt(N))
    v = trunc(sqrt(N))
    Y = vapply(1:u, function(x) mean(h[((x-1)*v):(x*v-1)] / q[((x-1)*v):(x*v-1)]),numeric(1))
    p_est = mean(h / q)
    var_est =  (v / (u - 1) * sum((Y - p_est)^2)) / (N)
    return(list(p = p_est, var = var_est))
  }
  return(mean(h / q))
}




