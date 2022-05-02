MHScoreSample <- function(N, target.str, alphabet.len, gamma,
                          start.str = sample(0:(alphabet.len - 1), str.len, replace = T)){
  str.len = length(target.str)
  curr.str = start.str
  pos.change = sample(1:str.len, N, replace = T)
  shift = sample(0:(alphabet.len - 1), N, replace = T)
  unif01 <- runif(N)
  
  return(MHScoreSampling(N, 
                         as.integer(target.str),
                         as.integer(alphabet.len),
                         gamma,
                         as.integer(start.str),
                         as.integer(pos.change),
                         as.integer(shift),
                         unif01))
}

N.batches <- function(N) {
  u <- floor(sqrt(N))
  while (N %% u != 0) u <- u - 1
  return(u)
}

MH.UC.Rcpp <- function(N, target.str, thr, alphabet.len, gamma, batch = F, 
                       start.str = sample(0:(alphabet.len - 1), str.len, replace = T)) {
  str.len = length(target.str)
  scores <- MHScoreSample(N + 50000, target.str, alphabet.len, gamma, start.str)
  p = alphabet.len^(-str.len)
  scores <- scores[-(1:50000)]
  weights = vapply(scores, function(x)  exp(-x * gamma), numeric(1))
  h = as.numeric(scores >= thr)
  if (batch){
    u = N.batches(N)
    v = N / u
    Y = vapply(1:u, 
               function(n) mean(h[((n - 1) * v + 1) : (n * v)] * weights[((n - 1) * v + 1) : (n * v)]) / mean(weights),
               numeric(1))
    #Y = vapply(1:u, 
    #           function(x) mean(h[((x-1)*v):(x*v-1)] * weights[((x-1)*v):(x*v-1)]) / mean(weights),
    #           numeric(1))
    p_est = sum(h*weights) / sum(weights)
    var_est = (v / (u - 1)) * sum((Y - p_est)^2) / N
    return(list(p = p_est, var = var_est))
  }
  else{
    return(sum(h * weights) / sum(weights))
  }
}
