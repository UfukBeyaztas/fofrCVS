penalty_path <- function(delta_hat, precision_delta,
                         group_size, number_groups,
                         min_ratio = 1e-3, length.out = 31) {
  gradient_at_zero <- 2 * as.numeric(precision_delta %*% delta_hat)
  group_norms <- numeric(number_groups)
  
  for (k in seq_len(number_groups)) {
    indices <- ((k - 1) * group_size + 1):(k * group_size)
    group_norms[k] <- sqrt(sum(gradient_at_zero[indices]^2))
  }
  
  zeta_max <- max(group_norms)
  
  if (!is.finite(zeta_max) || zeta_max <= .Machine$double.eps) {
    return(0)
  }
  
  exp(seq(log(zeta_max), log(zeta_max * min_ratio), length.out = length.out))
}
