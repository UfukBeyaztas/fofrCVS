default_penalty_grid <- function(data_matrix, penalty_matrix,
                                 lower_power = -6, upper_power = 2,
                                 length.out = 31) {
  data_scale <- sum(diag(data_matrix)) / nrow(data_matrix)
  penalty_scale <- sum(diag(penalty_matrix)) / nrow(penalty_matrix)
  
  if (!is.finite(penalty_scale) || penalty_scale <= .Machine$double.eps) {
    penalty_scale <- 1
  }
  
  reference <- max(data_scale / penalty_scale, .Machine$double.eps)
  reference * 10^seq(lower_power, upper_power, length.out = length.out)
}