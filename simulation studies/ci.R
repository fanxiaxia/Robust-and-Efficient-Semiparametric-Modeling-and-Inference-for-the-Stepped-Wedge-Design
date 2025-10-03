calculate_ci <- function(est, se, conf_level = 0.95) {
  # Calculate the z-value for the given confidence level
  z <- qnorm(1 - (1 - conf_level) / 2)
  
  # Calculate the lower and upper bounds of the confidence interval
  lower_bound <- est - z * se
  upper_bound <- est + z * se
  
  ci <- c(lower = lower_bound, upper = upper_bound)
  return(ci)
}