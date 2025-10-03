instant.model <- list(
  type = "spline",
  params = list(knots=c(0,0.1), slopes=10))

lagged.model <- list(
  type = "spline",
  params = list(knots=c(0,2,2.1), slopes=c(0,10))
)

curved.model <- list(
  type = "exp",
  params = list(d=1.5)
)

pconvex.model <- list(
  type = "spline",
  params = list(knots=c(0,2,4), slopes=c(0.1,0.4))
)

effect_curve <- function(x, type, params) {
  
  p <- params
  
  if (type=="exp") {
    y <- ifelse(x>6, 1, (1-exp(-x/p$d)) / (1-exp(-6/p$d)))
  }
  
  if (type=="spline") {
    
    if ((length(p$knots)-1)!=length(p$slopes)) {
      stop("Length of `knots` must equal length of `slopes` plus one")
    }
    
    y <- p$slopes[1] * pmax(0,x)
    
    for (i in 2:length(p$knots)) {
      
      if (i<length(p$knots)) {
        y <- y + (p$slopes[i] - p$slopes[i-1]) * pmax(0,x-p$knots[i])
      } else {
        y <- y - p$slopes[i-1] * pmax(0,x-p$knots[i])
      }
      
    }
    
  }
  
  if (type=="parabola") {
    
    y <- ifelse(x>6, 1, (p$a * x^2) + (p$b * x) + p$c)
    
  }
  
  if (type=="non-monotonic") {
    
    y <- ifelse(x>6, 0.75, (-1/16*x^2)+(1/2*x))
    
  }
  
  if (!(type %in% c("exp", "spline", "parabola", "non-monotonic"))) {
    stop("type must be in c('exp', 'spline', 'parabola', 'non-monotonic')")
  }
  
  return (y)
  
}