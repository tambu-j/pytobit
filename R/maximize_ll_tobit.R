library(nloptr)
library(stats)
library(survreg)

data('tobin')
plot(dlogis(seq(-5,5,0.1),0,1))
plot(plogis(seq(-5,5,0.1),0,1))

# steps:-
# we have the log likelihood:-
# ll = sum(log f(z) - log(scale)) + sum(log(CDF))
# z = (y - xTb) / scale
# start without vectorize b? so function is logistic with mean = 0 and scale = scale (mean is 0 as we 
# are assuming scale * epsiilon_i is distributed as logistic with mean = 0 and scale = scale)
y <- tobin$durable
x1 <- tobin$age
x2 <- tobin$quant
censored <- y <= 0

f_zi <- function(y_i, x1_i, x2_i, b0, b1, b2, l_scale) {
  z_i <- (y_i - (b0 + b1 * x1_i + b2 * x2_i))/l_scale
  return(z_i)
}

zi_dlogis <- function(y_i, x1_i, x2_i, b0, b1, b2, l_scale) {
  return(dlogis(f_zi(y_i, x1_i, x2_i, b0, b1, b2, l_scale), 0, 1))
}

zi_cdf <- function(y_i, x1_i, x2_i, b0, b1, b2, l_scale) {
  return(plogis(f_zi(y_i, x1_i, x2_i, b0, b1, b2, l_scale), 0, 1))
}

ll <- function(y, x1, x2, censored, b0, b1, b2, l_scale) {
  res <- c()
  for (i in 1:length(y)) {
    if (censored[i])
      res <- append(res, log(zi_cdf(y[i], x1[i], x2[i], b0, b1, b2, l_scale)))
    else
      res <- append(res, log(zi_dlogis(y[i], x1[i], x2[i], b0, b1, b2, l_scale)) - log(l_scale))
  }
  return (sum(res))
}

ll_wrap <- function(x) {
  b0 = x[1]
  b1 = x[2]
  b2 = x[3]
  l_scale = x[4]
  return(-ll(y,x1,x2,censored,b0,b1,b2,l_scale))
}

ll_wrap(cbind(b0,b1,b2,l_scale))

lbfgs(cbind(b0,b1,b2,l_scale), ll_wrap)

f_v_z <- function(Y, X, betas, l_scale) {
  z <- (Y - X %*% betas)/l_scale
  return(z)
}

z_v_dlogis <- function(z) {
  return(dlogis(z, 0, 1))
}

z_v_cdf <- function(z) {
  return(plogis(z, 0, 1))
}

ll_v <- function(Y, X, censored, betas, l_scale) {
  z <- f_v_z(Y, X, betas, l_scale)
  res <- c()
  res[censored] <- log(z_v_cdf(z[censored]))
  res[!censored] <- log(z_v_dlogis(z[!censored])) - log(l_scale)
  return (sum(res))
}


# OLS estimate as a starting point
ols_estimate <- as.vector(coefficients(lm('durable ~ age + quant', data=tobin)))

# scale param is key for algo to converge - mustn't be too low.
y_mean <- sum(y)/length(y)
y_var  <- sum((y-y_mean)^2)/length(y)
y_coef <- c(y_mean, y_var/3.2)
y_coef <- c(y_coef[1], log(4*y_coef[2])/2 )   # log(2*sqrt(variance)) = log(4*var)/2

# now fit an intercept only
Y = matrix(y)
X = matrix(c(rep(1,length(x1)),x1,x2), ncol=3)
betas = matrix(c(b0,b1,b2))
Y - X %*% betas

ll_v_wrap <- function(x) {
  betas = matrix(x[1:length(x)-1])
  l_scale = x[length(x)]
  return(-ll_v(Y, X, censored, betas, l_scale))
}

X = matrix(c(rep(1,length(x1))), ncol=1)
initial_fit <- lbfgs(y_coef, ll_v_wrap)
scale_guess <- initial_fit$par[2]
X = matrix(c(rep(1,length(x1)),x1,x2), ncol=3)
lbfgs(c(ols_estimate,scale_guess), ll_v_wrap)
