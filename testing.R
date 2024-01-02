library(pracma)
library(ctsmr)

N <- 1000

a <- 1.3

t <- linspace(0, 10, N)

y <- rnorm(N)

for (i in 2:100)  {
  y[i] <- y[i-1] + (t[i] - t[i-1]) * a * exp(-y[i-1]) + y[i]
}

plot(t, y)

data <- data.frame(
  t = t,
  y = y
)

model <- ctsm()

model$addSystem(dY ~ (a * exp(-Y))*dt + exp(p11)*dw1)

model$addObs(y ~ Y)

model$setVariance(y ~ exp(e11))

model$setParameter(Y = c(init = 0.01, lb = -20, ub = 20))

model$setParameter(a = c(init = 1.3, lb = -10, ub = 10))

model$setParameter(p11 = c(init = 1, lb = -10, ub = 10))
model$setParameter(e11 = c(init = -1, lb = -10, ub = 10))


fit = model$estimate(data)
summary(fit, extended=TRUE)

Pred <- predict(fit)

plot(data$t, Pred[[1]]$output$pred$y)
lines(data$t, data$y)
