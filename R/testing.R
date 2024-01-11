library(ctsmr)
library(R.matlab)

## LOAD DATA

w0 = 0.1;

mod_t = linspace(0,100,500)
mod_y = sin(w0*mod_t-pi/2)

val_t = linspace(100,150, 200)
val_y = sin(w0*val_t-pi/2)


plot(mod_t, mod_y)
lines(val_t, val_y)

moddata = data.frame(
  t = mod_t,
  Y = mod_y 
)

valdata = data.frame(
  t = val_t,
  Y = val_y 
)


## CTSM-R SETUP

model <- ctsm()

model$addSystem( dX1 ~ ( -w*X2 - 0.01*X1*( X1^2 + X2^2 - 1 ) )*dt + exp(p11)*dw1 )
model$addSystem( dX2 ~ (  w*X1 - 0.01*X2*( X1^2 + X2^2 - 1 ) )*dt + exp(p22)*dw2 )

model$addObs(Y ~ X2 )

model$setVariance(Y ~ exp(e11))

model$setParameter(X1 = c(init = -1, lb = -5, ub = 5))
model$setParameter(X2 = c(init = 0, lb = -5, ub = 5))

model$setParameter(w = c(init = 0.05, lb = 0, ub = 10))

model$setParameter(p11 = c(init = 1, lb = -9, ub = 10))
model$setParameter(p22 = c(init = 1, lb = -9, ub = 10))

model$setParameter(e11 = c(init = -1, lb = -9, ub = 10))

## MODEL DATA

fit = model$estimate(moddata)
summary(fit, extended=TRUE)

Pred <- predict(fit, n.ahead = 1)

plot(moddata$t, moddata$Y, type = 'l')
lines(moddata$t, Pred[[1]]$output$pred$Y, col="red")

residual = moddata$Y - Pred[[1]]$output$pred$Y

plot(moddata$t, residual, type = 'l')

## VALIDATION DATA

predlags = 1

Pred <- predict(fit, newdata = valdata, n.ahead = predlags)


plot(valdata$t, valdata$Y, type = 'l')
plot(valdata$t, Pred[[1]]$output$pred$Y, col="red")

residual = valdata$Y - Pred$output$pred$Y

plot(valdata$t, residual, type = 'l')


acf(residual)


