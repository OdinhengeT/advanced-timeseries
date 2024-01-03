library(ctsmr)
library(R.matlab)

## LOAD DATA

A = readMat("data/proj23.mat")

ElGeneina <- list()

ElGeneina$nvdi <- unlist(A$ElGeneina[1]) / 127.5 - 1
ElGeneina$nvdi_t <- unlist(A$ElGeneina[2])

ElGeneina$rain_org <- unlist(A$ElGeneina[5])
ElGeneina$rain_org_t <- unlist(A$ElGeneina[6])

plot(ElGeneina$rain_org_t, ElGeneina$rain_org)

plot(ElGeneina$nvdi_t, ElGeneina$nvdi)

data = data.frame(
  t = ElGeneina$rain_org_t,
  tYear = ElGeneina$rain_org_t %% 1,  # TIME OF YEAR BETWEEN 0 AND 1
  obsY = ElGeneina$rain_org,
  obsZ = log(ElGeneina$rain_org/0.75 + 1) # LAMPERTI TRANSFORMATION 
)

moddata = data[1:382,]

valdata = data[383:480,]


## CTSM-R SETUP

model <- ctsm()

model$addSystem( dR ~ ( mu*(r0 - R) )*dt + exp(p11)*dw1 )
model$addSystem( dTheta ~ ( w )*dt + exp(p22)*dw2 )

#model$addSystem( dZ ~ ( exp(R*(1+sin(Theta)) - Z) - 0.5 )*dt + exp(p33)*dw3 )
# ( mu*(r0-R)*(1+sin(Theta)) + R*w*cos(Theta))*

model$addObs(obsZ ~ a*(exp(R*(1+sin(Theta))) - 1.0) )

model$setVariance(obsZ ~ exp(e11))

model$setParameter(R = c(init = 6, lb = 0.01, ub = 10))
model$setParameter(Theta = c(init = -pi/2, lb = -pi, ub = pi))
model$setParameter(Z = c(init = 0, lb = -10, ub = 10))

model$setParameter(a = c(init = 1, lb = 0.1, ub = 10))
#model$setParameter(b = c(init = 0.75, lb = 0.01, ub = 1000))
model$setParameter(r0 = c(init = 5, lb = 0.005, ub = 10))
model$setParameter(mu = c(init = 0.3, lb = 0.005, ub = 0.8))
model$setParameter(w = c(init = 6.4, lb = 1, ub = 10))

model$setParameter(p11 = c(init = 1, lb = -9, ub = 10))
model$setParameter(p22 = c(init = 1, lb = -9, ub = 10))
model$setParameter(p33 = c(init = 3, lb = -9, ub = 10))

model$setParameter(e11 = c(init = -1, lb = -9, ub = 10))

## MODEL DATA

fit = model$estimate(moddata)
summary(fit, extended=TRUE)

Pred <- predict(fit, n.ahead = 1)

plot(moddata$t, moddata$obsZ, type = 'l')
lines(moddata$t, Pred[[1]]$output$pred$obsZ, col="red")

residual = moddata$obsZ - Pred[[1]]$output$pred$obsZ

plot(moddata$t, residual, type = 'l')

## VALIDATION DATA

predlags = 12

Pred <- predict(fit, newdata = data, n.ahead = predlags)

plot(valdata$t, valdata$obsZ, type = 'l')
lines( valdata$t, Pred$output$pred$obsZ[383:480], col="red")

plot(valdata$t, valdata$obsY)
lines(valdata$t, 0.75 * (exp(Pred$output$pred$obsZ[383:480]) - 1), col="red")


residual = valdata$obsY - 0.75 * (exp(Pred$output$pred$obsZ[383:480]) - 1)

plot(valdata$t, residual, type = 'l')



acf(residual)


