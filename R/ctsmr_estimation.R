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

reduction_factor = 0.75

data = data.frame(
  t = ElGeneina$rain_org_t,
  tYear = ElGeneina$rain_org_t %% 1,  # TIME OF YEAR BETWEEN 0 AND 1
  obsY = ElGeneina$rain_org,
  obsZ = log(ElGeneina$rain_org/reduction_factor + 1) # LAMPERTI TRANSFORMATION 
)

moddata = data[1:382,]

valdata = data[383:480,]


## CTSM-R SETUP

model <- ctsm()

model$addSystem( dX1 ~ ( -w*X2 - mu*X1*(X1^2 + X2^2 - 1) )*dt + exp(p11)*dw1 )
model$addSystem( dX2 ~ (  w*X1 - mu*X2*(X1^2 + X2^2 - 1))*dt + exp(p22)*dw2 )
#model$addSystem( dX3 ~ ( c*w*X1 * exp(c*(1+X2)) * A*exp(-2*c) )*dt + exp(p33)*dw3 )

#model$addObs(obsZ ~ X3 )
model$addObs(obsZ ~ A*exp(-2*c)*(exp(c*(1+X2)) - 1) )


model$setVariance(obsZ ~ exp(e11))

model$setParameter(X1 = c(init = 0.9, lb = -1.5, ub = 1.5))
model$setParameter(X2 = c(init = 0.11, lb = -1.5, ub = 1.5))
model$setParameter(X3 = c(init = 0.0001, lb = -0.0, ub = 10))

model$setParameter(A = c(init = 10, lb = 0.01, ub = 50))
model$setParameter(c = c(init = 4, lb = 0.1, ub = 30))
model$setParameter(w = c(init = 6.4, lb = 1, ub = 10))
model$setParameter(mu = c(init = 0.1, lb = 0.01, ub = 0.5))


model$setParameter(p11 = c(init = 1, lb = -9, ub = 10))
model$setParameter(p22 = c(init = 1, lb = -9, ub = 10))
model$setParameter(p33 = c(init = 1, lb = -9, ub = 10))

model$setParameter(e11 = c(init = -1, lb = -9, ub = 10))

## MODEL DATA

fit = model$estimate(moddata)
summary(fit, extended=TRUE)

predlags = 1
Pred <- predict(fit, newdata = data, n.ahead = predlags)

## MODELING DATA

plot(moddata$t, moddata$obsZ, type = 'l', main = "1-Step Prediction in Transformed Domain")
lines(moddata$t, Pred$output$pred$obsZ[1:382], col="red")

plot(moddata$t, moddata$obsY, type = 'l', main = "1-Step Prediction")
lines(moddata$t, reduction_factor*( exp(Pred$output$pred$obsZ[1:382]) - 1), col="red")

residual_tf = moddata$obsZ - Pred$output$pred$obsZ[1:382]
residual = moddata$obsY - reduction_factor*( exp(Pred$output$pred$obsZ[1:382]) - 1)
writeMat(con="mod_resid_tf.mat", x=as.matrix(residual_tf))
writeMat(con="mod_resid.mat", x=as.matrix(residual))

plot(moddata$t, residual_tf, type = 'l', main = "Residual in Transformed Domain")

plot(moddata$t, residual, type = 'l', main = "Residual")

## VALIDATION DATA

plot(valdata$t, valdata$obsZ, type = 'l', main = "1-Step Prediction in Transformed Domain")
lines(valdata$t, Pred$output$pred$obsZ[383:480], col="red")

plot(valdata$t, valdata$obsY, type = 'l', main = "1-Step Prediction")
lines(valdata$t, reduction_factor*( exp(Pred$output$pred$obsZ[383:480]) - 1), col="red")

residual_tf = valdata$obsZ - Pred$output$pred$obsZ[383:480]
residual = valdata$obsY - reduction_factor*( exp(Pred$output$pred$obsZ[383:480]) - 1)
writeMat(con="val_resid_tf.mat", x=as.matrix(residual_tf))
writeMat(con="val_resid.mat", x=as.matrix(residual))

plot(valdata$t, residual_tf, type = 'l', main = "Residual in Transformed Domain")

plot(valdata$t, residual, type = 'l', main = "Residual")



