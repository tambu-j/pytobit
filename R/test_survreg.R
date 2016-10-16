library(survreg)

tfit <- survreg(Surv(durable, durable>0, type='left') ~age + quant,
                data=tobin, dist='gaussian')
predict(tfit,type="response")
summary(tfit)

tfit_logistic <- survreg(Surv(durable, durable>0, type='left') ~age + quant,
                data=tobin, dist='logistic')
predict(tfit_logistic,type="response")
summary(tfit_logistic)

