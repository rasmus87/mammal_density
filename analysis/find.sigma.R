load("data/dens.DNA.tax.model.RData")

dens.DNA.tax.model

n <- nrow(dens.DNA.tax.model$data)
sigma <- sqrt(sum((predict(dens.DNA.tax.model) - dens.DNA.tax.model$data$log10_Density)^2)/(n-2))
sigma

s <- summary(dens.DNA.tax.model)
sqrt(mean(s$deviance.resid^2))
