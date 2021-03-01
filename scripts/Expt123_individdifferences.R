library(tidyverse)
library(caret)
library(ppcor)
library(psychometric)

setwd("~/Dropbox/writing/articles/2019_rtPrestim/results/")

data <- read.csv(file = 'Expt123_individdifferences.csv',header=TRUE)

boxplot(data$resperr_ltm)
boxplot(data$sustained)
boxplot(data$spatial)

#ltm explained by sustained attention
fit_sust <- lm(resperr_ltm ~ sustained, data=data)
summary(fit_sust)
confint(fit_sust)

#ltm explained by spatial attention
fit_space <- lm(resperr_ltm ~ spatial, data=data)
summary(fit_space)
confint(fit_space)

#sustained vs. spatial attention
fit_sustspace <- lm(sustained ~ spatial, data=data)
summary(fit_sustspace)
confint(fit_sustspace)

#correlation
rho = cor.test(data$sustained,data$spatial,method="spearman",conf.level=0.95)
n = length(data$spatial)
round(CIr(rho$estimate,n,level=0.95),digits=2)

#partial correlation
rho = pcor.test(data$resperr_ltm,data$sustained,data$spatial,method="spearman")
CIr(rho$estimate,n,level=0.95)
rho = pcor.test(data$resperr_ltm,data$spatial,data$sustained,method="spearman")
round(rho$estimate,digits=2)
round(CIr(rho$estimate,n,level=0.95),digits=2)

#ltm explained by sustained and spatial
fit1 <- lm(resperr_ltm ~ spatial + sustained, data=data)
summary(fit1)
confint(fit1)

#test model fits
anova(fit1,fit_sust)
anova(fit1,fit_space)

library(performance)
ol <- check_outliers(fit_sust, method = "cook")
fit_sust_cook <- lm(resperr_ltm ~ sustained, data=data[ol==FALSE, ])
summary(fit_sust_cook)
ol <- check_outliers(fit_sust, method = "mahalanobis")
fit_sust_mahal <- lm(resperr_ltm ~ sustained, data=data[ol==FALSE, ])
summary(fit_sust_mahal)
ol <- check_outliers(fit_space, method = "cook")
fit_space_cook <- lm(resperr_ltm ~ spatial, data=data[ol==FALSE, ])
summary(fit_space_cook)
ol <- check_outliers(fit_space, method = "mahalanobis")
fit_space_mahal <- lm(resperr_ltm ~ spatial, data=data[ol==FALSE, ])
summary(fit_space_mahal)
ol <- check_outliers(fit1, method = "cook")
fit1_cook <- lm(resperr_ltm ~ sustained + spatial, data=data[ol==FALSE, ])
summary(fit1_cook)
ol <- check_outliers(fit1, method = "mahalanobis")
fit1_mahal <- lm(resperr_ltm ~ sustained + spatial, data=data[ol==FALSE, ])
summary(fit1_mahal)

library(gvlma)
gvlma(x=fit_sust)
gvlma(x=fit_space)
gvlma(x=fit1)


#ltm explained by sustained and spatial and wm
#fit2 <- lm(resperr_ltm ~ spatial + sustained + resperr_wm, data=data)
#summary(fit2)
#confint(fit2)

data <- read.csv(file = 'Expt123_individdifferences_z.csv',header=TRUE)

fit_sust <- lm(resperr_ltm ~ sustained, data=data)
summary(fit_sust)
confint(fit_sust)
fit_sust.diag.metrics <- augment(fit_sust)
ggplot(fit_sust.diag.metrics, aes(sustained, resperr_ltm)) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  geom_segment(aes(xend = sustained, yend = .fitted), color = "red", size = 0.3)

highleverage <- function(fit) {
  p <- length(coefficients(fit))
  n <- length(fitted(fit))
  ratio <-p/n
  plot(hatvalues(fit), main="Index Plot of Ratio")
  abline(h=c(2,3)*ratio, col="red", lty=2)
  identify(1:n, hatvalues(fit), names(hatvalues(fit)))
}
highleverage(fit)
influencePlot(fit, id.method="identify", main="Influence Plot", sub="Circle size is proportional to Cook’s distance")

fit_space <- lm(resperr_ltm ~ spatial, data=data)
summary(fit_space)
confint(fit_space)


