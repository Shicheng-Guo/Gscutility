
# Usage of RefFreeEWAS
# Data: 2016/09/29

require(RefFreeEWAS)
tmp1 <- lm(t(Y) ~ pred)
dim3 <- EstDimRMT(cbind(t(coef(tmp1)), t(resid(tmp1))), FALSE)$dim
test <- RefFreeEwasModel(M2B(Y), cbind(1, as.numeric(pred)-1), dim3)
testBoot <- BootRefFreeEwasModel(test, 50)
res <- summary(testBoot)
est <- res[, 2, 1, 'mean']
est.sd <- res[, 2, 1, 'sd']
Q2 <- (est/est.sd)^2
median(Q2, na.rm = TRUE)/qf(0.5, 1, dof)  # Inflation factor 
P0 <- 2*exp(pt(-sqrt(Q2), dof, log.p=T)) # 
qqplot1.pvals(P0)
