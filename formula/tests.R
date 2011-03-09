# A function to do "nested" F-tests

f.test.lm <- function(R.lm, F.lm) {
   SSE.R <- sum(resid(R.lm)^2)
   SSE.F <- sum(resid(F.lm)^2)
   df.num <- R.lm$df - F.lm$df
   df.den <- F.lm$df
   F <- ((SSE.R - SSE.F) / df.num) / (SSE.F / df.den)
   p.value <- 1 - pf(F, df.num, df.den)
   return(data.frame(F, df.num, df.den, p.value))
}

# A function to perform a General Linear Hypothesis Test

glt.lm <- function(model.lm, C, h=0) {
   df.num <- nrow(C)
   df.den <- model.lm$df
   design <- model.matrix(model.lm)
   beta.hat <- model.lm$coef
   sigma.sq <- sum(resid(model.lm)^2) / model.lm$df
   var.C <- sigma.sq * C %*% solve(t(design) %*% design) %*% t(C)
   F <- t(C %*% beta.hat - h) %*% (solve(var.C) %*% (C %*% beta.hat - h)) / df.num
   print(F)
   p.value <- 1 - pf(F, df.num, df.den)
   return(data.frame(F, df.num, df.den, p.value))
}

