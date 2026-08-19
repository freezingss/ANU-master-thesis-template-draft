library(gllvm)

data("microbialdata")
microbialdata
summary(microbialdata)
Y <- microbialdata$Y      
Xenv <- microbialdata$Xenv
str(Xenv)
table(Xenv$Site) 

build_long_microbial <- function(Y, Xenv, covariates, group_var = "Site") {
  N <- nrow(Y)
  Q <- ncol(Y)
  group <- factor(Xenv[[group_var]])
  long <- data.frame(
    count = as.vector(Y),
    category  = factor(rep(seq_len(Q), each = N)),
    group = factor(rep(group, times = Q)),
    log_total = rep(log(rowSums(Y)), times = Q)
  )
  long$obs <- interaction(long$group, long$category, drop = TRUE)
  for (cov in covariates) {
    long[[cov]] <- rep(Xenv[[cov]], times = Q)
  }
  long
}

long_microbial <- build_long_microbial(
  Y, Xenv,
  covariates = c("pH", "Phosp"),  
  group_var  = "Site"
)

str(long_microbial)