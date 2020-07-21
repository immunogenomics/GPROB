risk <- read.delim(
  sep = "",
  text = "
RA         SLE
-1.9379420  0.83464674
0.3140075  0.26236426
-0.3915622 -0.18392284
-1.8325815 -0.08164399
-1.0216512  0.62700738
-1.5185740 -0.57856671
-1.6755777 -0.18700450
-2.0433025  0.54844506
")

prevalence <- c("RA" = 0.001, "SLE" = 0.001)

prob <- function(alpha, risk) {
  1 / (
    1 + exp(alpha - risk)
  )
}

alpha <- sapply(seq(ncol(risk)), function(i) {
  o <- optimize(
    f        = function(alpha, risk, prevalence) {
      ( mean(prob(alpha, risk)) - prevalence ) ^ 2
    },
    interval = c(-100, 100),
    risk = risk[,i],
    prevalence = prevalence[i]
  )
  o$minimum
})
alpha

prevalence[1]
mean(prob(alpha[1], risk[,1]))

prevalence[2]
mean(prob(alpha[2], risk[,2]))

alphas <- seq(-10, 20, length.out = 1000)
cost <- sapply(alphas, function(a) {
  ( mean(prob(a, risk[,1])) - prevalence[1] ) ^ 2
})
plot(alphas, log(cost))
abline(v = alpha[1], col = "red")

res1 <- nlminb(3, obj)
res1$par

res2 <- optimize(obj, c(1, 5))
res2$minimum

all.equal(res1$objective, res2$objective)

xs <- seq(0, 100, length.out = 1000)
plot(xs, obj(xs))
abline(v = res1$par, col = "blue")
abline(v = res2$minimum, col = "red", lty = 3)
