library(mbereg)
#set.seed(22677)
dat <- dgp(a0 = 0, a1 = 0.3, b = 0.5, n = 10000)
attach(dat)

b <- 0.5

m1 <- function(a) {
  (1 - a) - Tobs * (1 - z) / mean(1 - z)
}

m2 <- function(a) {
  (1 - a) - Tobs *  z / mean(z)
}

m3 <- function(a) {
  (y - b * Tobs / (1 - a)) * z
}

m4 <- function(a) {
  (y^2 - 1 - 2 * b * y * Tobs / (1 - a) + b^2 * Tobs / (1 - a)) * z
}

f1 <- function(a) {
  x <- m1(a)
  n <- length(x)
  (mean(x) < 0) * n * mean(x)^2 / var(x)
}

f2 <- function(a) {
  x <- m2(a)
  n <- length(x)
  (mean(x) < 0) * n * mean(x)^2 / var(x)
}

f3 <- function(a) {
  x <- m3(a)
  n <- length(x)
  n * mean(x)^2 / var(x)
}

f4 <- function(a) {
  x <- m4(a)
  n <- length(x)
  n * mean(x)^2 / var(x)
}

Q <- function(a) f1(a) + f2(a) + f3(a) + f4(a)
Qv <- Vectorize(Q)

a_seq <- seq(0, 0.9, 0.005)
plot(a_seq, Qv(a_seq), type = 'l')

f1v <- Vectorize(f1)
f2v <- Vectorize(f2)
f3v <- Vectorize(f3)
f4v <- Vectorize(f4)

#plot(a_seq, f1v(a_seq), type = 'l')
#plot(a_seq, f2v(a_seq), type = 'l')
#plot(a_seq, f3v(a_seq), type = 'l')
#plot(a_seq, f4v(a_seq), type = 'l')
