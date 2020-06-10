#' The bnROC Function
#'
#' This function computes the ROC curve and the AUC value
#' for a range of cutoff thresholds between [-3,3].
#'
#' @param score: a vector of test scores, continuous
#' @param status: a vector of disease statuses, binary (0/1)
#' @return A list containing important variables like AUC, TPR, and FPR.
#' @export
bnROC <- function (score, status) {
  D <- score
  Y <- status

  posIndex <- which(Y == 1)
  negIndex <- which(Y == 0)

  DY <- D[posIndex]
  DYbar <- D[negIndex]

  nY <- sum(Y)
  n <- length(D)
  nYbar <- n - nY

  posparam <- MLestimates(DY)
  negparam <- MLestimates(DYbar)

  A <- abs(posparam$mu - negparam$mu)/posparam$sigma
  B <- negparam$sigma/posparam$sigma
  Zx <- c(-Inf, seq(-3, 3, 0.01), Inf)
  Cutoff <- negparam$mu - negparam$sigma * Zx

  FPR <- pnorm(Zx)
  TPR <- pnorm(A + B * Zx)

  AUC <- pnorm(A/sqrt(1 + B^2))

  returnval <- list(method = "binormal", pos_count = nY,

                    neg_count = nYbar, pos_D = DY, neg_D = DYbar, AUC = AUC,

                    param = list(posparam = posparam, negparam = negparam),

                    TPR = TPR, FPR = FPR)
  return(returnval)
}

#' ML Estimates Function
#' Computes the mean and standard deviation of the input vector
#' @param x: a vector of numbers
#' @return a list containing the mean and standard deviation.
MLestimates <- function (x) {
  mu <- mean(x)
  sigma <- sqrt(mean((x - mu)^2))
  list(mu = mu, sigma = sigma)
}


# Example:
score=c(rnorm(100),rnorm(100,.6,1.2))
status=c(rep(0,100),rep(1,100))
res=bnROC(score,status)
res
