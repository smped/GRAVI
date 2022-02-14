#' @title Call Super-Enhancers
#'
#' @description Functions for identifying super-enhancers
#'
#' @details
#' These two functions are taken from
#' https://github.com/rakarnik/ROSE/blob/master/ROSE_callSuper.R
#' and used with only superficial modifications
#'
#' @param x a numeric vector
#' @param drawPlot logical. Draws plot if required
#' @param ... Passed to `plot()`
#'
calculate_cutoff <- function(x, drawPlot = TRUE, ...) {

  x <- sort(x)
  #set those regions with more control than ranking equal to zero
  x[x < 0] <- 0

  #This is the slope of the line we want to slide. This is the diagonal.
  slope <- (max(x) - min(x)) / length(x)
  #Find the x-axis point where a line passing through that point has
  # the minimum number of points below it. (ie. tangent)
  xPt <- floor(
    optimize(
      .numPts_below_line,
      lower = 1,
      upper = length(x),
      myVector = x,
      slope = slope
    )$minimum
  )
  #The y-value at this x point. This is our cutoff.
  y_cutoff <- x[xPt]

  if (drawPlot) {
    plot(seq_along(x), x, type = "l", ...)
    b <- y_cutoff - (slope * xPt)
    abline(v = xPt, h = y_cutoff, lty = 2, col = 8)
    points(xPt, y_cutoff, pch = 16, cex = 0.9, col = 2)
    abline(coef = c(b, slope), col = 2)
    text(
      x = 0.25 * length(x),
      y = 0.9 * max(x),
      labels = paste(
        "x=", xPt,
        "\ny=", signif(y_cutoff, 3),
        "\nFold over Median=", signif(y_cutoff / median(x), 3),
        "\nFold over Mean=", signif(y_cutoff / mean(x), 3), "x",
        "\nDetected=", length(x) - xPt + 1,
        sep = ""
      )
    )
    axis(
      1,
      sum(x == 0),
      sum(x == 0),
      col.axis = "pink", col = "pink"
    ) #Number of regions with zero signal
  }

  list(
    absolute = y_cutoff,
    overMedian = y_cutoff / median(x),
    overMean = y_cutoff / mean(x)
    )

}

.numPts_below_line <- function(myVector, slope, x) {
  yPt <- myVector[x]
  b <- yPt - (slope * x)
  xPts <- seq_along(myVector)

  sum(myVector <= (xPts * slope + b))
}