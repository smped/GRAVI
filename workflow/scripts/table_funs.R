#' A set of functions for decorating tables produced by `reactable`
#' These define colours based on 'up', 'down' and unchanged status,
#' as well as the values in the logFC and logCPM columns
sprint_pval <- function(p) {
  fmt <- ifelse(p < 0.01, "%.2e", "%.3f")
  zero <- p == 0
  p <- sprintf(fmt, p)
  p[zero] <- "<2e-16"
  p
}
comma_col <- colFormat(separators = TRUE, digits = 0)
percent_col <- colFormat(percent = TRUE, digits = 2)

#' JavaScript functions for column filtering using >= or <=
js_greater <- htmlwidgets::JS(
  "function(rows, columnId, filterValue) {
    return rows.filter(function(row) {
      return row.values[columnId] >= filterValue
    })
  }"
)
js_less <- htmlwidgets::JS(
  "function(rows, columnId, filterValue) {
    return rows.filter(function(row) {
      return row.values[columnId] <= filterValue
    })
  }"
)

up_col <- function(x) {
  if (is.na(x) | is.nan(x)) return("#ffffff")
  rgb(
    colorRamp(c("#ffffff", colours$direction[["increased"]]))(x), maxColorValue = 255
  )
}
down_col <- function(x) {
  if (is.na(x) | is.nan(x)) return("#ffffff")
  rgb(
    colorRamp(c("#ffffff", colours$direction[["decreased"]]))(x), maxColorValue = 255
  )
}
unch_col <- function(x) {
  if (is.na(x) | is.nan(x)) return("#ffffff")
  rgb(
    colorRamp(c("#ffffff", colours$direction[["unchanged"]]))(x),
    maxColorValue = 255
  )
}
lfc_col <- function(x){
  if (is.na(x) | is.nan(x)) return("#ffffff")
  rgb(
    colorRamp(c(colours$direction[["decreased"]], "#ffffff", colours$direction[["increased"]]))(x),
    maxColorValue = 255
  )
}
expr_col <- function(x){
  if (is.na(x) | is.nan(x)) return("#ffffff")
  rgb(colorRamp(hcl.colors(9, "TealRose"))(x), maxColorValue = 255)
}

#' The following enable the addition of bars within cells and the use of tooltips
bar_style <- function(width = 1, fill = "#e6e6e6", height = "75%", align = c("left", "right"), color = NULL, fontSize = c()) {
  align <- match.arg(align)
  if (align == "left") {
    position <- paste0(width * 100, "%")
    image <- sprintf("linear-gradient(90deg, %1$s %2$s, transparent %2$s)", fill, position)
  } else {
    position <- paste0(100 - width * 100, "%")
    image <- sprintf("linear-gradient(90deg, transparent %1$s, %2$s %1$s)", position, fill)
  }
  styles <- list(
    backgroundImage = image,
    backgroundSize = paste("100%", height),
    backgroundRepeat = "no-repeat",
    backgroundPosition = "center",
    color = color
  )
  if (!is.null(fontSize)) styles$fontSize = fontSize
  styles
}
with_tooltip <- function(value, width = 30) {
  tags$span(title = value, str_trunc(value, width))
}
