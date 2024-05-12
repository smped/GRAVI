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
#' Print a range when aggregating a set of rows
range_js <- htmlwidgets::JS(
  "function(values) {
      var min_val = Math.round(100 * Math.min(...values)) / 100
      var max_val = Math.round(100 * Math.max(...values)) / 100
      if (min_val == max_val) {
        var ret_val = min_val.toString()
      } else {
        var ret_val = '[' + min_val.toString() + ', ' + max_val.toString() + ']'
      }
      return ret_val
    }"
)
#' Aggregate across multiple GRanges (concatenated as characters)
granges_js <- htmlwidgets::JS(
  "function(values){
      var chrom = [];
      var rng = [];
      var start = [];
      var end = [];
      for (i = 0; i < values.length; i++) {
        chrom[i] = values[i].split(':')[0];
        rng[i] = values[i].split(':')[1];
        start[i] = rng[i].split('-')[0];
        end[i] = rng[i].split('-')[1];
      }

      min_start = Math.min(...start);
      max_end = Math.max(...end);
      var ret_val = [...new Set(chrom)] + ':' + min_start.toString() + '-' + max_end.toString();
      return ret_val
    }"
)
drop_filter <- function(values, name, id = NULL, sort = TRUE) {
  if (sort) values <- sort(values)
  tags$select(
    # Set to undefined to clear the filter
    onchange = sprintf(
      "Reactable.setFilter('%s', '%s', event.target.value || undefined)",
      id, name
    ),
    # "All" has an empty value to clear the filter, and is the default option
    tags$option(value = "", "All"),
    lapply(unique(values), tags$option),
    "aria-label" = sprintf("Filter %s", name),
    style = "width: 100%; height: 28px;"
  )
}

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

#' Setup the columns for enrichment tables
enrich_cols <- list(
  gs_name = colDef(
    "GeneSet", minWidth = 180,
    cell = function(value) htmltools::tags$a(
      href = gs_url[[value]],
      target = "_blank",
      str_replace_all(value, "_", " ")
    ),
    html = TRUE
  ),
  genome_fraction = colDef(show = FALSE),
  observed_region_hits = colDef(show = FALSE),
  mean_tss_dist = colDef(
    "Mean TSS Distance (kb)", cell = \(value) round(value / 1e3, 2),
    filterMethod = js_greater
  ),
  gene_set_size = colDef(
    "Gene Set Size", maxWidth = 100, filterMethod = js_greater
  ),
  observed_gene_hits = colDef(
    "Gene Hits", maxWidth = 100, filterMethod = js_greater
  ),
  fold_enrichment = colDef(
    name = "Fold Enrichment", format = colFormat(digits = 3),
    minWidth = 110, filterMethod = js_greater
  ),
  genes_with_hits = colDef(
    name = "Genes With Associated Peaks", minWidth = 200,
    cell = \(value) with_tooltip(value, width = 60)
  ),
  gene_id = colDef(show = FALSE),
  p = colDef(show = FALSE),
  adj_p = colDef(
    name = glue("P<sub>adj</sub>"),
    html = TRUE, cell = \(value) sprint_pval(value),
    maxWidth = 110, filterMethod = js_less
  )
)

#' columns commonly in motif enrichment results
motif_cols <- list(
  cluster = colDef("Cluster", maxWidth = 70),
  altname = colDef(
    "Motif", minWidth = 180, aggregate = "unique",
    style = list(borderLeft = "1px solid rgba(0, 0, 0, 0.1)")
  ),
  name = colDef(
    "Name", minWidth = 120, aggregate = "unique",
    style = list(borderRight = "1px solid rgba(0, 0, 0, 0.1)")
  ),
  matches = colDef(
    "Total", maxWidth = 60, format = comma_col, aggregate = "max",
    filterMethod = js_greater
  ),
  expected = colDef(
    "Expected", maxWidth = 80, format = comma_col, aggregate = "max",
    filterMethod = js_greater
  ),
  enrichment = colDef(
    "Enrichment", maxWidth = 90, format = colFormat(digits = 3),
    aggregate = "max", filterMethod = js_greater,
    style = list(borderRight = "1px solid rgba(0, 0, 0, 0.1)")
  ),
  p = colDef(show = FALSE),
  adj_p = colDef(
    glue("p<sub>{motif_params$adj}</sub>"), html = TRUE, maxWidth = 70,
    cell = \(value) sprint_pval(value), aggregate = "min",
    format = colFormat(digits = 3), filterMethod = js_less
  ),
  PWM = colDef(
    name = "IC Matrix", minWidth = 180, filterable = FALSE,
    cell = function(value) tags$img(src = value, height = '80px')
  ),
  start = colDef(
    "Start", aggregate = "min", maxWidth = 50, filterMethod = js_greater
  ),
  end = colDef(
    "End", aggregate = "max", maxWidth = 50, filterMethod = js_greater
  ),
  centre = colDef(
    "Centre", aggregate = "median", maxWidth = 60, show = FALSE
  ),
  width = colDef(
    "Width", aggregate = "median", maxWidth = 60, filterMethod = js_less
  ),
  total_matches = colDef(
    "Total", aggregate = "max", maxWidth = 60, format = comma_col,
    filterMethod = js_greater
  ),
  matches_in_region = colDef(
    "Region", aggregate = "max", maxWidth = 60, filterMethod = js_greater
  ),
  prop_total = colDef(
    "% In Region", aggregate = "max", format = percent_col, maxWidth = 60,
    filterMethod = js_greater
  ),
  odds_ratio = colDef("Odds Ratio", show = FALSE)
)
