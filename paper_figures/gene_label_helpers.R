# Shared gene-labeling helpers for the short Manhattan builders
# (replot_rate_manhattan_short.R and replot_ppom_overlay_manhattans_short.R).
#
# select_label_rows() decides which variant rows carry a gene-name label:
#   "gene_list" - label the single top-scoring variant in each named gene
#   "top_n"     - label the top-n unique genes by score (greedy), skipping MODIFIER
# repel_label_layer() is the ggrepel layer shared by the RATE and beta builders.
#
# Gene display-name remapping and annotation lookup live in
# replot_rate_faceted_cutpoints.R (annotate_genes / .apply_gene_display_names); the
# helpers here operate on an already-annotated frame with a `gene` column.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

# "folA" -> "italic('folA')" so plotmath renders the gene name in italics
.gl_italic_gene_expr <- function(x) paste0("italic('", gsub("'", "", x), "')")

# Remap annotation gene names to display names from a 2-col genes-of-interest list
# (col1 = annotation name, col2 = display name). Ported from
# replot_rate_faceted_cutpoints.R so all gene-labeling helpers live in one place.
.apply_gene_display_names <- function(gene_names, genes_of_interest) {
  if (is.null(genes_of_interest) || is.null(gene_names)) return(gene_names)
  # trim whitespace so "dyr, folA (dhfR)" style lists map to clean display names
  display_map <- stats::setNames(trimws(genes_of_interest[[2]]),
                                 trimws(genes_of_interest[[1]]))
  hit <- !is.na(gene_names) & gene_names %in% names(display_map)
  gene_names[hit] <- unname(display_map[gene_names[hit]])
  gene_names
}

# POS -> gene from a snpEff annotations TSV (read.delim mangles the header so the
# gene column becomes ANN....GENE). Unannotated positions become "MODIFIER", then
# display-name remapping is applied. `df` must have a `pos` column; returns it with
# an added `gene` column.
annotate_genes <- function(df, annotations, genes_of_interest = NULL) {
  ann <- read.delim(annotations, stringsAsFactors = FALSE)
  df$gene <- ann[["ANN....GENE"]][match(df$pos, ann$POS)]
  df$gene[is.na(df$gene)] <- "MODIFIER"
  goi <- if (!is.null(genes_of_interest)) {
    read.csv(genes_of_interest, header = FALSE,
             col.names = c("gene", "display_name"), stringsAsFactors = FALSE)
  } else {
    NULL
  }
  df$gene <- .apply_gene_display_names(df$gene, goi)
  df
}

# df          annotated frame (must have `gene_col` and `score_col`)
# mode        "gene_list" (restrict to `genes`, top variant each) or
#             "top_n" (top `n` unique genes by score)
# returns the rows to label, with an added `gene_expr` plotmath column.
select_label_rows <- function(df, gene_col = "gene", score_col,
                              mode = c("top_n", "gene_list"),
                              genes = NULL, n = 10L) {
  mode  <- match.arg(mode)
  gene  <- df[[gene_col]]
  score <- df[[score_col]]

  if (mode == "gene_list") {
    if (is.null(genes)) stop("mode='gene_list' requires `genes`")
    keep <- which(gene %in% genes & !is.na(score))
    if (length(keep) == 0L) return(df[0, , drop = FALSE])
    # one row (the max-scoring variant) per requested gene that is present
    idx <- vapply(split(keep, gene[keep]),
                  function(ix) ix[which.max(score[ix])], integer(1))
    out <- df[sort(unname(idx)), , drop = FALSE]
  } else {
    valid <- !is.na(gene) & gene != "MODIFIER" & !is.na(score)
    s <- score; s[!valid] <- -Inf
    ord <- order(s, decreasing = TRUE)
    seen <- character(0); pick <- integer(0)
    for (i in ord) {
      g <- gene[i]
      if (valid[i] && !(g %in% seen)) {
        seen <- c(seen, g); pick <- c(pick, i)
        if (length(seen) >= n) break
      }
    }
    out <- df[pick, , drop = FALSE]
  }

  out$gene_expr <- .gl_italic_gene_expr(out[[gene_col]])
  out
}

# Shared label layer. `label_df` must carry `xcol`, `ycol` and `gene_expr`. Labels
# sit directly on top of their point (the region's peak SNP): repel only along y
# so the label stays centered over the peak's x, nudged up, with no connector
# segment. Dark-grey text.
repel_label_layer <- function(label_df, size = 5, xcol = "pos", ycol = "value",
                              colour = "grey25", nudge_y = 0) {
  ggrepel::geom_text_repel(
    data        = label_df,
    mapping     = ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]],
                               label = .data[["gene_expr"]]),
    parse              = TRUE,
    size               = size,
    colour             = colour,
    inherit.aes        = FALSE,
    direction          = "y",
    nudge_y            = nudge_y,
    vjust              = 0,
    min.segment.length = Inf,   # no connector arrows/segments
    box.padding        = 0.15,
    point.padding      = 0.2,
    max.overlaps       = 30,
    seed               = 1
  )
}
