#!/usr/bin/env Rscript
# Faceted RATE Manhattan for a PPOM run: one facet per cutpoint (MIC breakpoint),
# RATE on the y-axis, genome coordinate on the x-axis, with gene names overlaid on
# each facet. This is the faceted RATE plot gwas_workflow writes
# (manhattan_all_cutpoints_faceted_RATE.png) plus the per-cutpoint gene labels it
# puts on manhattan_RATE_cutpoint{c}.png. Labels follow the same rule as the
# manhattan_median_effects_cutpoint1.png plots: the top 20 unique genes by |median|
# for that cutpoint, applied independently to each facet.
#
# Usage:
#   Rscript replot_rate_faceted_cutpoints.R --run-dir <inference output dir> \
#     --genes-of-interest <drug>_genesofinterest.txt \
#     [--annotations <snpEff fields .txt>] [--output-dir <dir>] [--ncol 2]
#
# Default --annotations is the SPN fields file; default --output-dir is
#   /nfs/research/jlees/jacqueline/thesis_results/paper_figures/<basename(run-dir)>/manhattan_short
#
# Output is manhattan_all_cutpoints_faceted_RATE_labeled.png.

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(ggrepel)
})

DEFAULT_PAPER_FIGURES_DIR <- "/nfs/research/jlees/jacqueline/thesis_results/paper_figures"
DEFAULT_ANNOTATIONS <- "/nfs/research/jlees/jacqueline/gwas_data/spn_pneumo/genotype/fields_filtered_maf05_multiallelic.txt"

parse_args <- function(argv) {
  out <- list(run_dir = NULL, output_dir = NULL,
              annotations = DEFAULT_ANNOTATIONS, genes_of_interest = NULL,
              ncol = 2L)
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a == "--run-dir") {
      out$run_dir <- argv[i + 1]; i <- i + 2
    } else if (a == "--output-dir") {
      out$output_dir <- argv[i + 1]; i <- i + 2
    } else if (a == "--annotations") {
      out$annotations <- argv[i + 1]; i <- i + 2
    } else if (a == "--genes-of-interest") {
      out$genes_of_interest <- argv[i + 1]; i <- i + 2
    } else if (a == "--ncol") {
      out$ncol <- as.integer(argv[i + 1]); i <- i + 2
    } else {
      stop("Unknown argument: ", a)
    }
  }
  if (is.null(out$run_dir)) stop("Required arg: --run-dir <dir>")
  if (is.null(out$output_dir)) {
    out$output_dir <- file.path(DEFAULT_PAPER_FIGURES_DIR,
                                basename(normalizePath(out$run_dir, mustWork = TRUE)),
                                "manhattan_short")
  }
  out
}

# phandango .plot is tab-separated with a BP (genome coordinate) column
read_positions <- function(phandango) {
  if (!file.exists(phandango)) stop("Missing: ", phandango)
  as.numeric(read.table(phandango, header = TRUE, sep = "\t")$BP)
}

# RATE_values_*.txt has #-prefixed header lines then "snp_id RATE KLD"; col 2 is RATE
read_rate <- function(rate_file) {
  if (!file.exists(rate_file)) stop("Missing: ", rate_file)
  as.numeric(read.table(rate_file, comment.char = "#")[, 2])
}

# Gene-label helpers, ported verbatim from gwas_workflow manhattan_plots.R so the
# rendered labels match (italicised gene names; display-name remapping).
.italic_gene_expr <- function(x) paste0("italic('", gsub("'", "", x), "')")

.apply_gene_display_names <- function(gene_names, genes_of_interest) {
  if (is.null(genes_of_interest) || is.null(gene_names)) return(gene_names)
  display_map <- stats::setNames(genes_of_interest[[2]], genes_of_interest[[1]])
  hit <- !is.na(gene_names) & gene_names %in% names(display_map)
  gene_names[hit] <- unname(display_map[gene_names[hit]])
  gene_names
}

# Top-10-unique-genes-by-RATE greedy selection (gwas_workflow rule, but ranked on
# RATE rather than |median|). Returns the row indices of `df` to label. Excludes
# NA / "MODIFIER" genes.
top_label_idx <- function(gene, score, n = 10L) {
  valid <- !is.na(gene) & gene != "MODIFIER"
  score[!valid] <- -Inf
  ord <- order(score, decreasing = TRUE)
  seen <- character(0); label_idx <- integer(0)
  for (i in ord) {
    g <- gene[i]
    if (valid[i] && !(g %in% seen)) {
      seen <- c(seen, g); label_idx <- c(label_idx, i)
      if (length(seen) >= n) break
    }
  }
  label_idx
}

build_df <- function(run_dir) {
  cutpoint_files <- sort(Sys.glob(file.path(run_dir, "cppRATE_results",
                                            "RATE_values_cutpoint*_depruned.txt")))
  if (length(cutpoint_files) == 0L) {
    stop("No RATE_values_cutpoint*_depruned.txt found; this is not a PPOM run: ", run_dir)
  }
  n_cutpoints <- length(cutpoint_files)

  # cutpoint -> breakpoint MIC label and per-cutpoint medians (file order matches
  # RATE/phandango row order for each cutpoint)
  effects <- read.csv(file.path(run_dir, "fitted_model", "depruned_variant_effects.csv"))
  has_mic <- "cutpoint_MIC" %in% names(effects) && !all(is.na(effects$cutpoint_MIC))
  if (has_mic) {
    map <- unique(effects[, c("cutpoint", "cutpoint_MIC")])
    map <- map[order(map$cutpoint), ]
    labels <- as.character(map$cutpoint_MIC)
    legend_title <- expression(atop("breakpoint", (mu * "g" %.% "mL"^{-1})))
  } else {
    labels <- as.character(seq_len(n_cutpoints))
    legend_title <- "cutpoint"
  }

  parts <- lapply(seq_len(n_cutpoints), function(c) {
    pos  <- read_positions(file.path(run_dir, "cppRATE_results",
                                     sprintf("phandango_cutpoint%d.plot", c)))
    rate <- read_rate(cutpoint_files[c])
    med  <- effects$median[effects$cutpoint == c]
    if (length(pos) != length(rate) || length(pos) != length(med)) {
      stop("Cutpoint ", c, " length mismatch: pos=", length(pos),
           " rate=", length(rate), " median=", length(med))
    }
    data.frame(pos = pos, rate = rate, median = med,
               cutpoint = c, cutpoint_label = labels[c])
  })
  df <- do.call(rbind, parts)
  df$cutpoint_label <- factor(df$cutpoint_label, levels = labels)
  list(df = df, n = n_cutpoints, legend_title = legend_title)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  run_dir <- normalizePath(args$run_dir, mustWork = TRUE)
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  message("Run dir:    ", run_dir)
  message("Output dir: ", args$output_dir)

  built <- build_df(run_dir)
  df <- built$df

  # POS -> gene from snpEff annotations; unannotated -> "MODIFIER". read.delim
  # mangles the header so the gene column becomes ANN....GENE.
  ann <- read.delim(args$annotations, stringsAsFactors = FALSE)
  df$gene <- ann[["ANN....GENE"]][match(df$pos, ann$POS)]
  df$gene[is.na(df$gene)] <- "MODIFIER"

  # Display-name remapping from the genes-of-interest list (col1 = annotation name,
  # col2 = display label), read exactly as gwas_workflow does.
  goi <- if (!is.null(args$genes_of_interest)) {
    read.csv(args$genes_of_interest, header = FALSE,
             col.names = c("gene", "display_name"), stringsAsFactors = FALSE)
  } else {
    NULL
  }
  df$gene <- .apply_gene_display_names(df$gene, goi)

  # Per-facet top-10-by-RATE labels, placed at the RATE value
  label_df <- do.call(rbind, lapply(split(df, df$cutpoint), function(d) {
    d[top_label_idx(d$gene, d$rate, n = 10L), ]
  }))
  label_df$gene_expr <- .italic_gene_expr(label_df$gene)

  single_color <- viridis::viridis(6, option = "plasma")[2]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = pos, y = rate)) +
    ggplot2::geom_point(alpha = 0.4, colour = single_color) +
    ggplot2::facet_wrap(~ cutpoint_label, ncol = args$ncol,
                        labeller = ggplot2::labeller(cutpoint_label = function(x) x)) +
    ggrepel::geom_text_repel(
      data  = label_df,
      ggplot2::aes(x = pos, y = rate, label = gene_expr),
      parse = TRUE,
      size  = 4.4,
      arrow = grid::arrow(length = grid::unit(0.01, "npc"), type = "open"),
      colour = "black",
      inherit.aes = FALSE,
      max.overlaps = 30
    ) +
    ggplot2::xlab("genome coordinate (bp)") +
    ggplot2::ylab("relative centrality (RATE)") +
    ggplot2::theme_minimal(base_size = 18)

  height <- 4.5 * ceiling(built$n / args$ncol)
  outfile <- paste0(basename(run_dir), "_faceted_RATE_labeled.png")
  path <- file.path(args$output_dir, outfile)
  ggplot2::ggsave(path, plot = p, width = 22, height = height, dpi = 300, limitsize = FALSE)
  message("Wrote ", path)
}

main()
