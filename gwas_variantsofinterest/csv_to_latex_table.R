#!/usr/bin/env Rscript
# Turn a GOI-by-position CSV (as written by build_pos_table.R, columns like run, gene,
# POS, best_cp, RATE, beta, ref_residue, variants, n_alleles_nt, n_alleles_aa, effects,
# all_HGVS_P) into a paper-ready LaTeX table, written to a .txt file.
#
# One booktabs table per input CSV. Numeric columns get scientific notation, gene names
# are italicised, text is LaTeX-escaped. Columns are picked automatically from whatever
# the file contains (minus a drop list); unknown columns fall through as escaped text.
#
# Usage:
#   mamba activate gwas_pipeline
#   Rscript csv_to_latex_table.R <in1.csv> [in2.csv ...] [out.txt] [--caption "..."] [--label prefix]
#
# Defaults: out.txt = <first input>_latextable.txt in the same directory.
#   With several inputs, each gets its own table (caption/label derived from its filename)
#   concatenated into the single out.txt.

args <- commandArgs(trailingOnly = TRUE)

# ---- argument parsing -------------------------------------------------------
caption_arg <- NULL; label_arg <- NULL; group_arg <- "gene,binning"; longtable <- FALSE; subblank_arg <- ""
i <- which(args == "--caption"); if (length(i)) { caption_arg <- args[i + 1]; args <- args[-c(i, i + 1)] }
i <- which(args == "--label");   if (length(i)) { label_arg   <- args[i + 1]; args <- args[-c(i, i + 1)] }
i <- which(args == "--group");   if (length(i)) { group_arg   <- args[i + 1]; args <- args[-c(i, i + 1)] }
i <- which(args == "--subblank"); if (length(i)) { subblank_arg <- args[i + 1]; args <- args[-c(i, i + 1)] }  # cols blanked on repeat within a block
i <- which(args == "--longtable"); if (length(i)) { longtable <- TRUE; args <- args[-i] }  # page-breaking landscape longtable
group_cols    <- trimws(strsplit(group_arg, ",")[[1]])   # one or more columns form the block key
subblank_cols <- if (nzchar(subblank_arg)) trimws(strsplit(subblank_arg, ",")[[1]]) else character(0)

inputs <- args[grepl("\\.csv$", args)]
outp   <- args[grepl("\\.txt$", args)]
if (length(inputs) == 0L) stop("give at least one input .csv")
out_file <- if (length(outp)) outp[1] else sub("\\.csv$", "_latextable.txt", inputs[1])

# ---- column spec ------------------------------------------------------------
# header label + column type per known column; drop = never emitted.
HEADERS <- c(gene = "Gene", gene_eggnog = "eggNOG gene", desc_eggnog = "eggNOG description",
             binning = "Binning", model = "Model", POS = "Position",
             best_cp = "cp", RATE = "RATE", beta = "$|\\beta|$", ref_residue = "Residue",
             variants = "Variants", n_alleles_nt = "$n_{\\mathrm{nt}}$",
             n_alleles_aa = "$n_{\\mathrm{aa}}$", effects = "Effects")
NUMCOLS  <- c("RATE", "beta")                                  # scientific-notation cols
INTCOLS  <- c("POS", "best_cp", "n_alleles_nt", "n_alleles_aa") # right-aligned integers
WRAP     <- c(effects = "3cm", desc_eggnog = "5cm")            # wrapped fixed-width columns
GENE_ATTACH <- c("gene_eggnog", "desc_eggnog")   # gene-level; placed by gene, shown once per gene
DROP     <- c("all_HGVS_P", "eggnog_pident", "best_cp")        # too long / uninformative

# ---- formatters -------------------------------------------------------------
tex_escape <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  for (ch in c("_", "%", "#", "&", "{", "}")) x <- gsub(ch, paste0("\\", ch), x, fixed = TRUE)
  x <- gsub("~", "\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("^", "\\textasciicircum{}", x, fixed = TRUE)
  x
}

fmt_sci <- function(x) {                       # 3 sig figs; sci notation outside [1e-3,1e4)
  vapply(x, function(v) {
    if (is.na(v)) return("")
    v <- as.numeric(v); if (v == 0) return("0")
    a <- abs(v)
    if (a >= 1e-3 && a < 1e4) return(format(signif(v, 3), scientific = FALSE, trim = TRUE))
    e <- floor(log10(a)); m <- v / 10^e
    sprintf("$%.2f\\mathrm{e}{%d}$", m, e)   # compact e-notation, e.g. 4.67e-5
  }, character(1))
}

format_col <- function(name, values) {
  if (name == "beta") values <- abs(as.numeric(values))     # tables report |beta|
  if (name %in% NUMCOLS) return(fmt_sci(values))
  if (name %in% INTCOLS) return(format(as.integer(values), trim = TRUE))
  if (name == "binning") {                              # keep >= / % as math, no escaping
    v <- gsub("%", "\\%", as.character(values), fixed = TRUE)
    return(gsub(">=", "$\\geq$", v, fixed = TRUE))
  }
  if (name == "effects") {                              # drop "_variant", "_"->" ", comma-join
    return(vapply(strsplit(as.character(values), ";"), function(v) {
      v <- gsub("_", " ", sub("_variant$", "", trimws(v)))
      paste(tex_escape(v), collapse = ", ")
    }, character(1)))
  }
  txt <- tex_escape(gsub(";", ", ", as.character(values)))   # ";"-lists -> ", "
  if (name == "gene") txt <- paste0("\\textit{", txt, "}")
  txt
}

# ---- table builder ----------------------------------------------------------
# Redundancy is reduced by grouping on one or more `group` columns (the block key). Rows
# are ordered so equal block keys are contiguous (first-appearance order of each column,
# so a pre-sorted CSV keeps its order), every group cell is blanked after the first row
# of its block, and a \midrule separates blocks -- the style of the reference table.
build_table <- function(csv, caption, label, group = "gene", longtable = FALSE, subblank = character(0)) {
  d <- read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE)
  cols <- setdiff(names(d), DROP)
  cols <- cols[cols %in% names(HEADERS)]                 # keep only known, in file order
  group <- group[group %in% cols]                        # ignore group cols not present
  if (length(group)) cols <- c(group, setdiff(cols, group))  # group columns leftmost
  # gene-level columns (eggNOG) sit right after `gene`, shown once per gene
  attach_cols <- intersect(GENE_ATTACH, cols)
  if (length(attach_cols) && "gene" %in% cols)
    cols <- c("gene", attach_cols, setdiff(cols, c("gene", attach_cols)))

  # order rows by the block key, respecting each group column's first-appearance order
  if (length(group)) {
    keys <- lapply(group, function(c) factor(d[[c]], levels = unique(d[[c]])))
    d <- d[do.call(order, keys), , drop = FALSE]
  }

  header <- paste(HEADERS[cols], collapse = " & ")
  align  <- vapply(cols, function(c)
    if (c %in% names(WRAP)) sprintf(">{\\raggedright\\arraybackslash}p{%s}", WRAP[[c]])
    else if (c %in% c(NUMCOLS, INTCOLS)) "r" else "l", character(1))

  mat <- do.call(cbind, lapply(cols, function(c) format_col(c, d[[c]])))
  colnames(mat) <- cols
  key  <- if (length(group)) apply(d[group], 1, paste, collapse = "\r") else seq_len(nrow(d))
  gval <- if ("gene" %in% names(d)) as.character(d[["gene"]]) else rep(NA_character_, nrow(d))

  sb <- intersect(subblank, cols)                        # blanked on repeat within a block
  body <- character(0); prev <- NULL; prevg <- NULL; prevsb <- NULL
  for (r in seq_len(nrow(mat))) {
    row <- mat[r, ]
    if (!is.null(prev) && key[r] != prev) body <- c(body, "\\midrule")   # block boundary
    else if (!is.null(prev)) {
      for (g in group) row[[g]] <- ""                                    # blank within block
      for (s in sb) if (!is.null(prevsb) && mat[r, s] == prevsb[[s]]) row[[s]] <- ""
    }
    if (!is.null(prevg) && !is.na(gval[r]) && gval[r] == prevg)          # eggNOG once per gene
      for (a in attach_cols) row[[a]] <- ""
    body <- c(body, paste0(paste(row, collapse = " & "), " \\\\"))
    prev <- key[r]; prevg <- gval[r]; prevsb <- as.list(mat[r, ])
  }

  if (longtable) {
    # A page-breaking longtable: the header repeats on every page and a "continued" line
    # marks the break. `\rowcolor` is attached inline to the header row (a standalone
    # \rowcolor line is typeset as its own header row in longtable). Needs longtable.
    n  <- length(cols)
    hd <- c("\\toprule", paste0("\\rowcolor{gray!15} ", header, " \\\\"), "\\midrule")
    return(c("\\begin{landscape}",
      "\\footnotesize", "\\setlength{\\tabcolsep}{4pt}",
      "\\renewcommand{\\arraystretch}{1.08}",
      sprintf("\\begin{longtable}{@{} %s @{}}", paste(align, collapse = " ")),
      sprintf("\\caption{%s}\\label{tab:%s}\\\\", caption, label),
      hd, "\\endfirsthead",
      sprintf("\\multicolumn{%d}{c}{\\footnotesize\\itshape \\tablename\\ \\thetable{} -- continued from previous page}\\\\", n),
      hd, "\\endhead",
      "\\midrule",
      sprintf("\\multicolumn{%d}{r}{\\footnotesize\\itshape Continued on next page}\\\\", n),
      "\\endfoot",
      "\\bottomrule", "\\endlastfoot",
      body, "\\end{longtable}", "\\end{landscape}", ""))
  }

  c("\\begin{table}[htbp]", "\\centering",
    "\\footnotesize", "\\setlength{\\tabcolsep}{4pt}",
    "\\renewcommand{\\arraystretch}{1.08}",
    "\\rotatebox{90}{%",                                 # graphicx: rotate table 90 deg
    sprintf("\\begin{tabular}{@{} %s @{}}", paste(align, collapse = " ")),
    "\\toprule", paste0("\\rowcolor{gray!15} ", header, " \\\\"), "\\midrule",
    body, "\\bottomrule", "\\end{tabular}%", "}",
    sprintf("\\caption{%s}", caption), sprintf("\\label{tab:%s}", label),
    "\\end{table}", "")
}

nice_label <- function(csv) gsub("[^A-Za-z0-9]+", "_", tools::file_path_sans_ext(basename(csv)))
nice_cap   <- function(csv) {
  s <- gsub("_", " ", tools::file_path_sans_ext(basename(csv)))
  paste0(toupper(substring(s, 1, 1)), substring(s, 2), ".")
}

lines <- unlist(lapply(seq_along(inputs), function(k) {
  cap <- if (!is.null(caption_arg) && length(inputs) == 1L) caption_arg else nice_cap(inputs[k])
  lab <- if (!is.null(label_arg)) paste0(label_arg, if (length(inputs) > 1L) paste0("_", k) else "")
         else nice_label(inputs[k])
  build_table(inputs[k], cap, lab, group = group_cols, longtable = longtable, subblank = subblank_cols)
}))

writeLines(lines, out_file)
cat("Wrote", out_file, "(", length(inputs), "table(s) )\n")
