#!/usr/bin/env Rscript
# Compute prior-based ROPE significance for a regularized horseshoe GWAS fit
# and produce labeled Manhattan plots in the style of gwas_workflow's
# manhattan_exp_median plots.
#
# Outputs go to <run-dir>/rope_pd_bf/.

suppressPackageStartupMessages({
  library(optparse)
  library(posterior)
  library(ggplot2)
  library(ggrepel)
})

# ---- CLI --------------------------------------------------------------------

option_list <- list(
  make_option("--run-dir",          type = "character", default = NULL,
              help = "Parent of fitted_model/ (required)"),
  make_option("--snpeff",           type = "character", default = NULL,
              help = "snpEff TSV (POS, ANN....GENE, ANN....HGVS_P, ANN....HGVS_C, [ANN....EFFECT])"),
  make_option("--variant-index",    type = "character", default = NULL,
              help = "<dataset>_variant_index.csv (required)"),
  make_option("--rope-prob",        type = "double",    default = 0.89,
              help = "Central probability of the prior null defining the ROPE [%default]"),
  make_option("--ci-prob",          type = "double",    default = 0.89,
              help = "Posterior CI probability used for the CI-significance criterion [%default]"),
  make_option("--n-prior-samples",  type = "integer",   default = 200000L,
              help = "Total prior null samples to pool [%default]"),
  make_option("--genes-of-interest", type = "character", default = NULL,
              help = "Optional 1- or 2-col file of gene names to highlight")
)

opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt[["run-dir"]]),
          !is.null(opt[["snpeff"]]),
          !is.null(opt[["variant-index"]]))

run_dir       <- opt[["run-dir"]]
snpeff_path   <- opt[["snpeff"]]
variant_index_path <- opt[["variant-index"]]
rope_prob     <- opt[["rope-prob"]]
ci_prob       <- opt[["ci-prob"]]
n_prior       <- opt[["n-prior-samples"]]
goi_path      <- opt[["genes-of-interest"]]

out_dir <- file.path(run_dir, "rope_pd_bf")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load the most recent fitted RDS ----------------------------------------

fitted_dir <- file.path(run_dir, "fitted_model")
rds_files <- list.files(fitted_dir, pattern = "\\.RDS$", full.names = TRUE)
if (length(rds_files) == 0L) stop("No .RDS in ", fitted_dir)
rds_path <- rds_files[which.max(file.mtime(rds_files))]
message("[rope] loading ", rds_path)
fit <- readRDS(rds_path)

# ---- Choose tau prior scale based on model type -----------------------------
# Logistic / continuous models use tau_0 = 0.1; POM / PPOM use tau_0 = 1.
# We only use this for documentation in rope_bounds.txt -- the ROPE itself is
# built from the *posterior* draws of tau and c2, so it doesn't depend on the
# nominal prior scale.

detect_model_type <- function(draw_names) {
  if (any(grepl("^cutpoints\\[", draw_names))) {
    if (any(grepl("^beta_variant_std\\[\\d+,\\d+\\]$", draw_names))) "PPOM" else "POM"
  } else if (any(grepl("^sigma$", draw_names)) && any(grepl("^beta_variant_std\\[", draw_names))) {
    "continuous"
  } else {
    "logistic"
  }
}

# Pull only the parameter columns we need - the full draws df can be huge.
need_vars <- c("tau", "c2", "beta_variant_std", "beta_variant")
draws_all <- tryCatch(
  posterior::as_draws_df(fit$draws(variables = need_vars)),
  error = function(e) posterior::as_draws_df(fit$draws())
)
draws_all <- as.data.frame(draws_all)

model_type <- detect_model_type(names(draws_all))
message("[rope] model type: ", model_type)

# ---- Variant index ----------------------------------------------------------

variant_index <- read.csv(variant_index_path, stringsAsFactors = FALSE)
stopifnot(all(c("variant_name", "position") %in% names(variant_index)))
V <- nrow(variant_index)
message("[rope] V = ", V)

# ---- Detect shape (vector vs PPOM matrix) and gather index info -------------

is_ppom <- model_type == "PPOM"
if (is_ppom) {
  # beta_variant_std columns look like "beta_variant_std[v,k]"
  bvs_cols <- grep("^beta_variant_std\\[\\d+,\\d+\\]$", names(draws_all), value = TRUE)
  idx <- regmatches(bvs_cols, regexec("^beta_variant_std\\[(\\d+),(\\d+)\\]$", bvs_cols))
  vk  <- do.call(rbind, lapply(idx, function(m) as.integer(m[-1])))
  K_m1 <- max(vk[, 2])
  message("[rope] PPOM: K-1 = ", K_m1)
} else {
  bvs_cols <- grep("^beta_variant_std\\[\\d+\\]$", names(draws_all), value = TRUE)
  K_m1 <- 1L
}

# ---- Build the prior null distribution from posterior tau, c2 ---------------
# beta_null_std = z * tau * sqrt( (c2 * lambda^2) / (c2 + tau^2 * lambda^2) )
# with z ~ N(0,1) and lambda ~ half-Cauchy(0, 2) (the model's prior on lambda_variant).

tau_post <- draws_all[["tau"]]
c2_post  <- draws_all[["c2"]]
S <- length(tau_post)
M <- max(1L, ceiling(n_prior / S))

message("[rope] sampling prior null: ", S, " posterior draws x ", M, " prior samples = ", S * M)
lambda_null <- abs(rcauchy(S * M, 0, 2))
z_null      <- rnorm(S * M)
tau_rep <- rep(tau_post, each = M)
c2_rep  <- rep(c2_post,  each = M)
lambda_tilde_null <- sqrt((c2_rep * lambda_null^2) /
                          (c2_rep + tau_rep^2 * lambda_null^2))
beta_null_std <- z_null * tau_rep * lambda_tilde_null

q_lo <- (1 - rope_prob) / 2
q_hi <- 1 - q_lo
rope <- unname(quantile(beta_null_std, probs = c(q_lo, q_hi)))
message(sprintf("[rope] standardized ROPE (%.0f%%): [%.5g, %.5g]",
                100 * rope_prob, rope[1], rope[2]))

# Free memory before the per-variant summary loop.
rm(lambda_null, z_null, tau_rep, c2_rep, lambda_tilde_null, beta_null_std)
gc(verbose = FALSE)

# ---- snpEff annotation ------------------------------------------------------

ann <- read.delim(snpeff_path, stringsAsFactors = FALSE)
ann_idx <- match(variant_index$position, ann$POS)

ann_get <- function(col, default = NA_character_) {
  if (!col %in% names(ann)) return(rep(default, V))
  out <- ann[[col]][ann_idx]
  out
}

gene_names <- ann_get("ANN....GENE", default = NA_character_)
gene_names[is.na(gene_names)] <- "MODIFIER"
hgvs_p     <- ann_get("ANN....HGVS_P")
hgvs_c     <- ann_get("ANN....HGVS_C")
effect_imp <- ann_get("ANN....EFFECT")

format_variant_labels <- function(hgvs_p, hgvs_c, fallback_id) {
  is_blank <- function(v) is.na(v) | !nzchar(v) | v == "."
  lbl <- character(length(fallback_id))
  has_p <- !is_blank(hgvs_p)
  has_c <- !is_blank(hgvs_c)
  lbl[has_p]           <- sub("^p\\.", "", hgvs_p[has_p])
  lbl[!has_p & has_c]  <- hgvs_c[!has_p & has_c]
  lbl[!has_p & !has_c] <- as.character(fallback_id[!has_p & !has_c])
  lbl
}

variant_labels <- format_variant_labels(hgvs_p, hgvs_c, variant_index$variant_name)

# Optional genes-of-interest display-name remap (col1 = ann name, col2 = display).
if (!is.null(goi_path) && file.exists(goi_path)) {
  goi <- tryCatch(read.delim(goi_path, header = FALSE, stringsAsFactors = FALSE),
                  error = function(e) NULL)
  if (!is.null(goi) && ncol(goi) >= 2) {
    map <- setNames(goi[[2]], goi[[1]])
    hit <- !is.na(gene_names) & gene_names %in% names(map)
    gene_names[hit] <- unname(map[gene_names[hit]])
  }
}

# ---- Per-variant posterior summary helper -----------------------------------

ci_lo_p <- (1 - ci_prob) / 2
ci_hi_p <- 1 - ci_lo_p

summarize_variant <- function(beta_std_mat, beta_mat) {
  # beta_*_mat: S x V matrices
  med_std <- matrixStats_colMedians(beta_std_mat)
  ci_std  <- matrixStats_colQuantiles(beta_std_mat, c(ci_lo_p, ci_hi_p))
  med_b   <- matrixStats_colMedians(beta_mat)
  ci_b    <- matrixStats_colQuantiles(beta_mat, c(ci_lo_p, ci_hi_p))
  list(
    median_std = med_std,
    ci_lo_std  = ci_std[, 1],
    ci_hi_std  = ci_std[, 2],
    median_beta = med_b,
    ci_low      = ci_b[, 1],
    ci_high     = ci_b[, 2]
  )
}

# Lightweight column median / quantile to avoid needing matrixStats.
matrixStats_colMedians <- function(m) apply(m, 2L, median)
matrixStats_colQuantiles <- function(m, p) t(apply(m, 2L, function(x) quantile(x, probs = p)))

# ---- Build per-variant rows -------------------------------------------------

build_row_df <- function(summary, cutpoint = NULL) {
  df <- data.frame(
    variant_name  = variant_index$variant_name,
    position      = variant_index$position,
    gene          = gene_names,
    hgvs_c        = hgvs_c,
    hgvs_p        = hgvs_p,
    effect_impact = effect_imp,
    median_beta   = summary$median_beta,
    ci_low        = summary$ci_low,
    ci_high       = summary$ci_high,
    stringsAsFactors = FALSE
  )
  df$signif_median <- summary$median_std < rope[1] | summary$median_std > rope[2]
  df$signif_ci     <- summary$ci_hi_std  < rope[1] | summary$ci_lo_std  > rope[2]
  if (!is.null(cutpoint)) df$cutpoint <- cutpoint
  df
}

# ---- Manhattan plot helper --------------------------------------------------
# Style mirrored from gwas_workflow/R/manhattan_plots.R (manhattan_exp_median),
# with two key differences:
#   * label every significant variant (no top-20 cap)
#   * subtitle states which criterion drove the plot
italic_gene_expr <- function(x) paste0("italic('", gsub("'", "", x), "')")

write_manhattan <- function(df, signif_col, criterion_label, out_path,
                            rope_bounds) {
  pdf <- data.frame(
    pos    = as.numeric(df$position),
    median = as.numeric(df$median_beta),
    signif = ifelse(df[[signif_col]], "significant", "not significant"),
    gene   = df$gene,
    hgvs   = df$hgvs_p,
    stringsAsFactors = FALSE
  )
  pdf$signif <- factor(pdf$signif, levels = c("not significant", "significant"))
  sig_colors <- c("not significant" = "grey60", "significant" = "red")

  label_df <- pdf[df[[signif_col]], , drop = FALSE]
  if (nrow(label_df) > 0L) {
    gene_for_lbl <- label_df$gene
    fallback     <- format_variant_labels(label_df$hgvs, df$hgvs_c[df[[signif_col]]],
                                          df$variant_name[df[[signif_col]]])
    use_gene <- !is.na(gene_for_lbl) & gene_for_lbl != "MODIFIER"
    lbl_text <- ifelse(use_gene, gene_for_lbl, fallback)
    label_df$gene_expr <- ifelse(use_gene,
                                 italic_gene_expr(lbl_text),
                                 paste0("'", gsub("'", "", lbl_text), "'"))
  }

  subtitle <- sprintf("%s outside ROPE [%.4g, %.4g]", criterion_label,
                      rope_bounds[1], rope_bounds[2])

  p <- ggplot(pdf, aes(x = pos, y = exp(median), colour = signif)) +
    geom_point(alpha = 0.4) +
    scale_colour_manual(values = sig_colors) +
    xlab("genome coordinate (bp)") +
    ylab(expression("e"^{tilde(beta)})) +
    labs(colour = "ROPE", subtitle = subtitle) +
    theme_minimal(base_size = 14)

  if (nrow(label_df) > 0L) {
    p <- p + geom_text_repel(
      data = label_df,
      aes(x = pos, y = exp(median), label = gene_expr),
      parse = TRUE,
      size = 3,
      arrow = grid::arrow(length = grid::unit(0.01, "npc"), type = "open"),
      colour = "black",
      inherit.aes = FALSE,
      max.overlaps = Inf
    )
  }

  ggsave(out_path, plot = p, width = 16, height = 6, dpi = 300)
}

# ---- Run summary + write outputs --------------------------------------------

extract_beta_matrix <- function(prefix, k = NULL) {
  if (is.null(k)) {
    cols <- paste0(prefix, "[", seq_len(V), "]")
  } else {
    cols <- paste0(prefix, "[", seq_len(V), ",", k, "]")
  }
  missing <- !cols %in% names(draws_all)
  if (any(missing)) stop("Missing columns: ", paste(head(cols[missing], 5), collapse=", "))
  as.matrix(draws_all[, cols, drop = FALSE])
}

write_rope_bounds <- function() {
  lines <- c(
    sprintf("model_type:           %s", model_type),
    sprintf("rope_prob:            %.4f", rope_prob),
    sprintf("ci_prob:              %.4f", ci_prob),
    sprintf("rope_lo (standardized): %.6g", rope[1]),
    sprintf("rope_hi (standardized): %.6g", rope[2]),
    sprintf("tau posterior median: %.6g", median(tau_post)),
    sprintf("c2  posterior median: %.6g", median(c2_post)),
    sprintf("n_posterior_draws:    %d", S),
    sprintf("n_prior_samples:      %d", S * M)
  )
  writeLines(lines, file.path(out_dir, "rope_bounds.txt"))
}

write_rope_bounds()

if (!is_ppom) {
  beta_std_mat <- extract_beta_matrix("beta_variant_std")
  beta_mat     <- extract_beta_matrix("beta_variant")
  summary <- summarize_variant(beta_std_mat, beta_mat)
  rm(beta_std_mat, beta_mat); gc(verbose = FALSE)

  df <- build_row_df(summary)
  write.csv(df, file.path(out_dir, "all_variants.csv"), row.names = FALSE)
  write.csv(df[df$signif_median, , drop = FALSE],
            file.path(out_dir, "significant_median.csv"), row.names = FALSE)
  write.csv(df[df$signif_ci, , drop = FALSE],
            file.path(out_dir, "significant_ci.csv"), row.names = FALSE)
  message(sprintf("[rope] median-significant: %d / %d", sum(df$signif_median), V))
  message(sprintf("[rope] CI-significant:     %d / %d", sum(df$signif_ci), V))

  write_manhattan(df, "signif_median", "median effect",
                  file.path(out_dir, "manhattan_median.png"), rope)
  write_manhattan(df, "signif_ci",     sprintf("%.0f%% CI", 100 * ci_prob),
                  file.path(out_dir, "manhattan_ci.png"), rope)

} else {
  # PPOM: per-cutpoint plus union.
  per_cutpoint_dfs <- vector("list", K_m1)
  for (k in seq_len(K_m1)) {
    beta_std_mat <- extract_beta_matrix("beta_variant_std", k)
    beta_mat     <- extract_beta_matrix("beta_variant",     k)
    summary <- summarize_variant(beta_std_mat, beta_mat)
    rm(beta_std_mat, beta_mat); gc(verbose = FALSE)

    dfk <- build_row_df(summary, cutpoint = k)
    per_cutpoint_dfs[[k]] <- dfk

    write.csv(dfk[dfk$signif_median, , drop = FALSE],
              file.path(out_dir, sprintf("significant_median_cutpoint%d.csv", k)),
              row.names = FALSE)
    write.csv(dfk[dfk$signif_ci, , drop = FALSE],
              file.path(out_dir, sprintf("significant_ci_cutpoint%d.csv", k)),
              row.names = FALSE)
    message(sprintf("[rope] cutpoint %d: median %d / CI %d significant",
                    k, sum(dfk$signif_median), sum(dfk$signif_ci)))

    write_manhattan(dfk, "signif_median",
                    sprintf("median effect (cutpoint %d)", k),
                    file.path(out_dir, sprintf("manhattan_median_cutpoint%d.png", k)),
                    rope)
    write_manhattan(dfk, "signif_ci",
                    sprintf("%.0f%% CI (cutpoint %d)", 100 * ci_prob, k),
                    file.path(out_dir, sprintf("manhattan_ci_cutpoint%d.png", k)),
                    rope)
  }

  # Union: variant significant at >=1 cutpoint, with count column and per-k flags.
  build_union <- function(crit) {
    sig_mat <- do.call(cbind, lapply(per_cutpoint_dfs, function(d) d[[crit]]))
    colnames(sig_mat) <- sprintf("%s_k%d", crit, seq_len(K_m1))
    n_signif <- rowSums(sig_mat)
    any_sig  <- n_signif > 0L

    # For median/CI summaries in the union, take the cutpoint with max |median_beta|.
    median_mat <- do.call(cbind, lapply(per_cutpoint_dfs, function(d) d$median_beta))
    pick_k <- max.col(abs(median_mat), ties.method = "first")
    pick_idx <- cbind(seq_len(V), pick_k)
    median_pick <- median_mat[pick_idx]
    ci_lo_mat <- do.call(cbind, lapply(per_cutpoint_dfs, function(d) d$ci_low))
    ci_hi_mat <- do.call(cbind, lapply(per_cutpoint_dfs, function(d) d$ci_high))
    ci_lo_pick <- ci_lo_mat[pick_idx]
    ci_hi_pick <- ci_hi_mat[pick_idx]

    union_df <- data.frame(
      variant_name  = variant_index$variant_name,
      position      = variant_index$position,
      gene          = gene_names,
      hgvs_c        = hgvs_c,
      hgvs_p        = hgvs_p,
      effect_impact = effect_imp,
      median_beta   = median_pick,
      ci_low        = ci_lo_pick,
      ci_high       = ci_hi_pick,
      cutpoint_of_max_abs_median = pick_k,
      n_signif_cutpoints = n_signif,
      stringsAsFactors = FALSE
    )
    union_df[[paste0("signif_", crit_short(crit))]] <- any_sig
    cbind(union_df, as.data.frame(sig_mat))
  }
  crit_short <- function(crit) sub("^signif_", "", crit)

  union_median <- build_union("signif_median")
  union_ci     <- build_union("signif_ci")

  write.csv(union_median[union_median$signif_median, , drop = FALSE],
            file.path(out_dir, "significant_median_union.csv"), row.names = FALSE)
  write.csv(union_ci[union_ci$signif_ci, , drop = FALSE],
            file.path(out_dir, "significant_ci_union.csv"), row.names = FALSE)
  message(sprintf("[rope] union (median): %d / %d", sum(union_median$signif_median), V))
  message(sprintf("[rope] union (CI):     %d / %d", sum(union_ci$signif_ci), V))

  # Union plots: subtitle says criterion + "any cutpoint"; labels gain count suffix.
  write_manhattan_union <- function(df, signif_col, criterion_label, out_path) {
    df_plot <- df
    use_gene <- !is.na(df_plot$gene) & df_plot$gene != "MODIFIER"
    base_lbl <- ifelse(use_gene, df_plot$gene,
                       format_variant_labels(df_plot$hgvs_p, df_plot$hgvs_c,
                                             df_plot$variant_name))
    df_plot$lbl_text <- sprintf("%s (k=%d)", base_lbl, df_plot$n_signif_cutpoints)
    pdf <- data.frame(
      pos    = as.numeric(df_plot$position),
      median = as.numeric(df_plot$median_beta),
      signif = ifelse(df_plot[[signif_col]], "significant", "not significant"),
      lbl    = df_plot$lbl_text,
      use_italic = use_gene,
      stringsAsFactors = FALSE
    )
    pdf$signif <- factor(pdf$signif, levels = c("not significant", "significant"))
    sig_colors <- c("not significant" = "grey60", "significant" = "red")

    label_df <- pdf[df_plot[[signif_col]], , drop = FALSE]
    if (nrow(label_df) > 0L) {
      label_df$gene_expr <- ifelse(
        label_df$use_italic,
        paste0("italic('", gsub("'", "", sub(" \\(k=.*$", "", label_df$lbl)), "')~'",
               sub("^[^ ]+ ", "", label_df$lbl), "'"),
        paste0("'", gsub("'", "", label_df$lbl), "'")
      )
    }

    subtitle <- sprintf("%s outside ROPE [%.4g, %.4g] (any cutpoint)",
                        criterion_label, rope[1], rope[2])

    p <- ggplot(pdf, aes(x = pos, y = exp(median), colour = signif)) +
      geom_point(alpha = 0.4) +
      scale_colour_manual(values = sig_colors) +
      xlab("genome coordinate (bp)") +
      ylab(expression("e"^{tilde(beta)})) +
      labs(colour = "ROPE", subtitle = subtitle) +
      theme_minimal(base_size = 14)
    if (nrow(label_df) > 0L) {
      p <- p + geom_text_repel(
        data = label_df,
        aes(x = pos, y = exp(median), label = gene_expr),
        parse = TRUE,
        size = 3,
        arrow = grid::arrow(length = grid::unit(0.01, "npc"), type = "open"),
        colour = "black",
        inherit.aes = FALSE,
        max.overlaps = Inf
      )
    }
    ggsave(out_path, plot = p, width = 16, height = 6, dpi = 300)
  }

  write_manhattan_union(union_median, "signif_median", "median effect",
                        file.path(out_dir, "manhattan_median_union.png"))
  write_manhattan_union(union_ci,     "signif_ci",
                        sprintf("%.0f%% CI", 100 * ci_prob),
                        file.path(out_dir, "manhattan_ci_union.png"))
}

message("[rope] done. outputs in: ", out_dir)
