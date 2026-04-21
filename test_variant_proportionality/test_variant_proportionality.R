#!/usr/bin/env Rscript
# Standalone comparison of variant effects under POM vs PPOM model specifications.
#
# Reads paired POM/PPOM pipeline output directories, computes per-variant
# Delta(beta) across adjacent cutpoints from PPOM posterior draws, checks
# containment of POM beta within each PPOM cutpoint's credible interval, and
# produces per-variant + genome-wide proportionality diagnostics.
#
# Operates on the de-pruned variant set: pruned variants inherit their
# representative's draws (same convention as depruned_variant_effects.csv).

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(data.table)
})

script_dir <- tryCatch({
  this_file <- sys.frame(1)$ofile
  if (is.null(this_file)) {
    args_cmd <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args_cmd, value = TRUE)
    if (length(file_arg) > 0L) sub("^--file=", "", file_arg[1]) else "."
  } else {
    this_file
  }
}, error = function(e) ".")
script_dir <- dirname(normalizePath(script_dir, mustWork = FALSE))
source(file.path(script_dir, "lib", "helpers.R"))

# ---- CLI -------------------------------------------------------------------

option_list <- list(
  make_option("--pom_dir",     type = "character", help = "POM pipeline output directory"),
  make_option("--ppom_dir",    type = "character", help = "PPOM pipeline output directory"),
  make_option("--out_dir",     type = "character", help = "Directory to write comparison outputs"),
  make_option("--ci_level",    type = "double", default = 0.89,
              help = "Credible interval level (default %default)"),
  make_option("--top_n_plots", type = "integer", default = 20,
              help = "Overlay plots for the top-N most PO-violating variants (default %default)"),
  make_option("--mic_csv",     type = "character", default = NULL,
              help = "Optional: use this depruned_variant_effects.csv to label cutpoint_MIC breakpoints. Defaults to <ppom_dir>/fitted_model/depruned_variant_effects.csv if present.")
)

opt <- parse_args(OptionParser(option_list = option_list))
for (req in c("pom_dir", "ppom_dir", "out_dir")) {
  if (is.null(opt[[req]])) stop("Missing required argument: --", req)
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$out_dir, "overlay_plots"), recursive = TRUE, showWarnings = FALSE)

probs <- ci_tail_probs(opt$ci_level)  # c(alpha/2, 0.5, 1-alpha/2)
lo_q_name <- paste0("q", probs[1])
md_q_name <- paste0("q", probs[2])
hi_q_name <- paste0("q", probs[3])

# ---- Load pruning metadata -------------------------------------------------

pom_ld_summary  <- file.path(opt$pom_dir,  "cppRATE_matrices", "ld_pruning_summary.csv")
ppom_ld_summary <- file.path(opt$ppom_dir, "cppRATE_matrices", "ld_pruning_summary.csv")

if (!file.exists(pom_ld_summary)) {
  stop("LD pruning summary not found for POM: ", pom_ld_summary)
}
if (!file.exists(ppom_ld_summary)) {
  stop("LD pruning summary not found for PPOM: ", ppom_ld_summary)
}

pom_final_rep  <- parse_ld_summary(pom_ld_summary)
ppom_final_rep <- parse_ld_summary(ppom_ld_summary)

pom_kept  <- read_kept_variants(opt$pom_dir)
ppom_kept <- read_kept_variants(opt$ppom_dir)

if (!identical(sort(pom_kept), sort(ppom_kept))) {
  warning("POM and PPOM kept-variant sets differ. Using PPOM's kept set as the canonical reference. ",
          "This usually means the two runs were fit with different LD thresholds; results may not be directly comparable.")
}

# Total variant count: representatives + unique pruned ids (1..V_total).
# Derive from PPOM side (canonical).
pruned_ppom <- unique(na.omit(
  as.vector(as.matrix(ppom_final_rep[, 2:ncol(ppom_final_rep)]))
))
V_total <- length(ppom_kept) + length(pruned_ppom)
# Sanity: variants are 1-indexed contiguous, so max id must equal V_total.
max_id <- max(max(ppom_kept), max(pruned_ppom, na.rm = TRUE))
if (max_id != V_total) {
  warning("Inferred V_total (", V_total, ") does not equal max variant id (",
          max_id, "). Using max id.")
  V_total <- max_id
}

# ---- Load draws ------------------------------------------------------------

pom_rds  <- find_fit_rds(opt$pom_dir)
ppom_rds <- find_fit_rds(opt$ppom_dir)

message("Loading POM draws from:  ", pom_rds)
pom_draws_kept <- load_beta_variant_draws(pom_rds)
message("Loading PPOM draws from: ", ppom_rds)
ppom_draws_kept <- load_beta_variant_draws(ppom_rds)

if (ncol(pom_draws_kept) != length(pom_kept)) {
  stop("POM draws have ", ncol(pom_draws_kept), " columns but ", length(pom_kept), " kept variants were expected.")
}
V_kept_ppom <- length(ppom_kept)
if (ncol(ppom_draws_kept) %% V_kept_ppom != 0L) {
  stop("PPOM draw columns (", ncol(ppom_draws_kept),
       ") not a multiple of kept-variant count (", V_kept_ppom, ").")
}
n_cutpoints <- ncol(ppom_draws_kept) %/% V_kept_ppom
message("Detected K - 1 = ", n_cutpoints, " cutpoints.")

# ---- De-prune draws to V_total space --------------------------------------

message("De-pruning POM draws to all ", V_total, " variants...")
pom_draws_full <- deprun_draws_single(
  draws_kept         = pom_draws_kept,
  kept_variants      = pom_kept,
  final_rep_variants = pom_final_rep,
  n_total_variants   = V_total
)

message("De-pruning PPOM draws to all ", V_total, " variants (x ", n_cutpoints, " cutpoints)...")
ppom_draws_full <- deprun_draws_ppom(
  ppom_draws         = ppom_draws_kept,
  kept_variants      = ppom_kept,
  final_rep_variants = ppom_final_rep,
  n_total_variants   = V_total
)  # list of length n_cutpoints, each [draws × V_total]

# ---- Per-variant summary stats --------------------------------------------

message("Computing POM per-variant credible intervals...")
pom_q <- safe_column_quantiles(pom_draws_full, probs = probs)
pom_median <- pom_q[, md_q_name]
pom_ci_lo  <- pom_q[, lo_q_name]
pom_ci_hi  <- pom_q[, hi_q_name]

message("Computing PPOM per-cutpoint credible intervals...")
ppom_qs <- lapply(ppom_draws_full, safe_column_quantiles, probs = probs)
# Matrices [V_total × n_cutpoints]
ppom_median <- do.call(cbind, lapply(ppom_qs, function(q) q[, md_q_name]))
ppom_ci_lo  <- do.call(cbind, lapply(ppom_qs, function(q) q[, lo_q_name]))
ppom_ci_hi  <- do.call(cbind, lapply(ppom_qs, function(q) q[, hi_q_name]))
colnames(ppom_median) <- paste0("ppom_median_k", seq_len(n_cutpoints))
colnames(ppom_ci_lo)  <- paste0("ppom_ci_lo_k",  seq_len(n_cutpoints))
colnames(ppom_ci_hi)  <- paste0("ppom_ci_hi_k",  seq_len(n_cutpoints))

# ---- Containment check ----------------------------------------------------

in_ci <- (ppom_ci_lo <= pom_median) & (pom_median <= ppom_ci_hi)
n_cis_covering <- rowSums(in_ci, na.rm = TRUE)
pom_in_all <- rowSums(in_ci, na.rm = TRUE) == n_cutpoints

# ---- Delta(beta) across adjacent cutpoints --------------------------------

message("Computing Delta(beta) draws and CIs for adjacent cutpoint pairs...")
n_pairs <- n_cutpoints - 1L
delta_median <- matrix(NA_real_, nrow = V_total, ncol = n_pairs)
delta_ci_lo  <- matrix(NA_real_, nrow = V_total, ncol = n_pairs)
delta_ci_hi  <- matrix(NA_real_, nrow = V_total, ncol = n_pairs)
po_violation <- matrix(FALSE,    nrow = V_total, ncol = n_pairs)
delta_abs_median <- matrix(NA_real_, nrow = V_total, ncol = n_pairs)

for (p in seq_len(n_pairs)) {
  diff_draws <- ppom_draws_full[[p + 1L]] - ppom_draws_full[[p]]
  q <- safe_column_quantiles(diff_draws, probs = probs)
  delta_median[, p] <- q[, md_q_name]
  delta_ci_lo[,  p] <- q[, lo_q_name]
  delta_ci_hi[,  p] <- q[, hi_q_name]
  po_violation[, p] <- (q[, lo_q_name] > 0) | (q[, hi_q_name] < 0)
  delta_abs_median[, p] <- abs(q[, md_q_name])
}

max_abs_delta <- apply(delta_abs_median, 1, max, na.rm = TRUE)
n_violations  <- rowSums(po_violation, na.rm = TRUE)

# ---- Representative metadata ----------------------------------------------

rep_map <- build_representative_map(ppom_final_rep, V_total)

# ---- Sanity check: PPOM medians must match depruned_variant_effects.csv ---

eff_csv_default <- file.path(opt$ppom_dir, "fitted_model", "depruned_variant_effects.csv")
eff_csv <- if (!is.null(opt$mic_csv)) opt$mic_csv else eff_csv_default
cutpoint_mic_labels <- setNames(as.character(seq_len(n_cutpoints)), seq_len(n_cutpoints))

if (file.exists(eff_csv)) {
  message("Cross-checking PPOM medians against ", eff_csv, "...")
  eff <- data.table::fread(eff_csv)
  # Pull MIC breakpoint labels if present.
  if ("cutpoint_MIC" %in% names(eff) && "cutpoint" %in% names(eff)) {
    by_cp <- unique(eff[, .(cutpoint, cutpoint_MIC)])
    by_cp <- by_cp[order(cutpoint)]
    if (nrow(by_cp) == n_cutpoints) {
      cutpoint_mic_labels <- setNames(as.character(by_cp$cutpoint_MIC), by_cp$cutpoint)
    }
  }
  # Verify median agreement for cutpoint 1 (representatives only, as a sanity check).
  eff_c1 <- eff[cutpoint == 1]
  comp <- data.frame(
    variant_id = eff_c1$variant_id,
    csv_median = eff_c1$median,
    our_median = ppom_median[eff_c1$variant_id, 1]
  )
  comp_nonrep <- comp[rep_map$is_representative[comp$variant_id], ]
  # Both should use the same draws; tolerate small numeric noise.
  max_diff <- max(abs(comp_nonrep$csv_median - comp_nonrep$our_median), na.rm = TRUE)
  if (is.finite(max_diff) && max_diff > 1e-6) {
    warning("Max abs difference between our PPOM median (cutpoint 1) and depruned_variant_effects.csv: ",
            format(max_diff, digits = 3),
            ". This may indicate a variant-indexing mismatch or different draws.")
  } else {
    message("  OK: max abs difference = ", format(max_diff, digits = 3))
  }
} else {
  message("No depruned_variant_effects.csv found for cross-check (looked at ", eff_csv, ")")
}

# ---- Assemble and write per-variant CSV ----------------------------------

message("Writing per_variant_comparison.csv...")
per_variant <- data.frame(
  variant_id        = rep_map$variant_id,
  is_representative = rep_map$is_representative,
  representative_id = rep_map$representative_id,
  pom_median        = pom_median,
  pom_ci_lo         = pom_ci_lo,
  pom_ci_hi         = pom_ci_hi
)
per_variant <- cbind(per_variant, ppom_median, ppom_ci_lo, ppom_ci_hi)
per_variant$n_ppom_cis_covering_pom <- n_cis_covering
per_variant$pom_in_all_ppom_cis     <- pom_in_all
per_variant$n_po_violations         <- n_violations
per_variant$max_abs_delta_beta      <- ifelse(is.infinite(max_abs_delta), NA_real_, max_abs_delta)

data.table::fwrite(per_variant,
                   file.path(opt$out_dir, "per_variant_comparison.csv"),
                   row.names = FALSE)

# ---- Delta(beta) long-format CSV -----------------------------------------

message("Writing delta_beta_summary.csv...")
pair_labels <- paste0(seq_len(n_pairs), "->", seq_len(n_pairs) + 1L)
delta_long <- data.frame(
  variant_id   = rep(seq_len(V_total), times = n_pairs),
  cutpoint_pair = rep(pair_labels, each = V_total),
  delta_median = as.vector(delta_median),
  delta_ci_lo  = as.vector(delta_ci_lo),
  delta_ci_hi  = as.vector(delta_ci_hi),
  po_violation = as.vector(po_violation)
)
data.table::fwrite(delta_long,
                   file.path(opt$out_dir, "delta_beta_summary.csv"),
                   row.names = FALSE)

# ---- Genome-wide summary --------------------------------------------------

message("Writing genome_summary.txt...")
frac_violation_per_pair <- colMeans(po_violation, na.rm = TRUE)
frac_any_violation <- mean(n_violations > 0, na.rm = TRUE)
frac_pom_in_all    <- mean(pom_in_all, na.rm = TRUE)
abs_deltas <- as.vector(delta_abs_median)
summary_lines <- c(
  "POM vs PPOM proportionality comparison",
  paste0("CI level: ", opt$ci_level),
  paste0("Total de-pruned variants: ", V_total),
  paste0("Representatives: ", sum(rep_map$is_representative)),
  paste0("Cutpoints (K - 1): ", n_cutpoints),
  paste0("Cutpoint-pairs analysed: ", n_pairs),
  "",
  "-- Per-pair violation rates (fraction of variants where CI of Delta(beta) excludes 0) --"
)
for (p in seq_len(n_pairs)) {
  summary_lines <- c(
    summary_lines,
    sprintf("  %s (%s -> %s):  %.4f",
            pair_labels[p],
            cutpoint_mic_labels[as.character(p)],
            cutpoint_mic_labels[as.character(p + 1L)],
            frac_violation_per_pair[p])
  )
}
summary_lines <- c(
  summary_lines,
  "",
  sprintf("Fraction of variants with >=1 pair violating PO: %.4f", frac_any_violation),
  sprintf("Fraction of variants where POM median lies inside ALL PPOM CIs: %.4f", frac_pom_in_all),
  "",
  "-- |Delta(beta)| distribution across all variants x pairs --",
  sprintf("  median: %.4f", stats::median(abs_deltas, na.rm = TRUE)),
  sprintf("  95th pct: %.4f", stats::quantile(abs_deltas, 0.95, na.rm = TRUE)),
  sprintf("  max:    %.4f", max(abs_deltas, na.rm = TRUE))
)
writeLines(summary_lines, file.path(opt$out_dir, "genome_summary.txt"))
message("")
writeLines(summary_lines)
message("")

# ---- Manhattan-style summary plot ----------------------------------------

message("Writing manhattan_po_violation.png...")
# Continuous violation score = sum over pairs of |delta_median| / width(CI)/2 (a z-like score).
delta_half_width <- (delta_ci_hi - delta_ci_lo) / 2
eps <- 1e-12
viol_score <- rowSums(delta_abs_median / (delta_half_width + eps), na.rm = TRUE)

manhattan_df <- data.frame(
  variant_id = seq_len(V_total),
  score      = viol_score,
  n_viol     = n_violations,
  is_rep     = rep_map$is_representative
)

p_manhattan <- ggplot(manhattan_df, aes(x = variant_id, y = score, colour = factor(n_viol))) +
  geom_point(size = 0.6, alpha = 0.8) +
  labs(x = "Variant ID (1..V, de-pruned)",
       y = "PO-violation score (sum |Delta(beta)_median| / CI half-width)",
       colour = "# pairs violating PO",
       title = "Genome-wide PO violation score",
       subtitle = sprintf("%.1f%% of variants violate PO at >=1 pair; %.1f%% have POM median in all PPOM CIs",
                          100 * frac_any_violation, 100 * frac_pom_in_all)) +
  theme_bw() +
  theme(legend.position = "right")

ggsave(file.path(opt$out_dir, "manhattan_po_violation.png"),
       plot = p_manhattan, width = 10, height = 4.5, dpi = 150)

# ---- Overlay plots for top-N violating variants --------------------------

message("Writing overlay plots for top ", opt$top_n_plots, " most PO-violating variants...")
ranked <- order(viol_score, decreasing = TRUE, na.last = TRUE)
ranked <- ranked[!is.na(viol_score[ranked]) & viol_score[ranked] > 0]
top_ids <- head(ranked, opt$top_n_plots)

mic_levels_vec <- cutpoint_mic_labels[as.character(seq_len(n_cutpoints))]
x_labels <- paste0("k=", seq_len(n_cutpoints), "\n", mic_levels_vec)

for (v in top_ids) {
  ppom_row <- data.frame(
    cutpoint = seq_len(n_cutpoints),
    median   = ppom_median[v, ],
    lo       = ppom_ci_lo[v, ],
    hi       = ppom_ci_hi[v, ]
  )
  p_v <- ggplot(ppom_row, aes(x = cutpoint, y = median)) +
    annotate("rect",
             xmin = 0.5, xmax = n_cutpoints + 0.5,
             ymin = pom_ci_lo[v], ymax = pom_ci_hi[v],
             alpha = 0.15, fill = "steelblue") +
    geom_hline(yintercept = pom_median[v], colour = "steelblue",
               linetype = "solid", linewidth = 0.7) +
    geom_pointrange(aes(ymin = lo, ymax = hi), colour = "firebrick", size = 0.5) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = "dotted") +
    scale_x_continuous(breaks = seq_len(n_cutpoints), labels = x_labels) +
    labs(
      x = "Cutpoint (PPOM)",
      y = expression(beta),
      title = sprintf("Variant %d: POM vs PPOM beta (%.0f%% CI)",
                      v, 100 * opt$ci_level),
      subtitle = sprintf("rep=%s, POM median=%.3f; %d/%d pairs violate PO",
                         if (rep_map$is_representative[v]) "self" else as.character(rep_map$representative_id[v]),
                         pom_median[v], n_violations[v], n_pairs)
    ) +
    theme_bw()
  ggsave(file.path(opt$out_dir, "overlay_plots", sprintf("variant_%d.png", v)),
         plot = p_v, width = 6, height = 4, dpi = 150)
}

message("Done. Results written to ", opt$out_dir)
