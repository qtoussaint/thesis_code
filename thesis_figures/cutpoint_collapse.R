# ------------------------------------------------------------------
# Cutpoint collapse vs. well-separated cutpoints in a cumulative
# (ordered logistic) model -- styled to match the parallel-lines
# (proportional odds) diagram so the two figures sit together in
# one paper.
#
# 2 x 2 layout:
#   columns : collapsed cutpoints (left) | well-separated cutpoints (right)
#   rows    : cutpoint posteriors on the latent scale (top)
#             implied category probabilities with 95% CrI (bottom)
#
# Style conventions shared with po_diagram():
#   * base graphics, bty = "l", las = 1
#   * hcl.colors(., "Dark 3") line palette
#   * colour-matched plotmath labels (italic(c)[k]) placed on the plot
#     itself rather than in a legend box
#   * faint grey dotted/dashed reference guides
#   * bold per-panel titles + one bold outer title via mtext
# ------------------------------------------------------------------

cutpoint_collapse_diagram <- function(collapse_draws, good_draws,
                                      xlim = c(-4.2, 4.2), n = 400) {
  
  stopifnot(ncol(collapse_draws) == ncol(good_draws))
  
  Jc   <- ncol(collapse_draws)            # number of cutpoints (K - 1)
  K    <- Jc + 1
  cols <- hcl.colors(Jc, palette = "Dark 3")
  
  # ---- densities of every cutpoint, both scenarios -----------------
  dens_list <- function(d) lapply(seq_len(Jc), function(k) density(d[[k]]))
  dens_col  <- dens_list(collapse_draws)
  dens_good <- dens_list(good_draws)
  
  # common y-range for the two top panels, so widths/heights of the
  # posteriors are visually comparable across columns
  dens_ymax <- max(vapply(c(dens_col, dens_good), function(d) max(d$y), 0))
  dens_ylim <- c(0, dens_ymax * 1.12)
  
  # latent standard-logistic backdrop, scaled to sit under the densities
  xs      <- seq(xlim[1], xlim[2], length.out = n)
  latent  <- dlogis(xs) * (0.45 * dens_ymax / dlogis(0))
  
  # spread label positions apart so neighbouring labels keep a minimum
  # gap (same helper as in po_diagram)
  spread_labels <- function(y, mingap) {
    ord <- order(y)
    ys  <- y[ord]
    for (i in seq_along(ys)[-1])
      if (ys[i] - ys[i - 1] < mingap) ys[i] <- ys[i - 1] + mingap
    out <- y
    out[ord] <- ys
    out
  }
  
  # ---- top-row panel: cutpoint posteriors ---------------------------
  dens_panel <- function(dl, main) {
    plot(NA, xlim = xlim, ylim = dens_ylim,
         xlab = expression("Latent scale" ~ (italic(z))),
         ylab = "Posterior density",
         main = main, las = 1, bty = "l")
    
    # grey latent backdrop + zero guide, echoing the guide style of the
    # parallel-lines figure
    polygon(c(xs, rev(xs)), c(latent, rep(0, n)),
            col = "grey88", border = NA)
    abline(v = 0, lty = 2, col = "grey80")
    
    for (k in seq_len(Jc)) {
      d <- dl[[k]]
      polygon(c(d$x, rev(d$x)), c(d$y, rep(0, length(d$y))),
              col = adjustcolor(cols[k], alpha.f = 0.35), border = NA)
      lines(d$x, d$y, col = cols[k], lwd = 2.4)
    }
    
    # colour-matched italic(c)[k] labels above each density peak,
    # nudged horizontally apart if two peaks nearly coincide
    peak_x <- vapply(dl, function(d) d$x[which.max(d$y)], 0)
    peak_y <- vapply(dl, function(d) max(d$y), 0)
    lab_x  <- spread_labels(peak_x, 0.05 * diff(xlim))
    text(lab_x, pmin(peak_y + 0.05 * dens_ymax, dens_ylim[2]),
         labels = parse(text = paste0("italic(c)[", seq_len(Jc), "]")),
         col = cols, cex = 1, xpd = NA)
  }
  
  # ---- bottom-row panel: implied category probabilities -------------
  prob_summary <- function(d) {
    cum <- cbind(0, plogis(as.matrix(d)), 1)
    p   <- cum[, -1, drop = FALSE] - cum[, -ncol(cum), drop = FALSE]
    q   <- apply(p, 2, quantile, probs = c(0.025, 0.975))
    list(mean = colMeans(p), lo = q[1, ], hi = q[2, ])
  }
  
  ps_col  <- prob_summary(collapse_draws)
  ps_good <- prob_summary(good_draws)
  prob_ymax <- max(ps_col$hi, ps_good$hi) * 1.22
  
  prob_panel <- function(ps, main) {
    mid <- barplot(ps$mean, names.arg = seq_len(K),
                   col = "grey45", border = NA,
                   ylim = c(0, prob_ymax),
                   xlab = expression("Ordinal category" ~ (italic(k))),
                   ylab = expression(plain(P) * group("(", italic(Y) == italic(k), ")")),
                   main = main, las = 1)
    box(bty = "l")
    arrows(mid, ps$lo, mid, ps$hi,
           angle = 90, code = 3, length = 0.03, lwd = 1.2)
    text(mid, ps$hi + 0.035 * prob_ymax,
         labels = sprintf("%.2f", ps$mean), cex = 0.8, xpd = NA)
  }
  
  op <- par(mfrow = c(2, 2), mar = c(4.2, 4.6, 3, 2.4), oma = c(0, 0, 2.5, 0),
            cex.main = 1.05, font.main = 2)
  on.exit(par(op))
  
  dens_panel(dens_col,  "Collapsed cutpoints (flat prior)")
  dens_panel(dens_good, "Well-separated cutpoints (informative priors)")
  prob_panel(ps_col,    "Implied category probabilities (flat prior)")
  prob_panel(ps_good,   "Implied category probabilities (informative priors)")
  mtext("Simulation of cutpoint collapse under a flat cutpoint prior",
        outer = TRUE, cex = 1.15, font = 2)
  invisible(NULL)
}

# ------------------------------------------------------------------
# Inputs: simulated posterior draws
# (swap in as_draws_df() output from real Stan fits -- columns c1..c4)
# ------------------------------------------------------------------
set.seed(42)
n_draws <- 4000

collapse_draws <- data.frame(
  c1 = rnorm(n_draws, -1.60, 0.22),
  c2 = rnorm(n_draws,  0.28, 0.45)
)
collapse_draws$c3 <- collapse_draws$c2 + abs(rnorm(n_draws, 0.03, 0.05))
collapse_draws$c4 <- rnorm(n_draws, 1.85, 0.24)

good_draws <- data.frame(
  c1 = rnorm(n_draws, -2.10, 0.15),
  c2 = rnorm(n_draws, -0.70, 0.13),
  c3 = rnorm(n_draws,  0.70, 0.13),
  c4 = rnorm(n_draws,  2.10, 0.15)
)

# draw to screen
cutpoint_collapse_diagram(collapse_draws, good_draws)

# save for the thesis (uncomment as needed)
# pdf("cutpoint_collapse.pdf", width = 11, height = 9)
# cutpoint_collapse_diagram(collapse_draws, good_draws); dev.off()
png("~/Desktop/cutpoint_collapse.png", width = 11, height = 9, units = "in", res = 600)
cutpoint_collapse_diagram(collapse_draws, good_draws)
dev.off()
