# ------------------------------------------------------------------
# Parallel-lines (proportional odds) diagram for an ordered model
# Model:  logit[P(Y <= k | x)] = alpha_k + beta_k * x
#
# 2 x 2 layout:
#   columns : proportional odds (left) | non-proportional odds (right)
#   rows    : log-odds scale (top)     | probability scale (bottom)
#
# Lines are labelled beta_k. Cutpoints c_k = alpha_k are drawn as faint
# horizontal guides on the LOG-ODDS panels only (a cutpoint is a threshold
# on the log-odds scale, so it is a horizontal reference at height alpha_k,
# not a vertical line through the x-axis). They are omitted on the
# probability panels, where alpha_k has no threshold interpretation.
#
# The two log-odds panels share a common y-range (computed from BOTH
# coefficient sets) so slopes are visually comparable across columns.
# ------------------------------------------------------------------

po_diagram <- function(alpha, beta_prop, beta_nonprop,
                       xlim = c(-3.2, 3.2), n = 400) {
  
  stopifnot(length(alpha) == length(beta_prop),
            length(alpha) == length(beta_nonprop))
  
  x     <- seq(xlim[1], xlim[2], length.out = n)
  J     <- length(alpha)
  cols  <- hcl.colors(J, palette = "Dark 3")
  expit <- function(z) 1 / (1 + exp(-z))
  
  # cumulative logits: returns an n x J matrix
  clogits <- function(beta) sapply(seq_len(J), function(k) alpha[k] + beta[k] * x)
  
  # common y-limits for the two log-odds panels, so the parallel-vs-fanning
  # contrast is not distorted by per-panel scaling
  logit_ylim <- range(clogits(beta_prop), clogits(beta_nonprop)) + c(-0.4, 0.4)
  
  # spread a set of y-positions apart so labels keep a minimum vertical gap,
  # preserving each element's identity (only the position is nudged upward)
  spread_labels <- function(y, mingap) {
    ord <- order(y)
    ys  <- y[ord]
    for (i in seq_along(ys)[-1])
      if (ys[i] - ys[i - 1] < mingap) ys[i] <- ys[i - 1] + mingap
    out <- y
    out[ord] <- ys
    out
  }
  
  # beta_k labels at the right-hand end of each line, colour-matched,
  # nudged apart vertically if two lines end too close together
  label_lines <- function(yvals) {
    yend   <- yvals[n, ]
    mingap <- 0.04 * diff(par("usr")[3:4])
    ypos   <- spread_labels(yend, mingap)
    text(x[n], ypos,
         labels = parse(text = paste0("beta[", seq_len(J), "]")),
         pos = 4, xpd = NA, col = cols, cex = 0.9)
  }
  
  # cutpoints c_k = alpha_k : horizontal guides on the log-odds scale only,
  # each running from the left margin to x = 0 (where the cumulative logit
  # equals its intercept alpha_k) and meeting a dashed vertical line at x = 0.
  # Sorted so the smallest alpha is labelled c_1. Label heights are spread
  # apart when neighbouring cutpoints sit too close, and lifted slightly upward.
  cutpoint_guides <- function() {
    ord      <- order(alpha)                 # smallest alpha -> c_1
    a_sorted <- alpha[ord]
    abline(v = 0, lty = 2, col = "grey80")    # x = 0: curves attain alpha_k here
    for (rank in seq_along(ord)) {
      a <- a_sorted[rank]
      segments(xlim[1], a, 0, a, lty = 3, col = "grey55", lwd = 0.8)
    }
    yr     <- diff(par("usr")[3:4])
    ylab   <- spread_labels(a_sorted, 0.045 * yr)
    yoff   <- 0.028 * yr                       # lift labels a bit higher
    xoff   <- xlim[1] + 0.04 * diff(xlim)
    for (rank in seq_along(ord))
      text(xoff, ylab[rank] + yoff,
           labels = parse(text = paste0("italic(c)[", rank, "]")),
           pos = 2, xpd = NA, cex = 1, col = "grey30")
  }
  
  one_panel <- function(beta, scale = c("logit", "prob"), main) {
    scale <- match.arg(scale)
    L <- clogits(beta)
    Y <- if (scale == "logit") L else expit(L)
    ylim <- if (scale == "logit") logit_ylim else c(0, 1)
    ylab <- if (scale == "logit")
      expression("logit" ~ plain(P) * group("(", Y <= k ~ "|" ~ x, ")"))
    else
      expression(plain(P) * group("(", Y <= k ~ "|" ~ x, ")"))
    
    plot(NA, xlim = xlim, ylim = ylim,
         xlab = expression("Predictor" ~ (italic(x))), ylab = ylab,
         main = main, las = 1, bty = "l")
    
    if (scale == "logit") {
      cutpoint_guides()
    } else {
      abline(h = 0.5, col = "grey80", lty = 2)
    }
    
    for (k in seq_len(J)) lines(x, Y[, k], col = cols[k], lwd = 2.4)
    label_lines(Y)
  }
  
  op <- par(mfrow = c(2, 2), mar = c(4.2, 4.6, 3, 2.4), oma = c(0, 0, 2.5, 0),
            cex.main = 1.05, font.main = 2)
  on.exit(par(op))
  
  one_panel(beta_prop,    "logit", "Proportional odds (log-odds of cumulative probability)")
  one_panel(beta_nonprop, "logit", "Non-proportional odds (log-odds of cumulative probability)")
  one_panel(beta_prop,    "prob",  "Proportional odds (cumulative probability)")
  one_panel(beta_nonprop, "prob",  "Non-proportional odds (cumulative probability)")
  mtext("Test of parallel lines for identification of non-proportional odds",
        outer = TRUE, cex = 1.15, font = 2)
  invisible(NULL)
}

# ------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------
alpha        <- c(-2.5, -1, 0, 1.5, 4)
beta_prop    <- c(0.90, 0.90, 0.90, 0.90, 0.90)
beta_nonprop <- c(0.15, 0.65, 1.20, 1.65, 0.85)

# draw to screen
po_diagram(alpha, beta_prop, beta_nonprop)

# save for the thesis (uncomment as needed)
# pdf("po_diagram.pdf", width = 11, height = 9)
# po_diagram(alpha, beta_prop, beta_nonprop); dev.off()
png("~/Desktop/parallel_lines.png", width = 11, height = 9, units = "in", res = 600)
po_diagram(alpha, beta_prop, beta_nonprop)
dev.off()
