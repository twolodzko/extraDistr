# Compare extraDistr::dnaka with nakagami::dnaka
#
# This script validates that the extraDistr implementation matches
# the CRAN nakagami package. Run after devtools::load_all().
#
# Note: nakagami uses (shape, scale), extraDistr uses (m, w) - same semantics.

if (!requireNamespace("nakagami", quietly = TRUE)) {
  stop("Install nakagami package: install.packages('nakagami')")
}

devtools::load_all()

# Test parameters
x <- seq(0.01, 5, by = 0.1)
shapes <- c(0.5, 1, 2, 5)
scales <- c(0.5, 1, 2)

cat("Comparing dnaka (density)...\n")
for (m in shapes) {
  for (w in scales) {
    extra <- dnaka(x, m = m, w = w)
    naka <- nakagami::dnaka(x, shape = m, scale = w)
    if (!isTRUE(all.equal(extra, naka, tolerance = 1e-10))) {
      stop(sprintf("dnaka mismatch at m=%g, w=%g", m, w))
    }
  }
}
cat("  OK\n")

cat("Comparing pnaka (CDF)...\n")
for (m in shapes) {
  for (w in scales) {
    extra <- pnaka(x, m = m, w = w)
    naka <- nakagami::pnaka(x, shape = m, scale = w)
    if (!isTRUE(all.equal(extra, naka, tolerance = 1e-10))) {
      stop(sprintf("pnaka mismatch at m=%g, w=%g", m, w))
    }
  }
}
cat("  OK\n")

cat("Comparing qnaka (quantile)...\n")
p <- seq(0.01, 0.99, by = 0.01)
for (m in shapes) {
  for (w in scales) {
    extra <- qnaka(p, m = m, w = w)
    naka <- nakagami::qnaka(p, shape = m, scale = w)
    if (!isTRUE(all.equal(extra, naka, tolerance = 1e-10))) {
      stop(sprintf("qnaka mismatch at m=%g, w=%g", m, w))
    }
  }
}
cat("  OK\n")

cat("Comparing rnaka (RNG distribution)...\n")
set.seed(42)
for (m in shapes) {
  for (w in scales) {
    extra <- rnaka(10000, m = m, w = w)
    naka <- nakagami::rnaka(10000, shape = m, scale = w)
    # Compare via KS test - both should follow same distribution
    ks <- ks.test(extra, naka)
    if (ks$p.value < 0.001) {
      stop(sprintf("rnaka distribution mismatch at m=%g, w=%g (p=%g)",
                   m, w, ks$p.value))
    }
  }
}
cat("  OK\n")

cat("\nAll comparisons passed!\n")
