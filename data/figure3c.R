combined_ttms_wgbs <- read.delim('combined_ttms_wgbs.txt')

smoothScatter(combined_ttms_wgbs$avg_meth.y, combined_ttms_wgbs$avg_meth.x, nrpoints = 1000,
              ylab = "TTMS", xlab = "WGBS", cex = 2)
abline(a = 0, b = 1, lty = 2)
