library(maftools)

# read the annotated MAF from Funcotator
maf <- read.maf(maf = "results/variants/somatic.funcotated.maf")

# --- Figure 1: Oncoplot (top mutated genes) ---
png("results/figures/oncoplot.png", width = 1200, height = 900, res = 150)
oncoplot(maf = maf, top = 20)
dev.off()

# --- Figure 2: MAF summary (variant types, counts) ---
png("results/figures/maf_summary.png", width = 1200, height = 900, res = 150)
plotmafSummary(maf = maf, addStat = "median", dashboard = TRUE)
dev.off()

# print a text summary to the log
cat("\n=== SAMPLE SUMMARY ===\n")
print(getSampleSummary(maf))
cat("\n=== TOP MUTATED GENES ===\n")
print(getGeneSummary(maf)[1:20])

cat("\nOncoplot written to results/figures/oncoplot.png\n")
