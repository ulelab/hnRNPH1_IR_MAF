#!/usr/bin/env Rscript
# One-shot Part 2: merge phastCons + MaxEnt into step1, plot, write final TSV.
# Run from hnRNPH1_IR_MAF:  Rscript scripts/run_v3_part2_merge_plot.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

out_dir <- "results"
step1_tsv <- file.path(out_dir, "decoy_intron_overlap_step1.tsv")
ph100_tsv <- file.path(out_dir, "v3_phastcons100_by_decoy.tsv")
ph470_tsv <- file.path(out_dir, "v3_phastcons470_by_decoy.tsv")
maxent_decoy_bed <- file.path(out_dir, "v3_maxent_decoy.bed")
maxent_canonical_bed <- file.path(out_dir, "v3_maxent_canonical.bed")
final_tsv <- file.path(out_dir, "decoy_intron_features_v3_final.tsv")
plot_png <- file.path(out_dir, "decoy_intron_v3_plot.png")

if (!file.exists(step1_tsv)) stop("Missing: ", step1_tsv)

merged <- fread(step1_tsv, sep = "\t", header = TRUE)
message("Reloaded rows: ", nrow(merged))

if (file.exists(ph100_tsv)) {
  ph100 <- fread(ph100_tsv, sep = "\t", header = FALSE, col.names = c("decoyID", "phastCons100_mean"))
  merged <- merge(merged, ph100, by = "decoyID", all.x = TRUE)
} else merged[, phastCons100_mean := NA_real_]

if (file.exists(ph470_tsv)) {
  ph470 <- fread(ph470_tsv, sep = "\t", header = FALSE, col.names = c("decoyID", "phastCons470_mean"))
  merged <- merge(merged, ph470, by = "decoyID", all.x = TRUE)
} else merged[, phastCons470_mean := NA_real_]

if (file.exists(maxent_decoy_bed)) {
  raw <- fread(maxent_decoy_bed, sep = "\t", header = FALSE)
  sc <- unique(
    data.table(decoyID = raw[[4]], MaxEnt_decoy = as.numeric(raw[[ncol(raw)]])),
    by = "decoyID"
  )
  merged <- merge(merged, sc, by = "decoyID", all.x = TRUE)
} else merged[, MaxEnt_decoy := NA_real_]

if (file.exists(maxent_canonical_bed)) {
  raw <- fread(maxent_canonical_bed, sep = "\t", header = FALSE)
  sc <- unique(
    data.table(decoyID = raw[[4]], MaxEnt_canonical = as.numeric(raw[[ncol(raw)]])),
    by = "decoyID"
  )
  merged <- merge(merged, sc, by = "decoyID", all.x = TRUE)
} else merged[, MaxEnt_canonical := NA_real_]

meta_cols <- c(
  "GENE", "EVENT", "COORD", "LENGTH", "FullCO", "COMPLEX",
  "chr", "intron_start", "intron_end", "intron_length"
)
feat <- c(
  "decoy_coord_1based", "decoy_distance_from_canonical_5ss",
  "start", "end", "decoyID", "SpliceAI_score", "strand",
  "MaxEnt_decoy", "MaxEnt_canonical", "phastCons470_mean", "phastCons100_mean"
)
psi_cols <- setdiff(names(merged), c(meta_cols, feat, "chr.1"))
psi_num <- psi_cols[sapply(merged[, ..psi_cols], is.numeric)]
merged[, psi_count_ge10 := rowSums(.SD >= 10, na.rm = TRUE), .SDcols = psi_num]

bed_cols <- c("chr", "start", "end", "decoyID", "SpliceAI_score", "strand")
pref <- c(
  "decoy_coord_1based", "decoy_distance_from_canonical_5ss",
  "MaxEnt_decoy", "MaxEnt_canonical",
  "phastCons470_mean", "phastCons100_mean", "psi_count_ge10"
)
pref <- pref[pref %in% names(merged)]
setcolorder(merged, c(bed_cols, pref, setdiff(names(merged), c(bed_cols, pref))))

fwrite(merged, final_tsv, sep = "\t", quote = FALSE)
message("Wrote: ", final_tsv)

plot_df <- merged[!is.na(decoy_distance_from_canonical_5ss) & !is.na(psi_count_ge10)]
plot_df[, log_distance_from_5ss := log10(abs(decoy_distance_from_canonical_5ss) + 1)]
highlight_id <- "HNRNPH1_179620582"

p <- ggplot(plot_df, aes(
  x = log_distance_from_5ss,
  y = psi_count_ge10,
  size = SpliceAI_score,
  color = phastCons470_mean
)) +
  geom_point(alpha = 0.5) +
  geom_point(data = plot_df[decoyID == highlight_id], color = "red", alpha = 0.9, stroke = 1.2) +
  geom_text(
    data = plot_df[decoyID == highlight_id],
    aes(label = decoyID),
    color = "red",
    vjust = -0.8,
    size = 3,
    show.legend = FALSE
  ) +
  labs(
    title = "Log distance from canonical 5'SS vs PSI columns >= 10",
    x = "log10(|distance| + 1)",
    y = "Count of PSI columns >= 10",
    size = "SpliceAI score",
    color = "phastCons 470-way (intron mean)"
  ) +
  scale_color_viridis_c(option = "C", na.value = "grey80") +
  theme_bw()

ggsave(plot_png, p, width = 10, height = 7, dpi = 150)
message("Wrote plot: ", plot_png)
