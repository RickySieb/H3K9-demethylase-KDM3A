## Track shot: ATAC WT vs KO, CTRL vs LPS signal at IFIH1 locus
## R. Siebeler

library(karyoploteR)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)  # for import of bigWig files

# ---- settings, edit before running ----
BW_DIR <- "/path/to/avgBigWig"   # directory containing averaged BigWig files
REGION <- "chr2:62597447-62646153"
REGION_PAD <- 10000

# ---- helper: pad a GRanges region on both sides ----
pad_region <- function(region, pad) {
  start(region) <- start(region) - pad
  end(region) <- end(region) + pad
  region
}

IFIH1.region <- pad_region(toGRanges(REGION), pad = REGION_PAD)

# ---- BigWig tracks used in this trackshot ----
bw.ATAC_WT     <- file.path(BW_DIR, "ATAC_WT_averaged.bw")
bw.ATAC_KO     <- file.path(BW_DIR, "ATAC_KO_averaged.bw")
bw.ATAC_WTLPS  <- file.path(BW_DIR, "ATAC_WT_LPS_averaged.bw")
bw.ATAC_KOLPS  <- file.path(BW_DIR, "ATAC_KO_LPS_averaged.bw")

tracks <- list(
  list(data = bw.ATAC_WT,    ymax = 10, label = "ATAC WT",     col = "#1e8239"),
  list(data = bw.ATAC_KO,    ymax = 10, label = "ATAC KO",     col = "#1e8239"),
  list(data = bw.ATAC_WTLPS, ymax = 10, label = "ATAC WT LPS", col = "#e0720d"),
  list(data = bw.ATAC_KOLPS, ymax = 10, label = "ATAC KO LPS", col = "#e0720d")
)
tracks <- rev(tracks)

# ---- gene model track ----
kp <- plotKaryotype(zoom = IFIH1.region, cex = 2)
genes.data <- makeGenesDataFromTxDb(TxDb.Mmusculus.UCSC.mm10.knownGene,
                                    karyoplot = kp,
                                    plot.transcripts = TRUE,
                                    plot.transcripts.structure = TRUE)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)

pp <- getDefaultPlotParams(plot.type = 1)
pp$leftmargin <- 0.15
pp$topmargin <- 15
pp$bottommargin <- 15
pp$ideogramheight <- 1
pp$data1inmargin <- 10
pp$data.panel.height <- 0.3

kp <- plotKaryotype(zoom = IFIH1.region, cex = 0.1, plot.params = pp)
kpPlotGenes(kp, data = genes.data, r0 = 0, r1 = 0.03, gene.name.cex = 1)

# ---- BigWig signal tracks ----
track_height <- 0.04
spacing <- 0.028
start_r <- 0.075  # vertical position where BigWig tracks start

for (i in seq_along(tracks)) {
  r0 <- start_r + (i - 1) * (track_height + spacing)
  r1 <- r0 + track_height

  tr <- tracks[[i]]

  kpPlotBigWig(kp, data = tr$data, ymax = tr$ymax, r0 = r0, r1 = r1, col = tr$col)
  kpAxis(kp, ymin = 0, ymax = tr$ymax, r0 = r0, r1 = r1, cex = .75)
  kpAddLabels(kp, labels = tr$label, r0 = r0, r1 = r1, cex = .75, label.margin = 0.07)
}
