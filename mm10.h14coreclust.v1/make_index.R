library(GenomicRanges)
library(AnnotationHub)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TFBSTools)
library(bacha)
library(Signac)

source("fetch.R")
options(timeout = 10000)
set.seed(1234)

# load open chromatin regions
fetch_file("https://downloads.wenglab.org/Registry-V4/mm10-cCREs.bed",
           destfile = "mm10-cCREs.bed")
ocr <- read.table("mm10-cCREs.bed", sep = "\t")
ocr <- makeGRangesFromDataFrame(ocr, seqnames.field = "V1", start.field = "V2", end.field = "V3")

# load regions to exclude
ah <- AnnotationHub()
excludeGR.mm10 <- ah[["AH107322"]]
repeat_regions <- ah[['AH99012']]
mask_regions <- c(excludeGR.mm10, repeat_regions)
mask_regions <- reduce(mask_regions)
strand(mask_regions) <- "*"

# load PWMs
if (!dir.exists("H14CORE-CLUSTERED")) {
  fetch_file("https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE-CLUSTERED/H14CORE-CLUSTERED_pwm.tar.gz",
             destfile = "H14CORE-CLUSTERED_pwm.tar.gz")
  untar("H14CORE-CLUSTERED_pwm.tar.gz")
}

pwm <- ReadPWM("H14CORE-CLUSTERED/pwm")

# restrict to TFs that exist in mouse; H14CORE covers human and mouse, and
# scanning the mouse genome for human-only TFs gives features with no mouse
# transcription factor behind them (see ../hocomoco_mouse.R)
source("hocomoco_mouse.R")
pwm <- subset_mouse_motifs(pwm)

# create Bacha index
BachaIndex(
  pwm = pwm,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  open_regions = ocr,
  mask_regions = mask_regions,
  tss_positions = promoters(TxDb.Mmusculus.UCSC.mm10.knownGene, upstream = 0, downstream = 1),
  match_features = c("GC", "entropy", "tssDist"),
  metadata = c(
    motif_source = "https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE-CLUSTERED/H14CORE-CLUSTERED_pwm.tar.gz (ReadPWM); restricted to TFs found in mouse (H14CORE_annotation.jsonl)",
    open_regions_source = "ENCODE cCREs v4, https://downloads.wenglab.org/Registry-V4/mm10-cCREs.bed",
    mask_regions_source = "AnnotationHub AH107322 (exclude ranges) + AH99012 (repeat regions)",
    tss_source = paste("TxDb.Mmusculus.UCSC.mm10.knownGene", packageVersion("TxDb.Mmusculus.UCSC.mm10.knownGene"), "transcript starts")
  ),
  output_path = "index/mm10.h14coreclust.v1.bi",
  motif_set_name = "H14CORE_CLUSTERED_MOUSE"
)
