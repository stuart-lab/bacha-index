library(GenomicRanges)
library(AnnotationHub)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TFBSTools)
library(bacha)
library(Signac)

source("fetch.R")
options(timeout = 10000)
set.seed(1234)

# load open chromatin regions
fetch_file("https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed",
           destfile = "GRCh38-cCREs.bed")
ocr <- read.table("GRCh38-cCREs.bed", sep = "\t")
ocr <- makeGRangesFromDataFrame(ocr, seqnames.field = "V1", start.field = "V2", end.field = "V3")

# load regions to exclude
ah <- AnnotationHub()
excludeGR.hg38 <- ah[["AH107305"]]
repeat_regions <- ah[['AH99003']]
mask_regions <- c(excludeGR.hg38, repeat_regions)
mask_regions <- reduce(mask_regions)
strand(mask_regions) <- "*"

# load PWMs
if (!dir.exists("H14CORE-CLUSTERED")) {
  fetch_file("https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE-CLUSTERED/H14CORE-CLUSTERED_pwm.tar.gz",
             destfile = "H14CORE-CLUSTERED_pwm.tar.gz")
  untar("H14CORE-CLUSTERED_pwm.tar.gz")
}

pwm <- ReadPWM("H14CORE-CLUSTERED/pwm")

BachaIndex(
  pwm = pwm,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  open_regions = ocr,
  mask_regions = mask_regions,
  tss_positions = promoters(TxDb.Hsapiens.UCSC.hg38.knownGene, upstream = 0, downstream = 1),
  match_features = c("GC", "entropy", "tssDist"),
  metadata = c(
    motif_source = "https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE-CLUSTERED/H14CORE-CLUSTERED_pwm.tar.gz (ReadPWM)",
    open_regions_source = "ENCODE cCREs v4, https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed",
    mask_regions_source = "AnnotationHub AH107305 (exclude ranges) + AH99003 (repeat regions)",
    tss_source = paste("TxDb.Hsapiens.UCSC.hg38.knownGene", packageVersion("TxDb.Hsapiens.UCSC.hg38.knownGene"), "transcript starts")
  ),
  output_path = "index/hg38.h14coreclust.v1.bi",
  motif_set_name = "H14CORE_CLUSTERED"
)
