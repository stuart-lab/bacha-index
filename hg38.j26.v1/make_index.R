library(GenomicRanges)
library(AnnotationHub)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Signac)
library(TFBSTools)
library(bacha)

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

fetch_file("https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
           destfile = "JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt")
pwm <- ReadJASPAR("JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
                  mode = "jaspar", keep = "all")

# BachaIndex() takes each feature's ID (the JASPAR accession) and display
# name (the TF) from the matrices themselves; TFs with several matrices
# are labelled "NAME (ID)" in the output.

# create Bacha index
BachaIndex(
  pwm = pwm,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  open_regions = ocr,
  mask_regions = mask_regions,
  tss_positions = promoters(TxDb.Hsapiens.UCSC.hg38.knownGene, upstream = 0, downstream = 1),
  match_features = c("GC", "entropy", "tssDist"),
  metadata = c(
    motif_source = "https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
    open_regions_source = "ENCODE cCREs v4, https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed",
    mask_regions_source = "AnnotationHub AH107305 (exclude ranges) + AH99003 (repeat regions)",
    tss_source = paste("TxDb.Hsapiens.UCSC.hg38.knownGene", packageVersion("TxDb.Hsapiens.UCSC.hg38.knownGene"), "transcript starts")
  ),
  output_path = "index/hg38.j26.v1.bi",
  motif_set_name = "JASPAR2026"
)
