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
mem.maxVSize(vsize = Inf)

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
fetch_file("https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/Vierstra-Human-Motifs.rds",
           destfile = "Vierstra-Human-Motifs.rds")

pwm <- readRDS("Vierstra-Human-Motifs.rds")

# The matrices carry empty ID and name slots; the identifiers are only in the
# list names (e.g. "AHR_1#AHR:bHLH:284:1"). Store the full identifier as the
# feature ID and the part before "#" (e.g. "AHR_1") as its display name.
pwm <- do.call(PWMatrixList, lapply(names(pwm), function(id) {
  m <- pwm[[id]]
  m@ID <- id
  m@name <- sub("#.*$", "", id)
  m
}))

BachaIndex(
  pwm = pwm,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  open_regions = ocr,
  mask_regions = mask_regions,
  tss_positions = promoters(TxDb.Hsapiens.UCSC.hg38.knownGene, upstream = 0, downstream = 1),
  match_features = c("GC", "entropy", "tssDist"),
  metadata = c(
    motif_source = "https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/Vierstra-Human-Motifs.rds",
    open_regions_source = "ENCODE cCREs v4, https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed",
    mask_regions_source = "AnnotationHub AH107305 (exclude ranges) + AH99003 (repeat regions)",
    tss_source = paste("TxDb.Hsapiens.UCSC.hg38.knownGene", packageVersion("TxDb.Hsapiens.UCSC.hg38.knownGene"), "transcript starts")
  ),
  output_path = "index/hg38.Vierstra.v1.bi",
  motif_set_name = "Vierstra"
)
