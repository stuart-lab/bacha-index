library(GenomicRanges)
library(AnnotationHub)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TFBSTools)
library(bacha)
library(rtracklayer)
library(R.utils)
library(Signac)

source("fetch.R")
options(timeout = 10000)
set.seed(1234)

# load open chromatin regions
fetch_file("https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed",
           destfile = "GRCh38-cCREs.bed")
ocr <- read.table("GRCh38-cCREs.bed", sep = "\t")
ocr <- makeGRangesFromDataFrame(ocr, seqnames.field = "V1", start.field = "V2", end.field = "V3")

# lift over to hg19 (ENCODE cCREs v4 are only released for GRCh38)
fetch_file("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
           destfile = "hg38ToHg19.over.chain.gz")
chain_file <- tempfile(fileext = ".chain")
gunzip("hg38ToHg19.over.chain.gz", destname = chain_file, remove = FALSE, overwrite = TRUE)
chain <- import.chain(chain_file)

ocr_hg19 <- unlist(liftOver(ocr, chain))
ocr_hg19 <- unique(ocr_hg19)
ocr_hg19 <- ocr_hg19[width(ocr_hg19) > 150]

# load regions to exclude
ah <- AnnotationHub()
excludeGR.hg19 <- ah[["AH107314"]]
repeat_regions <- ah[['AH99002']]
mask_regions <- c(excludeGR.hg19, repeat_regions)
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
  genome = BSgenome.Hsapiens.UCSC.hg19,
  open_regions = ocr_hg19,
  mask_regions = mask_regions,
  tss_positions = promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream = 0, downstream = 1),
  match_features = c("GC", "entropy", "tssDist"),
  metadata = c(
    motif_source = "https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE-CLUSTERED/H14CORE-CLUSTERED_pwm.tar.gz (ReadPWM)",
    open_regions_source = "ENCODE hg38 cCREs v4 (https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed) lifted over with hg38ToHg19.over.chain.gz, regions > 150 bp",
    mask_regions_source = "AnnotationHub AH107314 (exclude ranges) + AH99002 (repeat regions)",
    tss_source = paste("TxDb.Hsapiens.UCSC.hg19.knownGene", packageVersion("TxDb.Hsapiens.UCSC.hg19.knownGene"), "transcript starts")
  ),
  output_path = "index/hg19.h14coreclust.v1.bi",
  motif_set_name = "H14CORE_CLUSTERED"
)
