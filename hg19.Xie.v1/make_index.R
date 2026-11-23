library(GenomicRanges)
library(AnnotationHub)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TFBSTools)
library(bacha)
library(rtracklayer)
library(R.utils)
library(readxl)

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
fetch_file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-025-08844-z/MediaObjects/41586_2025_8844_MOESM5_ESM.xlsx",
           destfile = "composite_motifs.xlxs")

read_pwms_from_excel <- function(path, sheet = 1) {
  dat <- read_excel(path, sheet = sheet, col_names = FALSE, skip = 13, col_types = "text")
  lines <- trimws(unlist(dat))
  header_idx <- grep("^Base", lines)
  motifs <- list()
  
  for (i in seq_along(header_idx)) {
    start <- header_idx[i]
    next_header <- if (i < length(header_idx)) header_idx[i + 1] else length(lines)
    name <- dat[start, 2][[1]]
    mat <- as.matrix(dat[(start + 1):(start + 4), 2:ncol(dat)])
    # trim NA positions
    na_pos <- which(is.na(mat[1, ]))
    if (length(x = na_pos) > 0) {
      mat <- mat[, 1:na_pos[[1]]-1]
    }
    mat <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat))
    rownames(mat) <- c("A", "C", "G", "T")
    
    pfm <- PFMatrix(ID = name, name = name, profileMatrix = mat)
    pwm <- suppressWarnings(toPWM(pfm, pseudocounts = 1))
    motifs[[name]] <- pwm
  }
  
  do.call(PWMatrixList, motifs)
}

pwm <- read_pwms_from_excel("composite_motifs.xlxs")

# create Bacha index
BachaIndex(
  pwm = pwm,
  genome = BSgenome.Hsapiens.UCSC.hg19,
  open_regions = ocr_hg19,
  mask_regions = mask_regions,
  tss_positions = promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream = 0, downstream = 1),
  match_features = c("GC", "entropy", "tssDist"),
  metadata = c(
    motif_source = "Xie et al. 2025 composite motifs, https://doi.org/10.1038/s41586-025-08844-z (Supplementary Table, toPWM pseudocounts = 1)",
    open_regions_source = "ENCODE hg38 cCREs v4 (https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed) lifted over with hg38ToHg19.over.chain.gz, regions > 150 bp",
    mask_regions_source = "AnnotationHub AH107314 (exclude ranges) + AH99002 (repeat regions)",
    tss_source = paste("TxDb.Hsapiens.UCSC.hg19.knownGene", packageVersion("TxDb.Hsapiens.UCSC.hg19.knownGene"), "transcript starts")
  ),
  output_path = "index/hg19.Xie.v1.bi",
  motif_set_name = "Xie"
)
