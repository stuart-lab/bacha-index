library(GenomicRanges)
library(AnnotationHub)
library(BSgenome.Mmusculus.UCSC.mm39)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(TFBSTools)
library(bacha)
library(rtracklayer)
library(R.utils)
library(readxl)


source("fetch.R")
options(timeout = 10000)
set.seed(1234)

# load open chromatin regions
# the wenglab server drops large HTTP/2 transfers, so fetch over HTTP/1.1 with retries
fetch_file("https://downloads.wenglab.org/Registry-V4/mm10-cCREs.bed",
           destfile = "mm10-cCREs.bed")
ocr <- read.table("mm10-cCREs.bed", sep = "\t")
ocr <- makeGRangesFromDataFrame(ocr, seqnames.field = "V1", start.field = "V2", end.field = "V3")

# lift over to mm39
url <- "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz"
chain_file <- tempfile(fileext = ".chain.gz")
download.file(url, chain_file, mode = "wb")
gunzip(chain_file, remove = FALSE, overwrite = TRUE)
chain <- import.chain(sub("\\.gz$", "", chain_file))

ocr_mm39 <- unlist(liftOver(ocr, chain))
ocr_mm39 <- unique(ocr_mm39)
ocr_mm39 <- ocr_mm39[width(ocr_mm39) > 150]

# load regions to exclude
ah <- AnnotationHub()
excludeGR.mm39 <- ah[["AH107321"]]
repeat_regions <- ah[['AH99013']]
mask_regions <- c(excludeGR.mm39, repeat_regions)
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
  genome = BSgenome.Mmusculus.UCSC.mm39,
  open_regions = ocr_mm39,
  mask_regions = mask_regions,
  tss_positions = promoters(TxDb.Mmusculus.UCSC.mm39.knownGene, upstream = 0, downstream = 1),
  match_features = c("GC", "entropy", "tssDist"),
  metadata = c(
    motif_source = "Xie et al. 2025 composite motifs, https://doi.org/10.1038/s41586-025-08844-z (Supplementary Table, toPWM pseudocounts = 1)",
    open_regions_source = "ENCODE mm10 cCREs v4 (https://downloads.wenglab.org/Registry-V4/mm10-cCREs.bed) lifted over with mm10ToMm39.over.chain.gz, regions > 150 bp",
    mask_regions_source = "AnnotationHub AH107321 (exclude ranges) + AH99013 (repeat regions)",
    tss_source = paste("TxDb.Mmusculus.UCSC.mm39.knownGene", packageVersion("TxDb.Mmusculus.UCSC.mm39.knownGene"), "transcript starts")
  ),
  output_path = "index/mm39.Xie.v1.bi",
  motif_set_name = "Xie"
)
