# Shared download helper for the make_index.R scripts.
#
# Several index builds run concurrently from the repo root and share input
# files (cCRE beds, motif files). A plain `if (!file.exists(f)) download.file()`
# lets one job read a file another job is still writing. fetch_file() avoids
# that by:
#   * checking an existing file's size against the server's Content-Length,
#     and re-downloading if it differs;
#   * downloading to a job-private temp file in the same directory, verifying
#     its size, then renaming it into place (atomic on the same filesystem),
#     so readers only ever see a missing file or a complete one.

remote_size <- function(url) {
  hdr <- tryCatch(
    suppressWarnings(system2("curl", c("-sIL", "--http1.1", "--max-time", "60",
                                       shQuote(url)),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) character()
  )
  cl <- grep("^content-length:", hdr, ignore.case = TRUE, value = TRUE)
  if (length(cl) == 0) return(NA_real_)
  # last one wins when redirects are followed
  as.numeric(trimws(sub("^[^:]+:", "", tail(cl, 1))))
}

fetch_file <- function(url, destfile) {
  expected <- remote_size(url)

  if (file.exists(destfile)) {
    if (is.na(expected)) {
      message("fetch_file: could not get remote size for ", basename(destfile),
              "; using existing file unverified")
      return(invisible(destfile))
    }
    if (file.size(destfile) == expected) {
      return(invisible(destfile))
    }
    message("fetch_file: ", destfile, " is ", file.size(destfile),
            " bytes, expected ", expected, "; re-downloading")
  }

  tmp <- tempfile(pattern = paste0(".", basename(destfile), ".part."),
                  tmpdir = dirname(destfile))
  on.exit(unlink(tmp), add = TRUE)
  status <- download.file(url, destfile = tmp, method = "curl", quiet = TRUE,
                          extra = "--http1.1 --retry 5 --retry-all-errors --fail -L")
  if (status != 0) stop("fetch_file: download failed for ", url)
  got <- file.size(tmp)
  if (!is.na(expected) && got != expected) {
    stop("fetch_file: ", basename(destfile), " downloaded ", got,
         " bytes, expected ", expected)
  }
  if (!file.rename(tmp, destfile)) {
    stop("fetch_file: could not move download into place at ", destfile)
  }
  invisible(destfile)
}
