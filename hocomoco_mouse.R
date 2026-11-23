# Restrict a HOCOMOCO v14 motif set to TFs that exist in mouse.
#
# H14CORE is TF-centric and spans human and mouse: each entry's masterlist
# records which species the TF is found in. 298 of the 1107 core TFs are
# human-only, and there are no mouse-only TFs -- the masterlist is anchored on
# the human complement -- so this filter is only ever needed in one direction.
#
# Sourced by the mm10/mm39 HOCOMOCO index builds.

hocomoco_mouse_tfs <- function(
    annotation = "H14CORE_annotation.jsonl",
    url = paste0("https://hocomoco14.autosome.org/final_bundle/hocomoco14/",
                 "H14CORE/H14CORE_annotation.jsonl")) {
  fetch_file(url, annotation)
  ann <- lapply(readLines(annotation), jsonlite::fromJSON)
  tf <- vapply(
    ann,
    function(a) {
      if ("MOUSE" %in% names(a$masterlist_info$species)) {
        a$masterlist_info$tf
      } else {
        NA_character_
      }
    },
    character(1)
  )
  unique(tf[!is.na(tf)])
}

# Subset a PWM list to those motifs.
#
# The two readers name motifs differently, so match on the tf_id, which is
# everything before the first dot in either convention:
#
#   ReadJASPAR  AHR.H14CORE.0.P.B  (full HOCOMOCO id)  -> AHR
#   ReadPWM     ANDR, ANDR.1       (numeric suffix)    -> ANDR
#
# tf_ids never contain a dot, so stripping from the first one is unambiguous.
subset_mouse_motifs <- function(pwm, mouse_tfs = hocomoco_mouse_tfs(),
                                verbose = TRUE) {
  base <- sub("\\..*$", "", names(pwm))
  keep <- base %in% mouse_tfs
  if (verbose) {
    message(sprintf(
      "mouse subset: keeping %d of %d motifs (%d human-only TFs dropped)",
      sum(keep), length(keep), sum(!keep)
    ))
  }
  pwm[keep]
}
