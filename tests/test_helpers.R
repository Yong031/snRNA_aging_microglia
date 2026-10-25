source(file.path("R", "utils.R"))

stopifnot(identical(
  classify_state(c(-1, 1, -1, 1), c(-1, -1, 1, 1)),
  c("Homeostatic-like", "Activated/pathological", "Primed-like", "Primed-pathological")
))

meta <- data.frame(
  sample_id = c("S1", "S2"),
  condition = c("Sham_Young", "SNI_Aged")
)
validate_metadata(meta, c("sample_id", "condition"))

failed <- FALSE
tryCatch(
  validate_metadata(meta, c("sample_id", "cell_type")),
  error = function(e) failed <<- grepl("cell_type", conditionMessage(e), fixed = TRUE)
)
stopifnot(failed)

tmp <- tempfile("snrna-output-")
out <- ensure_output_dir(tmp)
stopifnot(dir.exists(out), normalizePath(out) == normalizePath(tmp))

cat("helper tests passed\n")
