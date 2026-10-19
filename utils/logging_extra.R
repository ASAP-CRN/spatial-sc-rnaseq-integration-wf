this_file <- function() {
  sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
}

# Tees calling script invocation (with a timestamp) and all subsequent stdout to
# `log_file`. R's message/stderr stream has no split-sink equivalent, so
# warnings/messages stay console-only, unlike the Python `log_run_command`.
log_run_command <- function(log_file) {
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  log_con <- file(log_file, open = "wt")
  cmd <- paste(c("Rscript", this_file(), commandArgs(trailingOnly = TRUE)), collapse = " ")
  cat(sprintf("[%s] %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), cmd), file = log_con)
  sink(log_con, split = TRUE)
}
