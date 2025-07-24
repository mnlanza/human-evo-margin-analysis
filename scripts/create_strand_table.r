source("scripts/utils.r")

create_strand_table <- function(ifn_tsv, ofn) {
  cat(sprintf("Reading TSV file %s\n", ifn_tsv))
  df <- read.delim(ifn_tsv)

  # Extract sequence IDs by removing _P/_M suffix
  base_ids <- unique(sub("_[PM]$", "", df$seq_id))
  rr <- data.frame(
    id = base_ids,
    plus = df$total_log_likelihood[match(paste0(base_ids, "_P"), df$seq_id)],
    minus = df$total_log_likelihood[match(paste0(base_ids, "_M"), df$seq_id)]
  )

  cat(sprintf("Saving to %s\n", ofn))
  write.table(rr, ofn, sep = "\t", row.names = FALSE, quote = FALSE)
}
