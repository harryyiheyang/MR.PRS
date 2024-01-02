combine_chrwise_profiles <- function(prs_file_dir, NAM, num_chr = 22) {
  prs_data_list <- list()

  for (chr in 1:num_chr) {
    chr_data_list <- list()
    for (i in 1:length(NAM)) {
      file_name <- glue("{prs_file_dir}/{NAM[i]}_{chr}.profile")
      if (file.exists(file_name)) {
        chr_data <- data.table::fread(file_name)
        setnames(chr_data, "SCORESUM", paste0(NAM[i], "_chr", chr))
        chr_data_list[[i]] <- chr_data
      }
    }

    if (length(chr_data_list) > 0) {
      combined_chr_data <- Reduce(function(x, y) merge(x, y, by = "ID"), chr_data_list)
      prs_data_list[[chr]] <- combined_chr_data
    }
  }
  combined_prs_data <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), prs_data_list)

  return(combined_prs_data)
}
