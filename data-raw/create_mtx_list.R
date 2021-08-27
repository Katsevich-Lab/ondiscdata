library(ondisc)
ondiscdata_dir <- "/Users/timbarry/research_code/ondiscdata/inst/extdata/mtx_list"

n_col <- c(2000, 3000, 2500)
cumsum_ncol <- c(0, cumsum(n_col))
n_row <- 10000
file_dirs <- paste0(ondiscdata_dir, paste0("/mtx_dir", seq(1, length(n_col))))

for (i in seq(1, length(n_col))) {
  if (!dir.exists(file_dirs[i])) dir.create(file_dirs[i])
  create_synthetic_data(n_row = n_row,
                        n_col = n_col[i],
                        logical_mat = FALSE,
                        write_as_mtx_to_disk = TRUE,
                        start_pos = cumsum_ncol[i],
                        write_as_h5_to_disk = FALSE,
                        file_dir = file_dirs[i])
  }
