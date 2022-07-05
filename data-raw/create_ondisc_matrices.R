ondiscdata_dir <- "/Users/timbarry/research_code/ondiscdata/"
library(ondisc)

# Load the gene data and create an odm
gene_mat <- readRDS(paste0(ondiscdata_dir, "inst/extdata/r_matrix/gene/matrix.rds"))
gene_barcodes <- readRDS(paste0(ondiscdata_dir, "inst/extdata/r_matrix/gene/barcodes.rds"))
gene_feature_df <- readRDS(paste0(ondiscdata_dir, "inst/extdata/r_matrix/gene/features.rds"))
to_save_odm <- paste0(ondiscdata_dir, "inst/extdata/odm/gene/matrix.odm")
to_save_metadata <- paste0(ondiscdata_dir, "inst/extdata/odm/gene/metadata.rds")
odm <- create_ondisc_matrix_from_R_matrix(r_matrix = gene_mat, barcodes = gene_barcodes,
                                          features_df = gene_feature_df, odm_fp = to_save_odm,
                                          metadata_fp = to_save_metadata)

# load the gRNA data and create an odm
gene_mat <- readRDS(paste0(ondiscdata_dir, "inst/extdata/r_matrix/gRNA/matrix.rds"))
gene_barcodes <- readRDS(paste0(ondiscdata_dir, "inst/extdata/r_matrix/gRNA/barcodes.rds"))
gene_feature_df <- readRDS(paste0(ondiscdata_dir, "inst/extdata/r_matrix/gRNA/features.rds"))
to_save_odm <- paste0(ondiscdata_dir, "inst/extdata/odm/gRNA/matrix.odm")
to_save_metadata <- paste0(ondiscdata_dir, "inst/extdata/odm/gRNA/metadata.rds")
odm <- create_ondisc_matrix_from_R_matrix(r_matrix = gene_mat, barcodes = gene_barcodes,
                                          features_df = gene_feature_df, odm_fp = to_save_odm,
                                          metadata_fp = to_save_metadata)

# synthesize gRNA count data
set.seed(4)
n_cells <- 10000
n_gRNAs <- 100
gRNA_exp <- replicate(n = n_cells, expr = {
  v <- integer(n_gRNAs)
  rand_idxs <- sample(x = seq(1, n_gRNAs), size = sample(x = seq(1,4), size = 1), replace = FALSE)
  v[rand_idxs] <- rpois(n = length(rand_idxs), lambda = 100)
  v
}, simplify = TRUE)
gRNA_ids <- paste0("gRNA_", seq(1, n_gRNAs))
cell_ids <- paste0("cell_", seq(1, n_cells))
gRNA_groups <- rep(paste0("gRNA_group_", seq(1, n_gRNAs/2)), each = 2)
to_save_odm <- paste0(ondiscdata_dir, "inst/extdata/odm/gRNA_expression/matrix.odm")
to_save_metadata <- paste0(ondiscdata_dir, "inst/extdata/odm/gRNA_expression/metadata.rds")

odm <- create_ondisc_matrix_from_R_matrix(r_matrix = gRNA_exp,
                                          barcodes = cell_ids,
                                          features_df = data.frame(gRNA_ids),
                                          odm_fp = to_save_odm,
                                          metadata_fp = to_save_metadata)
odm <- odm |> ondisc::mutate_feature_covariates(gRNA_group = gRNA_groups)
save_odm(odm = odm, metadata_fp = to_save_metadata)
