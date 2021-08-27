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
