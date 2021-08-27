# creates example R matrices containing gene and gRNA single-cell data. Converts these matrices into ODM format.
library(ondisc)
library(magrittr)
library(Matrix)
xie_offsite <- paste0(.get_config_path("LOCAL_XIE_2019_DATA_DIR"), "processed/")
ondiscdata_dir <- "~/research_code/ondiscdata/"

gene_odm <- read_odm(odm_fp = paste0(xie_offsite, "gene/expression_matrix.odm"),
                     metadata_fp = paste0(xie_offsite, "gene/metadata.rds"))
gRNA_odm <- read_odm(odm_fp = paste0(xie_offsite, "gRNA/binary_grouped.odm"),
                     metadata_fp = paste0(xie_offsite, "gRNA/binary_grouped_metadata.rds"))

# extract first batch
set.seed(10)
barcodes_to_extract <- sort(sample(x = seq(1, ncol(gene_odm)), size = 2000, replace = FALSE))

gene_odm <- gene_odm[,barcodes_to_extract]
gRNA_odm <- gRNA_odm[,barcodes_to_extract]

# get gene barcodes and features
gene_mat <- gene_odm[[1:nrow(gene_odm),]]
gene_mat <- as(gene_mat, "dgTMatrix")
gene_feature_df <- data.frame(gene_id = get_feature_ids(gene_odm),
                              gene_name = get_feature_names(gene_odm))
gene_barcodes <- get_cell_barcodes(gene_odm)
saveRDS(object = gene_mat, paste0(ondiscdata_dir, "inst/extdata/r_matrix/gene/matrix.rds"))
saveRDS(object = gene_barcodes, paste0(ondiscdata_dir, "inst/extdata/r_matrix/gene/barcodes.rds"))
saveRDS(object = gene_feature_df, paste0(ondiscdata_dir, "inst/extdata/r_matrix/gene/features.rds"))

# get gRNA barcodes and features
gRNA_mat <- as.matrix(gRNA_odm[[1:nrow(gRNA_odm),]])
gRNA_barcodes <- get_cell_barcodes(gRNA_odm)
gRNA_feature_df <- data.frame(gRNA_id = get_feature_ids(gRNA_odm))
saveRDS(object = gRNA_mat, paste0(ondiscdata_dir, "inst/extdata/r_matrix/gRNA/matrix.rds"))
saveRDS(object = gRNA_barcodes, paste0(ondiscdata_dir, "inst/extdata/r_matrix/gRNA/barcodes.rds"))
saveRDS(object = gRNA_feature_df, paste0(ondiscdata_dir, "inst/extdata/r_matrix/gRNA/features.rds"))
