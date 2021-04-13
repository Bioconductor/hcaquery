#' @title Helper functions for the processing of .loom and .h5ad files into
#' relational database tables
#'
#' @rdname sce
#'
#' @param loom_filepath character() absolute path to downloaded .loom file
#'
#' @importFrom LoomExperiment import
#'
#' @return SingleCellLoomExperiment object
loom_to_sce <- function(loom_filepath) {
    loom_sce <- LoomExperiment::import(loom_filepath,
                                       type ="SingleCellLoomExperiment")
    loom_sce
}

#' @rdname sce
#'
#' @param h5ad_filepath character() absolute path to downloaded .h5ad file
#'
#' @importFrom zellkonverter readH5AD
#'
#' @return SingleCellExperiment object
h5ad_to_sce <- function(h5ad_filepath) {
    h5ad_sce <- zellkonverter::readH5AD(h5ad_filepath)
    h5ad_sce
}

#' @rdname sce
#'
#' @description  function to change SingleCellLoomExperiment object's assa
#'  matrix into sparse representation (`dgCMatrix`)
#'
#' @param loom_filepath character() absolute path to downloaded .loom file
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom methods as
#'
#' @return dgCMatrix sparse matrix representation
loom_sparse_matrix <- function(loom_filepath) {
    loom_sce_obj <- loom_to_sce(loom_filepath)
    sparse_loom_matrix <- as(SummarizedExperiment::assay(loom_sce_obj), "dgCMatrix")
    sparse_loom_matrix
}

#' @rdname sce
#'
#' @description function to change extract the existing sparse matrix
#' representation (`dgRMatrix`) of the SingleCellExperiment object's assay
#'
#' @param h5ad_filepath character() absolute path to downloaded .h5ad file
#'
#' @importFrom SummarizedExperiment assay
#'
#' @return dgRMatrix sparse matrix representation
h5ad_sparse_matrix <- function(h5ad_filepath) {
    h5ad_sce_obj <- h5ad_to_sce(h5ad_filepath)
    sparse_h5ad_matrix <- SummarizedExperiment::assay(h5ad_sce_obj)
    sparse_h5ad_matrix
}

#' @rdname sce
#'
#' @description convert a sparse matrix representation of the assay into a
#' reformatted tibble
#'
#' @param sparse_matrix dgCMatrix sparse matrix representation
#'
#' @importFrom tibble tibble
#'
#' @return tibble with row_index, col_index, and value of each non-zero entry
#' in the assay matrix
loom_sparse_mtx_to_assay_tbl <- function(sparse_matrix) {
    stopifnot(
        ## sparse_matrix must inherit from class dgCMatrix
        `'sparse_matrix =' must inherit from class 'dgCMatrix'` =
            inherits(sparse_matrix, "dgCMatrix")
    )

    ## handling of sparse matrices depends on the initial file type
    ## help: https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    ##dim(sparse_matrix)
    ##length(sparse_matrix@i)
    ##length(sparse_matrix@p)
    ##length(sparse_matrix@x)

    ## increment each element of row_index by one to match R's 1-based indexing
    row_index <- sparse_matrix@i + 1
    values <- sparse_matrix@x

    ## diff between each consecutive integer in the sequence
    ## sparse_matrix@p gives number of non-zero values in each column
    ## as the sequence is made by reading the matrix
    ## bottom to top, then left to right

    loom_n_per_column <-  diff(sparse_matrix@p)
    col_index <-  rep(seq_along(loom_n_per_column), loom_n_per_column)
    ##str(loom_col_index)

    # tibble is the row, column, and value for the entire matrix
    loom_assay_tbl <- tibble(row_index, col_index, values)
    loom_assay_tbl
}

#' @rdname sce
#'
#' @description convert a sparse matrix representation of the assay into a
#' reformatted tibble
#'
#' @param sparse_matrix dgRMatrix sparse matrix representation
#'
#' @importFrom tibble tibble
#'
#' @return tibble with row_index, col_index, and value of each non-zero entry
#' in the assay matrix
h5ad_sparse_mtx_to_assay_tbl <- function(sparse_matrix) {
    stopifnot(
        ## sparse_matrix must inherit from class dgRMatrix
        `'sparse_matrix =' must inherit from class 'dgRMatrix'` =
            inherits(sparse_matrix, "dgRMatrix")
    )

    ## handling of sparse matrices depends on the initial file type
    ## same as with .loom matrices, just compressed row wise vs column wise
    ##dim(sparse_matrix)
    ##length(sparse_matrix@j)
    ##length(sparse_matrix@p)
    ##length(sparse_matrix@x)

    ## increment each element of col_index by one to match R's 1-based indexing
    col_index <- sparse_matrix@j + 1
    values <- sparse_matrix@x

    ## diff between each consecutive integer in the sequence
    ## sparse_matrix@p gives number of non-zero values in each row

    h5ad_n_per_column <-  diff(sparse_matrix@p)
    row_index <-  rep(seq_along(h5ad_n_per_column), h5ad_n_per_column)
    ##str(h5ad_col_index)

    # tibble is the row, column, and value for the entire matrix
    h5ad_assay_tbl <- tibble(row_index, col_index, values)
    h5ad_assay_tbl
}

#' @rdname sce
#'
#' @description function to change SingleCellExperiment or
#' SingleCellLoomExperiment object's `rowData` into a gene annotation table
#'
#' @param sce SingleCellExperiment or SingleCellLoomExperiment object
#'
#' @importFrom dplyr as_tibble
#' @importFrom tibble add_column
#' @importFrom SummarizedExperiment rowData
#'
#' @return tibble of gene annotations
sce_rowdata_to_gene_tbl <- function(sce) {
    stopifnot(
        ## sce must be a SingleCellExperiment object
        `'sce =' must inherit from class 'SingleCellExperiment'` =
        inherits(sce, "SingleCellExperiment")
    )

    ## add row_index as a column
    gene_tbl <- dplyr::as_tibble(rowData(sce))
    gene_tbl %>% tibble::add_column(row_index = 1:nrow(gene_tbl),
                            .before = colnames(gene_tbl)[1])
    gene_tbl
}

#' @rdname sce
#'
#' @description function to change SingleCellExperiment or
#' SingleCellLoomExperiment object's `colData` into a cell annotation table
#'
#' @param sce SingleCellExperiment or SingleCellLoomExperiment object
#'
#' @importFrom dplyr as_tibble
#' @importFrom tibble add_column
#' @importFrom SummarizedExperiment colData
#'
#' @return tibble of cell annotations
sce_coldata_to_cell_tbl <- function(sce) {
    stopifnot(
        ## sce must be a SingleCellExperiment object
        `'sce =' must inherit from class 'SingleCellExperiment'` =
            inherits(sce, "SingleCellExperiment")
    )

    ## add col_index as a column
    cell_tbl <- dplyr::as_tibble(colData(sce))
    cell_tbl %>% tibble::add_column(col_index = 1:nrow(cell_tbl),
                                    .before = colnames(cell_tbl)[1])
    cell_tbl
}
